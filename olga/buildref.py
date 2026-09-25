"""Build a cell-type reference pack from any annotated h5ad atlas.

A reference pack is the DERIVED statistics the attribution needs, so raw
single-cell atlases (GBs, redistribution-restricted) never ship with the tool:

  tau_matrix.tsv       tau specificity + argmax state (continuous attribution)
  markers_strict.tsv   mutually exclusive marker assignment (one state per gene)
  expr_by_cluster.tsv  mean normalised expression, gene x cell state
  manifest.json        provenance, sizes, gene space, method, licence note

Handles the two conventions seen in public atlases: Ensembl var_names with a
``feature_name``/``gene_symbols`` var column, and auto-detected cluster columns.
``anndata`` is imported lazily so the base install stays light.
"""
from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

import numpy as np

MIN_CELLS = 10
MIN_LOG2FC = 0.5
MIN_MARGIN = 0.25
CLUSTER_CANDIDATES = ["cell_type", "Cluster", "celltype", "cluster",
                      "cell_type_group", "cell_name_detailed", "leiden", "louvain"]


def resolve_symbols(adata):
    var = adata.var
    for col in ("feature_name", "gene_symbols"):
        if col in var.columns:
            syms = [str(x) for x in var[col]]
            if sum(s.startswith("ENSG") for s in syms) < 0.5 * len(syms):
                return syms, col
    return [str(x) for x in adata.var_names], "var_names"


def pick_cluster_col(adata, requested):
    if requested:
        return requested
    for c in CLUSTER_CANDIDATES:
        if c in adata.obs.columns:
            return c
    raise SystemExit(f"no cluster column found; tried {CLUSTER_CANDIDATES}")


def build_reference(atlases: list[tuple[str, str | None]], name: str, out_path: str,
                    cluster_col: str | None = None) -> dict:
    """atlases: list of (h5ad_path, state_label|None). Multiple atlases are
    merged on the gene UNION; genes absent from an atlas contribute zero
    expression there, and its states take the given prefix (e.g. "Epi") —
    the same convention that produced the bundled uc_gut_smillie pack.
    """
    import anndata as ad

    out = Path(out_path)
    out.mkdir(parents=True, exist_ok=True)

    gene_index: dict[str, int] = {}
    state_names: list[str] = []
    means = []  # (genes list, M matrix) per atlas
    n_cells = 0
    sym_source = None
    ccol_used = cluster_col
    for path, label in atlases:
        adata = ad.read_h5ad(path)
        n_cells += int(adata.n_obs)
        syms, sym_source = resolve_symbols(adata)
        ccol = pick_cluster_col(adata, ccol_used)
        ccol_used = ccol_used or ccol
        for g in syms:
            gene_index.setdefault(g, len(gene_index))
        labels = [str(x) for x in adata.obs[ccol]]
        counts = Counter(labels)
        keep = [c for c in sorted({s for s in labels})
                if counts[c] >= MIN_CELLS]
        X = adata.X.tocsr().copy()
        X.data = np.expm1(X.data)
        states = [f"{label}|{c}" if label else c for c in keep]
        state_names.extend(states)
        M = np.zeros((len(syms), len(states)), dtype=np.float64)
        lab = np.array(labels)
        for j, c in enumerate(keep):
            mask = lab == c
            if mask.any():
                M[:, j] = np.asarray(X[mask].mean(axis=0)).ravel()
        means.append((syms, states, M))
        del adata, X

    G, C = len(gene_index), len(state_names)
    if C == 0:
        raise SystemExit(
            f"no cell state reached MIN_CELLS={MIN_CELLS} in any input atlas; "
            "nothing to build a reference from")
    E = np.zeros((G, C), dtype=np.float64)
    col = 0
    for syms, states, M in means:
        rows = [gene_index[g] for g in syms]
        E[np.ix_(rows, np.arange(col, col + M.shape[1]))] = M
        col += M.shape[1]
    genes = [None] * G
    for g, i in gene_index.items():
        genes[i] = g

    tau = np.full(G, np.nan)
    valid = E.max(axis=1) > 0
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(E.max(axis=1, keepdims=True) > 0,
                         E / np.where(E.max(axis=1, keepdims=True) > 0,
                                      E.max(axis=1, keepdims=True), 1.0), 0.0)
        if C > 1:
            tau[valid] = ((1.0 - ratio[valid]).sum(axis=1)) / (C - 1)
    top = np.argmax(E, axis=1)

    strict = []
    order = np.argsort(-E, axis=1)
    for i, g in enumerate(genes):
        best = order[i, 0]
        second = order[i, 1] if C > 1 else best   # single-state pack: no contrast
        val, val2 = E[i, best], E[i, second]
        if val < MIN_LOG2FC or (val - val2) < MIN_MARGIN:
            continue
        strict.append({"gene": g, "cell_state": state_names[best],
                       "score": f"{val:.4f}", "margin": f"{val - val2:.4f}"})

    with (out / "tau_matrix.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene", "tau", "top_cell_state", "top_expr", "n_states"])
        for i, g in enumerate(genes):
            w.writerow([g, f"{tau[i]:.4f}" if np.isfinite(tau[i]) else "",
                        state_names[top[i]], f"{E[i, top[i]]:.4f}", C])
    with (out / "expr_by_cluster.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene"] + state_names)
        for i, g in enumerate(genes):
            w.writerow([g] + [f"{v:.4f}" for v in E[i]])
    with (out / "markers_strict.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t",
                           fieldnames=["gene", "cell_state", "score", "margin"])
        w.writeheader()
        w.writerows(sorted(strict, key=lambda r: (r["cell_state"], r["gene"])))
    manifest = {
        "name": name,
        "source": "; ".join(f"{p}" for p, _ in atlases),
        "n_cells": n_cells, "n_genes": G,
        "n_cell_states": C, "cluster_column": ccol_used,
        "symbol_source": sym_source, "min_cells": MIN_CELLS,
        "cell_states": state_names,
        "method": {"tau": "sum_i(1-x_i/x_max)/(n-1) on mean expm1(log-normalised X), "
                          "multi-atlas gene-union merge with zero fill",
                   "markers_strict": f"pseudobulk argmax proxy, expr>={MIN_LOG2FC}, "
                                     f"margin>={MIN_MARGIN} (bundled packs used one-vs-rest "
                                     f"Wilcoxon markers; tau is identical)"},
        "licence_note": "derived statistics only; cite the original atlas",
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest
