"""Output writers: chain tables and publication-ready figures.

Deliberately NOT a narrative report: the tool emits machine-readable tables plus
standalone figures that drop straight into a manuscript or downstream analysis.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .attribution import super_lineage_of  # noqa: E402

CHAINS_FIELDS = ["gene", "trait_id", "locus_id", "lead_snp", "cell_state", "lineage",
                 "super_lineage", "evidence_tier", "n_refs_covering", "n_agree",
                 "tau_mean", "per_pack_detail"]


def write_chains(results: list[dict], out_tsv: Path) -> None:
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=CHAINS_FIELDS)
        w.writeheader()
        for r in results:
            detail = "; ".join(f"{k}:{v['cell_state']}(tau={v['tau']})"
                               for k, v in r["per_pack"].items()) or ""
            w.writerow({
                "gene": r["gene"],
                "trait_id": r.get("trait_id", ""),
                "locus_id": r.get("locus_id", ""),
                "lead_snp": r.get("lead_snp", ""),
                "cell_state": r["consensus_state"] or "",
                "lineage": r["consensus_lineage"] or "",
                "super_lineage": r.get("consensus_super_lineage") or "",
                "evidence_tier": r["evidence_tier"],
                "n_refs_covering": r["n_refs_covering"],
                "n_agree": r["n_agree"],
                "tau_mean": r["tau_mean"] if r["tau_mean"] is not None else "",
                "per_pack_detail": detail,
            })


def plot_gene_celltype(r: dict, out_png: Path, out_pdf: Path | None = None) -> None:
    """Per-gene cell-state attribution across reference packs."""
    packs = list(r["per_pack"].keys())
    if not packs:
        return
    states = [r["per_pack"][p]["cell_state"] for p in packs]
    taus = [r["per_pack"][p]["tau"] or 0.0 for p in packs]
    consensus_super = r.get("consensus_super_lineage")
    colors = ["#2b6cb0" if super_lineage_of(s) == consensus_super else "#a0aec0"
              for s in states]
    fig, ax = plt.subplots(figsize=(max(4, 1.4 * len(packs) + 2), 3.2), dpi=150)
    bars = ax.bar(range(len(packs)), taus, color=colors)
    ax.set_xticks(range(len(packs)))
    ax.set_xticklabels([f"{p}\n{s}" for p, s in zip(packs, states)],
                       fontsize=7, rotation=20, ha="right")
    ax.set_ylabel("tau specificity")
    ax.set_ylim(0, 1.05)
    ax.set_title(f"{r['gene']}  |  consensus: {r['consensus_state']} "
                 f"({r['evidence_tier']})", fontsize=9)
    for b, t in zip(bars, taus):
        ax.text(b.get_x() + b.get_width() / 2, t + 0.02, f"{t:.2f}",
                ha="center", fontsize=6)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png)
    if out_pdf:
        fig.savefig(out_pdf)
    plt.close(fig)


def plot_evidence_tiers(results: list[dict], out_png: Path) -> None:
    from collections import Counter
    c = Counter(r["evidence_tier"] for r in results)
    order = ["cell_type_confirmed", "cell_type_consensus_partial", "cell_type_single",
             "cell_type_discordant", "gene_not_in_reference"]
    labels = [k for k in order if k in c]
    vals = [c[k] for k in labels]
    fig, ax = plt.subplots(figsize=(5.5, 3), dpi=150)
    ax.barh(range(len(labels)), vals, color="#2b6cb0")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("genes")
    ax.set_title("cell-type attribution evidence tiers", fontsize=9)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png)
    plt.close(fig)


def write_manifest(d: dict, out_json: Path) -> None:
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(d, indent=2))
