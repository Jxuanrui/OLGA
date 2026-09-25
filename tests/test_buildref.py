"""build-reference edge cases.

Regression locks for two crashes found in the pre-manuscript audit:
a pack reduced to ONE cell state (IndexError) and a pack with ZERO states
(opaque numpy zero-size error).
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from olga import buildref  # noqa: E402


def _write_h5ad(path, X, labels, genes):
    import anndata as ad
    import scipy.sparse as sp
    adata = ad.AnnData(X=sp.csr_matrix(np.log1p(X)), obs={"cell_type": labels})
    adata.var_names = genes
    adata.write_h5ad(path)
    return path


def test_tau_all_zero_gene_writes_empty(tmp_path):
    """A gene with zero expression everywhere -> empty tau, not a crash."""
    X = np.vstack([np.tile([5.0, 0.0], (10, 1)),   # 10 A cells: G1 only
                   np.zeros((10, 2))])             # 10 B cells: silent
    _write_h5ad(tmp_path / "tiny.h5ad", X, ["A"] * 10 + ["B"] * 10, ["G1", "G2"])
    m = buildref.build_reference([(str(tmp_path / "tiny.h5ad"), None)],
                                 "tiny", str(tmp_path / "pack"))
    tau = {}
    with (tmp_path / "pack" / "tau_matrix.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tau[r["gene"]] = r["tau"]
    assert float(tau["G1"]) == pytest.approx(1.0)   # expressed in one state only
    assert tau["G2"] == ""                          # silent gene -> empty
    assert m["n_cell_states"] == 2 and m["n_genes"] == 2


def test_build_reference_single_state_pack(tmp_path):
    """Regression: one state under MIN_CELLS leaves a single-state pack,
    which used to crash on the second-best-state index."""
    X = np.vstack([np.tile([1.0], (14, 1)), np.tile([2.0], (2, 1))])
    _write_h5ad(tmp_path / "t.h5ad", X, ["A"] * 14 + ["B"] * 2, ["G1"])
    m = buildref.build_reference([(str(tmp_path / "t.h5ad"), None)],
                                 "t", str(tmp_path / "p"))
    assert m["n_cell_states"] == 1                  # B dropped (< MIN_CELLS)


def test_build_reference_zero_states_clear_error(tmp_path):
    """Regression: every state under MIN_CELLS used to die inside numpy;
    it must stop with a clear message instead."""
    X = np.ones((4, 1))
    _write_h5ad(tmp_path / "t.h5ad", X, ["A", "A", "B", "B"], ["G1"])
    with pytest.raises(SystemExit, match="MIN_CELLS"):
        buildref.build_reference([(str(tmp_path / "t.h5ad"), None)],
                                 "t", str(tmp_path / "p"))


def test_multi_atlas_merge_prefixes_states(tmp_path):
    X = np.tile([3.0], (12, 1))
    p1 = _write_h5ad(tmp_path / "a.h5ad", X, ["c1"] * 12, ["G1"])
    p2 = _write_h5ad(tmp_path / "b.h5ad", X, ["c2"] * 12, ["G1"])
    m = buildref.build_reference([(str(p1), "Epi"), (str(p2), "LP")],
                                 "two", str(tmp_path / "pack"))
    assert m["cell_states"] == ["Epi|c1", "LP|c2"]
