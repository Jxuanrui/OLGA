"""LD-matrix / summary-statistic SNP-order alignment (protocol level).

Historical bug class this locks: an LD matrix written in reference-BIM order
consumed by code assuming sumstat order. Production consumers realign BY NAME
(`R <- R[rs, rs]` in the coloc.susie / ColocBoost runners). This test encodes
that protocol: given files whose SNP orders disagree, name-based alignment
must recover the true correlation matrix exactly, and any order-loss must be
caught by the set-equality guard rather than silently misaligned.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def _toeplitz_corr(n, rho):
    idx = np.arange(n)
    return rho ** np.abs(idx[:, None] - idx[None, :])


def _read_tsv_matrix(path):
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            rows.append([float(x) for x in line.rstrip("\n").split("\t")[1:]])
    return header[1:], np.array(rows)


def _read_snp_column(path):
    with open(path) as fh:
        fh.readline()
        return [line.split("\t", 1)[0].strip() for line in fh if line.strip()]


def test_name_alignment_recovers_true_ld(tmp_path):
    rng = np.random.default_rng(7)
    n = 40
    rsids = sorted(f"rs{i}" for i in range(1000, 1000 + n))
    C = _toeplitz_corr(n, 0.9)

    gwas_order = list(rng.permutation(rsids))      # sumstats order
    ld_order = list(rng.permutation(rsids))        # LD file order (e.g. BIM)

    idx = {r: i for i, r in enumerate(rsids)}
    with (tmp_path / "case.ld.tsv").open("w") as fh:
        fh.write("\t".join(["snp"] + ld_order) + "\n")
        for r in ld_order:
            fh.write("\t".join([r] + [f"{C[idx[r], idx[c]]:.6f}" for c in ld_order]) + "\n")
    with (tmp_path / "case.gwas.tsv").open("w") as fh:
        fh.write("snp\tbeta\tse\tN\tmaf\n")
        for r in gwas_order:
            fh.write(f"{r}\t0.10\t0.10\t1000\t0.30\n")

    ld_cols, R = _read_tsv_matrix(tmp_path / "case.ld.tsv")
    ld_rows = _read_snp_column(tmp_path / "case.ld.tsv")
    g = _read_snp_column(tmp_path / "case.gwas.tsv")

    # consumer contract: same variant set, row names == column names
    assert set(ld_cols) == set(ld_rows) == set(g)
    assert ld_cols != g or ld_rows == g            # orders may differ; sets must not

    pos = {r: i for i, r in enumerate(ld_cols)}
    aligned = R[[pos[r] for r in g]][:, [pos[r] for r in g]]   # R[rs, rs]
    truth = C[[idx[r] for r in g]][:, [idx[r] for r in g]]
    assert np.allclose(aligned, truth, atol=1e-6)


def test_set_mismatch_is_detected_not_misaligned(tmp_path):
    """A sumstat file with a SNP absent from the LD file must be caught by
    the set guard (the R runners' stopifnot equivalent), never silently
    dropped or misaligned."""
    ld_rs = [f"rs{i}" for i in range(30)]
    gwas_rs = ld_rs[:-1] + ["rsABSENT"]
    with (tmp_path / "m.ld.tsv").open("w") as fh:
        fh.write("\t".join(["snp"] + ld_rs) + "\n")
        for r in ld_rs:
            fh.write("\t".join([r] + ["0.0"] * 30) + "\n")
    with (tmp_path / "m.gwas.tsv").open("w") as fh:
        fh.write("snp\tbeta\tse\tN\tmaf\n")
        for r in gwas_rs:
            fh.write(f"{r}\t0.1\t0.1\t1000\t0.3\n")

    ld_cols, _ = _read_tsv_matrix(tmp_path / "m.ld.tsv")
    g = _read_snp_column(tmp_path / "m.gwas.tsv")
    assert set(g) - set(ld_cols) == {"rsABSENT"}   # guard fires here
