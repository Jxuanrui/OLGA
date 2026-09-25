"""CLI input parsing and output schema tests."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from olga.cli import main, read_effector_records  # noqa: E402
from olga.outputs import CHAINS_FIELDS  # noqa: E402


def test_read_effector_records_modes(tmp_path):
    eff = tmp_path / "eff.tsv"
    eff.write_text("trait_id\tlocus_id\tlead_snp\teffector_gene\n"
                   "T1\tchr1_x\trs1\tFUT2\n"
                   "T1\tchr2_y\trs2\t\n"          # empty gene -> skipped
                   "T1\tchr3_z\trs3\tNA\n")        # NA -> skipped
    recs, skipped = read_effector_records(str(eff))
    assert skipped == 2 and len(recs) == 1
    assert recs[0][0] == "FUT2"
    assert recs[0][1] == {"trait_id": "T1", "locus_id": "chr1_x", "lead_snp": "rs1"}

    plain = tmp_path / "plain.txt"
    plain.write_text("# comment\nFUT2\nMCM6\n\nLCT\n")
    recs2, skipped2 = read_effector_records(str(plain))
    assert [g for g, _ in recs2] == ["FUT2", "MCM6", "LCT"] and skipped2 == 0

    trait_only = tmp_path / "trait.tsv"
    trait_only.write_text("trait_id\tfoo\nT1\t1\n")
    with pytest.raises(SystemExit):
        read_effector_records(str(trait_only))


def test_chains_schema_constant():
    assert CHAINS_FIELDS == ["gene", "trait_id", "locus_id", "lead_snp",
                             "cell_state", "lineage", "super_lineage",
                             "evidence_tier", "n_refs_covering", "n_agree",
                             "tau_mean", "per_pack_detail"]


def test_run_without_genes_is_a_clear_error(capsys):
    assert main(["run", "--genes", "", "--out", "/tmp/olga_no_genes"]) == 2
    assert "provide effector genes" in capsys.readouterr().err
