"""Lightweight end-to-end smoke test on the bundled packs (CI-safe)."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from olga.cli import main  # noqa: E402


def test_run_readme_example(tmp_path):
    out = Path(tmp_path) / "results"
    assert main(["run", "--genes", "FUT2,MCM6,LCT", "--out", str(out)]) == 0
    rows = {r["gene"]: r for r in
            csv.DictReader((out / "chains.tsv").open(), delimiter="\t")}
    assert rows["FUT2"]["evidence_tier"] == "cell_type_confirmed"
    assert rows["FUT2"]["super_lineage"] == "epithelial"
    assert rows["MCM6"]["super_lineage"] == "proliferative"
    assert (out / "run_manifest.json").exists()
    assert (out / "celltypes" / "FUT2.celltype.png").exists()


def test_run_whole_chain_table_carries_context(tmp_path):
    eff = Path(tmp_path) / "eff.tsv"
    eff.write_text("trait_id\tlocus_id\tlead_snp\teffector_gene\n"
                   "MBG_X\tchr19_0049000000_0049999999\trs492602\tFUT2\n"
                   "MBG_Y\tchr2_0136000000_0136999999\trs1\t\n")
    out = Path(tmp_path) / "results"
    assert main(["run", "--genes-file", str(eff), "--out", str(out)]) == 0
    rows = list(csv.DictReader((out / "chains.tsv").open(), delimiter="\t"))
    assert len(rows) == 1                          # NA/empty gene skipped
    assert rows[0]["trait_id"] == "MBG_X"
    assert rows[0]["locus_id"] == "chr19_0049000000_0049999999"
    assert rows[0]["lead_snp"] == "rs492602"
