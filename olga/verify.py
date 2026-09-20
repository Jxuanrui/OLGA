"""Self-check of the installed tool: reference-pack integrity, lineage rules and
bundled golden attributions. ``olga verify`` is a research-tool function, not a
test suite: it ships with the package so any user (or CI) can confirm their
install produces the documented behaviour on the bundled reference packs.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

from .attribution import attribute_gene, lineage_of
from .reference import load_packs, list_packs

GOLDEN = Path(__file__).resolve().parent / "data" / "verify_golden.json"

LINEAGE_CHECKS = {
    "Epi|Enterocytes": "epithelial_absorptive",
    "LP|Cycling T": "proliferative",
    "enterocyte of epithelium proper of jejunum": "epithelial_absorptive",
    "tuft cell of colon": "epithelial_secretory",
    "plasma cell": "immune_b",
}


def verify() -> tuple[int, list[str]]:
    checks: list[tuple[bool, str]] = []

    packs_meta = list_packs()
    checks.append((len(packs_meta) >= 4,
                   f"reference packs found: {len(packs_meta)} (expect >= 4 bundled)"))
    packs = load_packs()
    for m in packs_meta:
        man_genes = m.get("n_genes")
        tau_rows = sum(1 for _ in open(Path(m["_path"]) / "tau_matrix.tsv")) - 1
        checks.append((man_genes == tau_rows,
                       f"{m['name']}: manifest n_genes={man_genes} == tau rows={tau_rows}"))

    for state, want in LINEAGE_CHECKS.items():
        got = lineage_of(state)
        checks.append((got == want, f"lineage_of({state!r}) == {want!r} (got {got!r})"))

    golden = json.loads(GOLDEN.read_text())
    for g in golden["genes"]:
        res = attribute_gene(g["gene"], packs)
        ok = (res["evidence_tier"] == g["tier"]
              and res.get("consensus_super_lineage") == g["super_lineage"]
              and res["consensus_state"] == g["consensus_state"])
        checks.append((ok, f"golden {g['gene']}: {g['consensus_state']}"
                           f"/{g['super_lineage']}/{g['tier']}"))

    n_fail = sum(0 if ok else 1 for ok, _ in checks)
    for ok, msg in checks:
        print(("PASS " if ok else "FAIL ") + msg)
    print(f"olga verify: {len(checks) - n_fail}/{len(checks)} checks passed")
    return (0 if n_fail == 0 else 1)
