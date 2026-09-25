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


def verify() -> int:
    checks: list[tuple[bool, str]] = []

    import matplotlib  # noqa: F401  (outputs hard-fail without it)
    import numpy
    checks.append((True, f"core imports ok (numpy {numpy.__version__}, "
                         f"matplotlib {matplotlib.__version__})"))

    packs_meta = list_packs()
    checks.append((len(packs_meta) >= 4,
                   f"reference packs found: {len(packs_meta)} (expect >= 4 bundled)"))
    packs = load_packs()
    for m in packs_meta:
        man_genes = m.get("n_genes")
        tau_rows = sum(1 for _ in open(Path(m["_path"]) / "tau_matrix.tsv")) - 1
        checks.append((man_genes == tau_rows,
                       f"{m['name']}: manifest n_genes={man_genes} == tau rows={tau_rows}"))
        # markers_strict is optional on disk but load-bearing for the
        # marker_state attribution field — a pack shipped without it must not
        # pass silently
        mk = Path(m["_path"]) / "markers_strict.tsv"
        mk_rows = (sum(1 for _ in open(mk)) - 1) if mk.exists() else -1
        checks.append((mk_rows > 0,
                       f"{m['name']}: markers_strict rows={mk_rows} (>0 required)"))
        expr = Path(m["_path"]) / "expr_by_cluster.tsv"
        if expr.exists():
            with expr.open() as fh:
                n_cols = len(fh.readline().rstrip("\n").split("\t")) - 1
            checks.append((n_cols == m.get("n_cell_states"),
                           f"{m['name']}: expr_by_cluster states={n_cols} "
                           f"== manifest {m.get('n_cell_states')}"))
        else:
            checks.append((False, f"{m['name']}: expr_by_cluster.tsv missing"))

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
