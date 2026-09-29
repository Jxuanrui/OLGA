"""Cell-type attribution tests on the bundled reference packs."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # checkout fallback

from olga.attribution import attribute_gene, lineage_of, super_lineage_of  # noqa: E402
from olga.reference import Pack, load_packs  # noqa: E402


@pytest.fixture(scope="module")
def packs():
    p = load_packs()
    assert len(p) >= 4, "expected the four bundled reference packs"
    return p


def test_gene_absent_from_every_reference(packs):
    r = attribute_gene("NOTAGENE123", packs)
    assert r["n_refs_covering"] == 0
    assert r["evidence_tier"] == "gene_not_in_reference"
    assert r["consensus_state"] is None


def test_positive_controls_recover(packs):
    r = attribute_gene("FUT2", packs)
    assert r["n_refs_covering"] >= 2
    assert super_lineage_of(r["consensus_state"]) == "epithelial"
    assert r["evidence_tier"] in ("cell_type_confirmed",
                                  "cell_type_consensus_partial")
    m = attribute_gene("MCM6", packs)
    assert super_lineage_of(m["consensus_state"]) == "immune"  # v0.1.3: corrected from proliferative (word-boundary fix)


def test_lineage_rules_ordering_and_unknown():
    assert lineage_of("Epi|Enterocytes") == "epithelial_absorptive"
    assert lineage_of("tuft cell of colon") == "epithelial_secretory"
    assert lineage_of("plasma cell") == "immune_b"
    assert lineage_of("something totally new") == "other"
    # documented quirk: rules are first-match, so a state containing both
    # "enterocyte" and "best4" resolves as absorptive; super-lineage unaffected
    assert lineage_of("BEST4+ enterocytes") == "epithelial_absorptive"
    assert super_lineage_of("BEST4+ enterocytes") == "epithelial"


def test_pack_tau_and_marker_reading(packs):
    uc = next(p for p in packs if p.name == "uc_gut_smillie")
    a = uc.attribute("FUT2")
    assert a is not None and a["cell_state"]
    assert 0.0 <= (a["tau"] or 0.0) <= 1.0
    assert uc.attribute("NOTAGENE123") is None


def test_tau_nan_encoding_is_treated_as_missing(packs):
    """Packs built by the research pipeline write literal 'nan'; the tool
    writes ''. Both must behave as no-tau, never poison tau_mean."""
    r = attribute_gene("SOX1", packs)          # 'nan' tau in uc_gut_smillie
    assert isinstance(r["tau_mean"], float) or r["tau_mean"] is None
    if isinstance(r["tau_mean"], float):
        assert r["tau_mean"] == r["tau_mean"]  # not NaN


def test_atlas_disagreement_downgrades_tier(tmp_path):
    """Two single-state packs disagreeing at super level -> discordant."""
    import json

    def fake_pack(name, state, tau):
        d = tmp_path / name
        d.mkdir()
        (d / "tau_matrix.tsv").write_text(
            "gene\ttau\ttop_cell_state\ttop_expr\tn_states\n"
            f"GENE\t{tau}\t{state}\t1.0\t1\n")
        (d / "markers_strict.tsv").write_text("gene\tcell_state\tscore\tmargin\n")
        (d / "manifest.json").write_text(json.dumps(
            {"name": name, "n_cell_states": 1, "n_genes": 1}))
        return Pack(d)

    p1 = fake_pack("p_epi", "Epi|Enterocytes", 0.9)
    p2 = fake_pack("p_imm", "LP|T cell", 0.8)
    assert attribute_gene("GENE", [p1, p2])["evidence_tier"] == "cell_type_discordant"
    assert attribute_gene("GENE", [p1])["evidence_tier"] == "cell_type_single"


def test_word_boundary_lineage_regressions():
    """v0.1.3 regressions: state names containing 'ta ', 't cell' or 'm cell'
    as substrings must not be misclassified (rule set frozen after the
    134-state audit; see repository history for the audit table)."""
    cases = {
        "CD4-positive, alpha-beta T cell": "immune_t",      # was proliferative
        "gamma-delta T cell": "immune_t",                   # was proliferative
        "naive thymus-derived CD4-positive, alpha-beta T cell": "immune_t",
        "myofibroblast cell": "stromal",                    # was immune_t
        "neural crest cell": "neural",                      # was immune_t
        "mast cell": "myeloid",                             # was immune_t
        "intestinal crypt stem cell of colon": "stem_progenitor",  # was secretory
        "germ cell": "other",                               # must NOT be secretory
        "TACIT cell": "other",                              # must NOT be immune_t
        "Epi|TA 1": "proliferative",                        # anchored 'ta' still hits
        "Epi|Secretory TA": "proliferative",                # adjudicated
        "Epi|M cells": "epithelial_secretory",              # anchored 'm cell' still hits
        "cDC2": "myeloid",                                  # substring 'dc' retained
    }
    for state, want in cases.items():
        assert lineage_of(state) == want, state


def test_silent_gene_guard(packs):
    """v0.1.3 regression: a gene is dropped only when specificity is missing
    AND expression is zero (AND, not OR)."""
    uc = next(p for p in packs if p.name == "uc_gut_smillie")
    tc = next(p for p in packs if p.name == "ts_colon")
    assert uc.attribute("AANAT") is None            # silent: tau nan + expr 0
    apo = tc.attribute("APOA4")                     # finite tau, expr prints 0.0000
    assert apo is not None and apo["tau"] == 1.0    # low-expression specific: KEPT
