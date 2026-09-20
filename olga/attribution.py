"""Multi-atlas cell-type attribution with lineage consensus and evidence tiers.

Different atlases name the same biology differently ("enterocyte" vs
"enterocyte of epithelium" vs "colon epithelial cell"), so consensus is computed
at the level of coarse LINEAGES derived from the state name, not by string
equality. The evidence tier then records how much independent atlas support the
attribution has, which is how the known limitation (expression specificity is a
proxy for the causal cell type) is made explicit rather than hidden.
"""
from __future__ import annotations

from collections import Counter

from .reference import Pack

# keyword -> lineage (ordered: first match wins for the broad classes)
_LINEAGE_RULES = [
    ("enterocyte", "epithelial_absorptive"),
    ("paneth", "epithelial_secretory"),
    ("goblet", "epithelial_secretory"),
    ("tuft", "epithelial_secretory"),
    ("enteroendocrine", "epithelial_secretory"),
    ("m cell", "epithelial_secretory"),
    ("best4", "epithelial_secretory"),
    ("colon epithelial", "epithelial_other"),
    ("epithelial", "epithelial_other"),
    ("stem", "stem_progenitor"),
    ("progenitor", "stem_progenitor"),
    ("cycling", "proliferative"),
    ("ta ", "proliferative"),
    ("t cell", "immune_t"),
    ("treg", "immune_t"),
    ("cd4", "immune_t"),
    ("cd8", "immune_t"),
    ("nk", "immune_nk"),
    ("b cell", "immune_b"),
    ("plasma", "immune_b"),
    ("monocyte", "myeloid"),
    ("macrophage", "myeloid"),
    ("dc", "myeloid"),
    ("mast", "myeloid"),
    ("neutrophil", "myeloid"),
    ("follicular", "immune_other"),
    ("gc", "immune_other"),
    ("fibroblast", "stromal"),
    ("myofibroblast", "stromal"),
    ("wnt", "stromal"),
    ("rspo", "stromal"),
    ("pericyte", "stromal"),
    ("mesothelial", "stromal"),
    ("endothelial", "endothelial"),
    ("microvascular", "endothelial"),
    ("venule", "endothelial"),
    ("glia", "neural"),
    ("neuron", "neural"),
    ("neural", "neural"),
    ("smooth muscle", "muscle"),
    ("myocyte", "muscle"),
    ("mesodermal", "mesoderm"),
    ("erythro", "blood"),
    ("hematopoietic", "immune_other"),
]


def lineage_of(cell_state: str) -> str:
    s = (cell_state or "").lower()
    for key, lin in _LINEAGE_RULES:
        if key in s:
            return lin
    return "other"


# coarse level at which cross-atlas consensus is computed: atlases disagree on
# epithelial subtypes (enterocyte vs tuft vs paneth) far more than on the broad
# compartment, so agreement is judged here and the fine state is reported.
_SUPER = {
    "epithelial_absorptive": "epithelial",
    "epithelial_secretory": "epithelial",
    "epithelial_other": "epithelial",
    "stem_progenitor": "epithelial",
    "proliferative": "proliferative",
    "immune_t": "immune", "immune_b": "immune", "immune_nk": "immune",
    "myeloid": "immune", "immune_other": "immune",
    "stromal": "stromal", "endothelial": "endothelial",
    "neural": "neural", "muscle": "muscle", "mesoderm": "other", "blood": "immune",
}


def super_lineage_of(cell_state: str) -> str:
    return _SUPER.get(lineage_of(cell_state), "other")


def attribute_gene(gene: str, packs: list[Pack]) -> dict:
    per_pack = {}
    for p in packs:
        a = p.attribute(gene)
        if a:
            per_pack[p.name] = a
    if not per_pack:
        return {"gene": gene, "n_refs_covering": 0, "per_pack": {},
                "consensus_lineage": None, "consensus_state": None,
                "n_agree": 0, "tau_mean": None,
                "evidence_tier": "gene_not_in_reference"}
    supers = [super_lineage_of(a["cell_state"]) for a in per_pack.values()]
    counts = Counter(supers)
    cons_super, n_agree = counts.most_common(1)[0]
    same = [a for a in per_pack.values() if super_lineage_of(a["cell_state"]) == cons_super]
    rep = max(same, key=lambda a: (a["tau"] or 0.0))
    taus = [a["tau"] for a in per_pack.values() if a["tau"] is not None]
    if len(per_pack) == 1:
        tier = "cell_type_single"
    elif n_agree == len(per_pack):
        tier = "cell_type_confirmed"
    elif n_agree >= 2:
        tier = "cell_type_consensus_partial"
    else:
        tier = "cell_type_discordant"
    return {"gene": gene, "n_refs_covering": len(per_pack), "per_pack": per_pack,
            "consensus_lineage": lineage_of(rep["cell_state"]),
            "consensus_super_lineage": cons_super,
            "consensus_state": rep["cell_state"],
            "n_agree": n_agree, "tau_mean": round(sum(taus) / len(taus), 4) if taus else None,
            "evidence_tier": tier}
