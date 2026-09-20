# OLGA

**O**ligogenic **L**ocus-resolved **G**enetic **A**ttribution — attributing
microbiome GWAS signals to effector genes and host cell types.

## Why this exists

Microbial GWAS are **oligogenic**: 17 of 131 MiBioGen taxa carry a genome-wide
significant locus and 16 of those carry exactly one independent locus. Methods
that aggregate the top-1000 genes (scBPS / scDRS / scPagwas style) therefore
dilute a 1–3 gene signal among ~1000 noise genes, and their enrichment tests are
null on real microbiome data. OLGA is locus-resolved instead: it defines the
1–3 independent loci, anchors the effector gene by eQTL colocalisation, and
attributes it to a cell type — then decides with calibrated statistics and
independent replication rather than a colocalisation score.

## Install

```bash
pip install .
# raw single-cell atlases are NOT required; derived reference packs ship with the tool
olga list-reference
```

Four gut reference packs ship as derived statistics (~46 MB total):
`uc_gut_smillie` (51 states), `fetal_gut_developing` (21), `ts_colon` (28),
`ts_small_intestine` (34). Build your own pack from any annotated h5ad:

```bash
olga build-reference --atlas atlas.h5ad --name my_gut --out packs/my_gut
olga attribute --genes FUT2 --refs packs --out out/
```

After installing, run `olga verify` to self-check the reference packs, lineage
rules and golden attributions.

## Use

```bash
# attribute effector genes (usually from the colocalisation step) to cell types
olga run --genes LCT,MCM6,FUT2 --out results/

# inspect a subset of reference packs
olga attribute --genes CHST1 --reference uc_gut_smillie fetal_gut_developing --out out/
```

Outputs (tables + publication figures, **not** a narrative report):

```
results/
├── chains.tsv                  gene, cell_state, lineage, evidence_tier,
│                               n_refs_covering, n_agree, tau_mean, per-pack detail
├── celltypes/<gene>.celltype.{png,pdf}   tau per reference pack, consensus highlighted
├── evidence_tiers.png
└── run_manifest.json
```

## Evidence tiers

Cross-atlas consensus is computed at the coarse lineage level (atlases disagree
on epithelial subtypes far more than on compartments), and reported honestly:

| tier | meaning |
|---|---|
| `cell_type_confirmed` | all covering references agree on the lineage |
| `cell_type_consensus_partial` | majority agree, at least one disagrees |
| `cell_type_single` | only one reference covers the gene |
| `cell_type_discordant` | no majority |
| `gene_not_in_reference` | gene absent from every reference |

## Known limits (stated, not hidden)

* Cell-type attribution uses **expression specificity (tau)**, which is a proxy
  for — not proof of — the causal cell type. Fixing this properly needs
  cell-type-specific eQTL, which does not yet exist at scale for gut.
* Coverage is bounded by the reference. LCT, for example, is absent from the
  adult UC atlas (adults do not express lactase) and is only recovered by adding
  the fetal-gut and small-intestine references.
* **Colocalisation is not truth.** A locus can reach PP.H4 > 0.97 and still fail
  independent replication; the binding constraint is the truth of the GWAS signal
  (winner's curse), not the colocalisation statistic. Always run a replication
  cohort when one exists.

## Citation

A methods manuscript is in preparation; this repository is the reference
implementation. Until then, please cite the underlying data resources:
MiBioGen (Kurilshikov et al., *Nat Genet* 2021), Smillie et al. (*Cell* 2019),
Tabula Sapiens (Tabula Sapiens Consortium, *Science* 2022), SCENIC+ (Bravo
González-Blas et al., *Nat Methods* 2023).

## License

MIT — see [LICENSE](LICENSE).
