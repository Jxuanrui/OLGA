# olga

**OLGA framework** — Oligogenic Locus-resolved Genetic Attribution: attributing
microbiome GWAS signals to effector genes and host cell types.

**olga software package** (this repository) — the pip-installable implementation
of the framework's cell-type attribution stage.

Microbial GWAS are oligogenic — in MiBioGen, 17 of 131 taxa carry an
independent genome-wide significant locus and 16 of those carry exactly one.
Methods that score the top-1000 genes dilute a 1-3 gene signal among noise.
The OLGA framework works at the locus instead: define the independent loci,
anchor the effector gene by eQTL colocalisation, attribute it to a cell type,
and arbitrate by independent replication rather than by a colocalisation score.

## Scope

This package implements **one stage** of that framework. The stages split as
follows (all are described in the manuscript):

| Framework stage | Where it lives |
|---|---|
| Locus resolution (LD clumping), effector-gene inference | research pipeline (plink2/MAGMA workflows; not in this package) |
| QTL colocalisation (coloc.abf / SuSiE / ColocBoost, SMR) | research pipeline (R/Python workflows over GTEx eQTL Catalogue; not in this package) |
| Independent replication & evidence arbitration | research pipeline (cohort-matching workflows; not in this package) |
| **Cell-type attribution** (tau specificity, multi-atlas consensus, evidence tiers) | **this package** |
| **Reference-pack construction** from any annotated h5ad atlas | **this package** |

Practically: you hand `olga run` a list of effector genes (or an effector
table produced by the upstream framework stages), and it attributes them to
cell types across the bundled gut reference packs. It does not perform GWAS,
clumping, or colocalisation. `locus_id` values in whole-chain mode use the
framework convention (`chr<N>_<10-digit start>_<10-digit end>`, GRCh37 1 Mb
bins) and are carried through unchanged.

## Install

```bash
pip install .
olga verify
```

`olga verify` checks the bundled reference packs, lineage rules and four
golden attributions on a fresh install. Runtime dependencies are numpy,
matplotlib and anndata (anndata is only needed by `build-reference`).

## Usage

Attribute a list of gene symbols (HGNC symbols; a gene absent from every
reference pack is reported as `gene_not_in_reference`, not an error):

```bash
olga list-reference
olga run --genes FUT2,MCM6,LCT --out results/
```

Attribute every locus in an OLGA effector table (columns `trait_id`,
`locus_id`, `lead_snp`, `effector_gene`); the chain columns carry through
to the output. A plain two-column-free text file also works: one gene symbol
per line, `#` comments allowed, and rows with an empty or `NA` effector gene
are skipped (the count is printed to stderr):

```bash
olga run --genes-file effector_genes.tsv --out results/
```

Use a subset of reference packs:

```bash
olga attribute --genes FUT2 --reference uc_gut_smillie ts_colon --out results/
```

Outputs, per run:

- `chains.tsv` — gene, trait, locus, lead SNP, consensus cell state,
  lineage, evidence tier, per-pack detail
- `celltypes/<gene>.celltype.{png,pdf}` — tau per reference pack
- `evidence_tiers.png`, `run_manifest.json`

### Reading the output

Each gene is attributed to the cell state with its highest mean expression in
each reference pack; **tau** (Yanai et al. 2005) scores how specific that
profile is: `tau = sum_i(1 - x_i/x_max)/(n-1)`, 0 = uniformly expressed across
all states, 1 = expressed in exactly one state. Consensus across packs is
judged at the lineage level (atlases disagree on epithelial subtypes far more
than on compartments); `n_refs_covering` / `n_agree` count how many packs
cover the gene and agree with the consensus, and `per_pack_detail` shows each
pack's state and tau. A gene covered by no pack reports `gene_not_in_reference`
— absence of evidence, not evidence of absence.

## Reference packs

Four packs ship with the package (derived statistics, about 46 MB total;
raw atlases are not redistributed):

| pack | states | source |
|---|---|---|
| `uc_gut_smillie` | 51 | Smillie et al., Cell 2019 |
| `fetal_gut_developing` | 21 | Elmentaite et al., Cell 2021 |
| `ts_colon` | 28 | Tabula Sapiens, Science 2022 |
| `ts_small_intestine` | 34 | Tabula Sapiens, Science 2022 |

Build a pack from any annotated h5ad (multiple atlases merge on the gene
union with per-atlas state prefixes):

```bash
olga build-reference --atlas epi.h5ad --atlas lp.h5ad --label Epi --label LP \
    --name my_gut --out packs/my_gut   # --cluster-col auto-detected if omitted
olga run --genes FUT2 --refs packs --out results/
```

Point `--refs` at a directory of pack folders, or set the `OLGA_REFS`
environment variable to make it the default reference root. A pack folder
must contain `manifest.json` and `tau_matrix.tsv` (`markers_strict.tsv` and
`expr_by_cluster.tsv` ship with the bundled packs; other directories in the
refs root are ignored).

## Evidence tiers

Consensus is judged at the lineage level, because atlases disagree on
epithelial subtypes far more than on compartments:

| tier | meaning |
|---|---|
| `cell_type_confirmed` | every covering reference agrees on the lineage |
| `cell_type_consensus_partial` | at least two references agree on the lineage, at least one disagrees |
| `cell_type_single` | only one reference covers the gene |
| `cell_type_discordant` | fewer than two references agree |
| `gene_not_in_reference` | gene absent from every reference |

## Limitations

- The package covers the attribution stage only (see Scope). Framework
  stages that need summary statistics, LD references or QTL data run in the
  research pipeline; the manuscript is their reference description.
- Attribution uses expression specificity (tau), a proxy for the causal
  cell type. Cell-type-specific eQTL would settle it; for gut this has only
  recently become available (IBDverse, Alegbe et al., Nature 2026) and is
  not yet wired in.
- Coverage follows the reference packs. LCT is absent from the adult UC
  atlas (adults do not express lactase) and is recovered only through the
  fetal and small-intestine packs.
- The framework's locus and colocalisation stages use GRCh37 summary
  statistics with ancestry-matched LD references (mismatched reference
  ancestry inflates gene-level signals; demonstrated for EUR references on
  East Asian cohorts). This package itself uses no LD or ancestry data.
- A locus can reach PP.H4 above 0.97 and still fail replication. The truth
  of the GWAS signal, not the colocalisation statistic, is the binding
  constraint; run a replication cohort whenever one exists.

## Troubleshooting

- `no reference packs found` — pass `--refs DIR` or set `OLGA_REFS`; a pip
  install includes the four bundled packs, an editable checkout needs the
  packaged data.
- `no cluster column found` (`build-reference`) — pass `--cluster-col` with
  the obs column holding cell states (auto-detection tries `cell_type`,
  `Cluster`, `celltype`, `cluster`, ...).
- Gene always `gene_not_in_reference` — check the symbol is a current HGNC
  symbol matching the atlas's gene annotation (Ensembl var_names are
  resolved via `feature_name`/`gene_symbols` when present).

## Citation

A methods manuscript is in preparation. Until then, cite the data resources:
MiBioGen (Kurilshikov et al., Nat Genet 2021); Smillie et al., Cell 2019;
Elmentaite et al., Cell 2021; Tabula Sapiens Consortium, Science 2022.

## License

MIT.
