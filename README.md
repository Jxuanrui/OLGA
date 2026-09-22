# OLGA

Oligogenic Locus-resolved Genetic Attribution: attributing microbiome GWAS
signals to effector genes and host cell types.

Microbial GWAS are oligogenic — in MiBioGen, 17 of 131 taxa carry an
independent genome-wide significant locus and 16 of those carry exactly one. Methods
that score the top-1000 genes dilute a 1-3 gene signal among noise. OLGA
works at the locus instead: define the independent loci, anchor the
effector gene by eQTL colocalisation, attribute it to a cell type, and
arbitrate by independent replication rather than by a colocalisation score.

## Install

```bash
pip install .
olga verify
```

`olga verify` checks the bundled reference packs, lineage rules and four
golden attributions on a fresh install.

## Usage

Attribute a list of genes:

```bash
olga list-reference
olga run --genes FUT2,MCM6,LCT --out results/
```

Attribute every locus in an OLGA effector table (columns `trait_id`,
`locus_id`, `lead_snp`, `effector_gene`); the chain columns carry through
to the output:

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

- Attribution uses expression specificity (tau), a proxy for the causal
  cell type. Cell-type-specific eQTL would settle it; for gut this has only
  recently become available (IBDverse, Alegbe et al., Nature 2026) and is
  not yet wired in.
- Coverage follows the reference packs. LCT is absent from the adult UC
  atlas (adults do not express lactase) and is recovered only through the
  fetal and small-intestine packs.
- A locus can reach PP.H4 above 0.97 and still fail replication. The truth
  of the GWAS signal, not the colocalisation statistic, is the binding
  constraint; run a replication cohort whenever one exists.

## Citation

A methods manuscript is in preparation. Until then, cite the data resources:
MiBioGen (Kurilshikov et al., Nat Genet 2021); Smillie et al., Cell 2019;
Elmentaite et al., Cell 2021; Tabula Sapiens Consortium, Science 2022.

## License

MIT.
