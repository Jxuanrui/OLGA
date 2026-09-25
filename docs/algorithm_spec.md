# The OLGA framework — method specification (public summary)

This document specifies the **OLGA framework** (Oligogenic Locus-resolved
Genetic Attribution) at the level of detail needed to understand what the
`olga` package implements and what the remaining stages do. The manuscript is
the authoritative description of the framework; internal analysis scripts are
described at the stage level only.

**Terminology.** The *OLGA framework* is the full research method. The *olga
software package* (this repository) implements its cell-type attribution stage
and the reference-pack tooling. Nothing in this package performs GWAS,
clumping, colocalisation, or replication.

## Stage overview

```
GWAS summary statistics (per taxon)
  -> Stage 1  QC + gene-level statistics (sentinel filter, region blacklist)
  -> Stage 2  locus resolution (LD clumping, preregistered parameters)
  -> Stage 3  effector-gene inference (highest gene-level Z in the locus window)
  -> Stage 4  allele harmonisation across GWAS / QTL / LD (three-way, by name)
  -> Stage 5  QTL colocalisation evidence (coloc.abf primary; SuSiE /
              ColocBoost / SMR / pQTL secondary & validation lines)
  -> Stage 6  cell-type attribution (tau specificity + multi-atlas consensus)  [THIS PACKAGE]
  -> Stage 7  independent replication (label-matched paired tests, direction +
              significance + heterogeneity) — the sole arbiter of locus truth
  -> Stage 8  evidence integration (reporting; colocalisation scores never
              arbitrate locus truth)
```

## Stage parameters (frozen)

- **Locus resolution**: plink2 clumping, index P<5e-8, secondary P<1e-5,
  r²<0.1, 500 kb; loci are 1 Mb GRCh37 bins anchored at the lead SNP.
- **Sentinel filter**: gene-level rows with the un-computable sentinel
  (|Z|≈8.1259 with P≥0.99) and a telomere/centromere/HLA region blacklist are
  excluded everywhere.
- **Effector gene**: the highest-Z valid gene overlapping the locus window;
  zero-candidate loci are reported as such, never imputed.
- **Colocalisation**: Wakefield ABF (coloc defaults: p1=p2=1e-4, p12=1e-5),
  ≥30 shared allele-matched variants, per-gene tests (never aggregated);
  SuSiE/ColocBoost runs align LD and summary statistics by variant NAME after
  three-way allele harmonisation.
- **Cell attribution** (this package): tau specificity on mean de-logged
  normalised expression; attribution state = argmax mean expression per pack;
  consensus at the lineage level across packs; five evidence tiers
  (`confirmed / consensus_partial / single / discordant / gene_not_in_reference`).
- **Replication**: a locus replicates when the lead SNP shows the same effect
  direction with P<0.05 in a label-matched independent cohort trait; cross-trait
  best hits are exploratory only.

## Scientific boundaries (enforced in all OLGA writing)

1. "Oligogenic" describes the **observed detectable architecture** at current
   cohort sizes, not a claim about the unobservable true architecture.
2. A high colocalisation posterior is consistency with a shared-variant
   hypothesis, **not** proof of a causal gene; the framework's negative
   controls (loci with H4>0.97 that failed replication) exist to keep this
   honest.
3. LCT/MCM6 and FUT2 are **positive controls** (method validation), not de
   novo discoveries.
4. LD references must match cohort ancestry (mismatched references inflate
   gene-level signals; demonstrated empirically).
5. `gene_not_in_reference` is absence of evidence, not evidence of absence.

## Provenance expectations

Any result reported by the framework carries: input checksums, tool versions,
frozen parameters, and random seeds where sampling occurs. The package itself
is deterministic (byte-identical outputs on repeat runs).
