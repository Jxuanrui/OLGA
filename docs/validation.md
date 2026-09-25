# Validation coverage

## What CI covers (this repository, every push/PR)

- `tests/test_attribution.py` — attribution semantics on the bundled packs:
  absent genes, positive-control recovery (FUT2/MCM6), lineage-rule table,
  tau NaN encoding, cross-atlas disagreement downgrading.
- `tests/test_buildref.py` — `build-reference` edge cases, including the two
  crash regressions found in the pre-manuscript audit (single-state and
  zero-state packs).
- `tests/test_cli.py` — effector-table parsing (empty/NA skipping, plain gene
  lists, trait-table misuse), chains.tsv schema, clear error on empty input.
- `tests/test_ld_alignment.py` — the LD-ordering protocol lock: name-based
  realignment recovers the true matrix from shuffled files, and set mismatch
  is caught by guards instead of silently misaligning.
- `tests/test_smoke_run.py` — README example end-to-end on bundled packs.
- `olga verify` — 23 checks: pack integrity (tau/markers/expression vs
  manifest), core imports, lineage rules, four golden attributions.

## Extended validation (research environment only)

Some regression tests intentionally live with the research pipeline because
they import private-pipeline modules or need cohort data / R / plink2 / the
1000G reference. They are documented here so their existence is public:

- per-variant discovery-statistic direction (same lead SNP in many trait
  GWAS with opposite effects must use the locus's own trait);
- end-to-end shuffle-invariance of the actual coloc.susie R runner
  (byte-identical output when only LD-file SNP order changes);
- bidirectional normalised-label matching for cohort replication pairing
  (both containment directions; empty labels never match);
- coloc.abf math (shared/distinct/no-signal posteriors; minimum-shared guard);
- MAGMA sentinel-row filter edges;
- matched-locus sampler strata/calipers;
- full four-cohort replication rerun + byte-comparison against frozen
  artifacts, and determinism (repeat runs byte-identical, including figures).

These run in the research repository's environment; their absence from CI is
a data-access constraint, not a coverage decision.
