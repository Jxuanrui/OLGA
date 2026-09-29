# Changelog

## v0.1.3 (2026-09-29)

Fixes (found by an independent pre-manuscript review; 134-state audit of all
bundled packs, zero unexplained rows):

- **Lineage word-boundary fix (`attribution.py`)**: three short keywords
  (`ta`, `m cell`, `t cell`) were substring-matched and misclassified real
  atlas state names. ERRATUM for v0.1.2 and earlier: state names containing
  the substrings "ta " (e.g. "CD4-positive, alpha-**beta T cell**"),
  "t cell" (e.g. "myofibroblas**t cell**", "neural cres**t cell**",
  "mas**t cell**") or "m cell" (e.g. "ste**m cell**") received the wrong
  lineage, which could change multi-atlas consensus and evidence tiers.
  Most visibly, MCM6's bundled-pack attribution changes from
  proliferative / `LP|Cycling T` to immune / `regulatory T cell`
  (tier unchanged: consensus_partial). "Epi|Secretory TA" is now classified
  proliferative (adjudicated with the audit). Statistical layers (MAGMA,
  colocalisation, replication) are entirely unaffected.
- **Silent-gene guard (`reference.py`)**: genes with no measurable
  specificity AND zero expression are no longer returned as "covered"
  attributions (previously they inherited an arbitrary argmax state and
  inflated `n_refs_covering`). Genes with finite tau are always kept, even
  when `top_expr` prints as 0.0000.
- `olga verify`: 23 -> 28 checks (five word-boundary lineage regressions).
- Tests: 17 -> 19 public cases (word-boundage regressions, silent-gene
  guard); CI now runs an end-to-end `olga run` smoke step asserting the
  FUT2 golden tier.

Known limitations deferred to v0.1.4 (documented, not fixed): 13 atlas
states have no lineage rule and map to `other` (natural killer cell,
myeloid dendritic cell, innate lymphoid cell, transit amplifying cell,
interstitial cell of Cajal, enterochromaffin-like cell, ILCs, MT-hi, myeloid
leukocyte); `build-reference` expects a log1p-normalised (sparse) matrix.

## v0.1.2 (2026-09-25)

Audit-hardened release: build-reference single/zero-state crash fixes,
NaN-tau handling, robust pack discovery, verify 23 checks, dependency trim
(pysam removed), README scope split (framework vs package).

## v0.1.1 (2026-09-24)

Ship `verify_golden.json` in package data; sync `__version__`.

## v0.1.0 (2026-09-24)

Initial release: locus-resolved attribution tool, four gut reference packs,
`olga run/verify/build-reference`.
