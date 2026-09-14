# GLOWr 0.2.0

`annotate_favor()` is rewritten around a corrected matching rule and gains a second database
backend. The returned table, the aGDS nodes and the provenance change shape, hence the minor
version. The arguments of 0.1.x keep their names and positions.

- **Corrected variant matching.** Under `match_method = "flexible"` the lookup key of an SNV is
  normalized against FAVOR's reference base, by swapping REF and ALT or complementing both
  alleles. Every input row records its outcome in `match_tier` (`exact`, `swapped`,
  `flipped`, `flipped_swapped`, `position_only`, or a reason: `rsid_conflict`, `unmatched_alt`,
  `unmatched_ref`, `uncovered`, `no_database`, `unsupported_allele`) and the database row it
  used in `favor_key`. The 0.1.x flexible matching had no strand tier: on genotyping-chip data
  half the SNVs fell through to the mean over the three alternative alleles at their position
  and stored a value that belongs to no variant. A matched variant now always carries its own
  row's values. Indels match by the exact key only. The cohort's alleles, genotypes and keys are
  never modified.
- **The rsID as evidence.** For every matched row the input's rsID (the GDS `annotation/id`,
  or an `rsID` column of a data.frame input) is compared with the FAVOR row's rsID and recorded
  in `rsid_check` (`same`, `differs`, `chip_none`, `favor_none`). The new argument
  `rsid_policy` is `"require"` by default: a swapped or flipped match that the rsID does not
  confirm is withheld as `rsid_conflict`, so a rewritten key is accepted only when dbSNP
  confirms it or the input carries no rsID to check. Exact matches are never withheld.
  `"record"` keeps every match and only records the check. The rsID is never used to find a
  row.
- **The annotation catalog.** Three entries of the internal name catalog that maps STAAR-style
  names to aGDS nodes pointed at nodes the annotator does not write. `aPC.Protein` now names
  `apc_protein_function_v3`, `aPC.Conservation` names `apc_conservation_v2` and
  `aPC.LocalDiversity` names `apc_local_nucleotide_diversity_v3`, the versioned nodes the FAVOR
  databases carry.
- **FAVOR 2.0 Parquet backend.** `favor_db_format = "parquet"` (detected under the default
  `"auto"`) reads the FAVOR team's per-chromosome Parquet files directly with the `arrow`
  package (now in Suggests), one row group at a time and only the leaf columns requested. It
  serializes GeneHancer and the GENCODE information fields into the v1 string layout. Both
  backends return identical outcomes for identical content. A field the database lacks is
  skipped with a message. FAVOR 2.0 has no `apc_conservation`,
  `apc_local_nucleotide_diversity`, `apc_local_nucleotide_diversity_v2` or
  `apc_proximity_to_coding`.
- **Position-only keys** (`CHR-POS-NA-NA`) are averaged over the SNV rows at the position only;
  a position with no SNV row is `uncovered`.
- **aGDS nodes and provenance.** The aGDS carries `annotation/info/favor_match_tier` and
  `annotation/info/favor_rsid_check` beside the `FunctionalAnnotation` folder, whose attributes
  record the database format and files (with their SHA-256 when a `SHA256SUMS` file sits
  beside them), the `favor_release` label, the GLOWr version, the matching settings, the rsID
  policy and source, and the counts per outcome and per check. The writers align rows to the
  GDS through `variant.id`.
- **Interface.** New arguments `favor_db_format`, `favor_release` and `rsid_policy`, placed after
  `verbose` so that positional calls written for 0.1.x keep their meaning. GDS input returns
  `variant_id` and `rsID` beside `VarInfo`. A VarInfo position must be a positive integer in
  decimal digits. A fractional or zero position used to be truncated onto another coordinate.
  The CSV backend fetches the rows at the input positions by an xsv join on `position`, or by
  `fread()` with a filter. xsv file paths are shell-quoted.
- **Documentation.** The new vignette `favor-annotation` explains the matching rule, what a tier
  claims, the rsID evidence, position-only keys, the optional PLINK normalization of chip data
  and the measured cost. `inst/scripts/annotate_favor_batch.R` exposes `--favor-db-format`,
  `--favor-release` and `--rsid-policy`.

# GLOWr 0.1.1

- Fixed a multi-chunk bug in `annotate_favor()`'s xsv fast path (exact
  matching): variants falling in the second and later FAVOR chunk files of a
  chromosome silently lost their annotations (the per-chunk left-joins made the
  downstream deduplication keep an empty row instead of the real match). The
  per-chunk joins are now inner joins, and every `xsv` invocation's exit status
  is checked, so a failed join aborts instead of writing a silently
  under-annotated aGDS.
- `annotate_favor()`'s default `features` now covers the full FAVOR Essential DB
  annotation content (30 columns, previously the 20 purely numeric ones), so an
  aGDS built with the defaults carries the categorical features (e.g. GENCODE
  categories) that the built-in coding masks filter on. String features are
  first-class on the R join path: flexible-match aggregation is type-aware
  (numeric features average across matches; character features take the first
  non-empty value), and `na_handling = "zero"` no longer writes `0` into
  character columns.
- New `annotation_predicate()` filter grammar for `variant_filter()`: a clause
  condition may now test `nonempty`/`empty`, `in`/`not_in`, or numeric
  `gt`/`ge`/`lt`/`le` alongside the existing set-membership conditions. Missing
  annotation values never satisfy a positive test.
- `calibrate_pvalues()` gains `inflation_only`: when `TRUE`, a calibration
  factor below 1 is clamped to 1, reproducing the field-standard one-sided
  genomic-control convention (correct inflation, leave deflation alone).
  Default `FALSE` keeps the previous two-sided behavior.
- The compiled code now links BLAS/LAPACK explicitly (`src/Makevars` and
  `src/Makevars.win`); the Rtools toolchain does not auto-link them, so this
  fixes the Windows `R CMD check`.
- Documentation: `PI` is described throughout as the variant-importance score,
  and the methodology reference now cites the GLOW methods paper (Zhang, Liu,
  Landers, and Wu, Annals of Applied Statistics, in revision).

# GLOWr 0.1.0

Initial public release.

- GLOW (inteGrative anaLysis using Optimized Weights) variant-set association tests:
  Burden, SKAT, Fisher-combination, and the Omnibus test that combines them,
  for both continuous and binary phenotypes with covariate adjustment.
- Optimal data-adaptive weighting: estimation of the effect-size model `B` and
  the variant-importance score `PI`, and their combination into per-variant
  weights used by the set-level tests.
- Single-variant (per-variant) score statistics with saddlepoint (SPA)
  calibration for binary traits.
- p-value calibration and genomic-inflation diagnostics, including LD-score
  regression utilities and Cauchy-combination (CCT) aggregation.
- Region definition and variant-set extraction from GDS/aGDS inputs, with FAVOR
  functional-annotation support.
- Built on the GFisher backend for the generalized Fisher-type combination
  p-values.
