# GLOWr 0.2.3

Changes for whole-genome-sequencing scans, after Alberto Brusati's Terra runs of chromosomes
with eleven million variant records. Results are unchanged wherever the previous version
completed.

- **A variant index for region lookup.** `build_variant_index()` reads a GDS file's
  chromosome, position and id vectors once. `extract_variant_set()` and `glow_region_test()`
  take it as `variant_index` and find a region's records by binary search, setting the GDS
  filter by file position, instead of re-reading the chromosome's variant table and matching
  ids for every region. On a chromosome of one million records a 2 kb window cost 0.40 s in
  the lookup before and about 0.02 s after; the saving grows with the chromosome's record
  count. A call without an index builds one for itself and selects exactly as before.
  `count_index_records()` counts the records of an index inside genomic spans.
- **A faster LD step.** `filter_variants_ld()` computes the correlation of its NA-free
  genotype matrix as one BLAS cross-product of the standardized columns instead of
  `cor(use = "pairwise.complete.obs")`, and its greedy prune keeps each variant's count of
  above-threshold partners up to date instead of recomputing the whole threshold mask after
  every removal. The selection rule and the tie-breaking are unchanged, so the kept variants
  are the same. A gene of 2,018 rare SNVs took 169 s in this step before. A matrix with
  missing values keeps the previous path.
- **A linear-dependence step that does not pay for the columns it drops.** After the LD
  prune, `filter_variants_ld()` drops the columns that are exact linear combinations of the
  columns before them. It did this with LINPACK `qr()`, which cycles every dropped column to
  the end of the matrix by shifting the remaining columns, so a gene with far more rare
  variants than samples (14,038 columns, 3,195 kept, at 3,202 samples) spent 367 s there.
  The step is now a blocked in-order Gram-Schmidt with `qr()`'s rule and tolerance (a column
  is dropped when its residual against the kept columns before it is below 1e-07 of its
  norm), in BLAS calls, stopping once the kept count reaches the number of samples carrying a
  variant. The kept columns are the same, checked on every chromosome 22 gene of the 1000
  Genomes cohort; that gene's step takes a few seconds. A matrix with missing values now
  stops with a message, where `qr()` stopped with an opaque one.
- **A fallback where SPA returns no p-value.** `SPAtest::ScoreTest_SPA()` skips a column whose
  `min(sum(g), sum(2 - g))` is below its minimum, a guard written for dosages in [0, 2]. A
  collapsed burden column (the row sum of many ultra-rare variants) with a mean above 1 fails
  it, came back NA, and `glow_test()` stopped on the gene. SPAtest's own computation uses the
  normal approximation for such a column past the guard, so `getZ_marg_score_binary_SPA()` now
  gives it the standard score Z, warns with the count, and returns the count as
  `n_spa_fallback`. This matches the GLOW methodology paper, which applied SPA to the
  non-collapsed variants only.
- **`marginal_scan()` gains `n_cores`** (default 1): the variant chunks are spread over forked
  workers, each with its own GDS handle, and the result table is the same as with one core.
  Not available on Windows.

# GLOWr 0.2.2

- **A quick-start vignette for the FAVOR annotator.** `vignette("favor-annotation-quick-start")`
  walks a new user from when to use `annotate_favor()`, through what to have ready and which
  settings to choose, to what the call returns and how to check it. It points into the full
  vignette for the evidence, the cost table and the matching procedure, which are unchanged.
  Written from user feedback. No code changed.

# GLOWr 0.2.1

- **Missing QC values no longer stop a region.** `SeqArray::seqVCF2GDS()` stores a VCF
  `FILTER` of `"."` (filters not applied) as `NA` in `annotation/filter`. The variant
  filter compared that value with `"PASS"`, obtained `NA`, and the region stopped with
  "missing value where TRUE/FALSE needed". A missing QC value now counts as not
  passing, so such variants are excluded like any other non-passing variant, and
  `?variant_filter` says so. A GDS converted from PLINK carries `PASS` for every
  variant and is unaffected. Reported by Alberto Brusati on an aGDS built from a VCF.
- **PI evaluation.** `evaluate_PI_models()` and `plot_PI_roc()` score tied predictions
  correctly: the AUC is the rank-based statistic with ties at half credit, and the ROC curve
  takes one point per distinct threshold. Before, an ensemble member whose LASSO kept no
  feature (a constant prediction) scored an AUC of 1.0 whenever the cases were listed first,
  which is what the earlier ALS evaluations' "AUC = 1" models were.

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
