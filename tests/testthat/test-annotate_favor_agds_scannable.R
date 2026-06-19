# ==============================================================================
# Regression tests: annotate_favor()'s aGDS output is directly scannable
# ==============================================================================
#
# These tests are the guard for the bug fixed on 2026-06-11: annotate_favor()'s
# three aGDS writers used to store /annotation/info/FunctionalAnnotation as a
# single numeric MATRIX node, which (a) exposes no per-feature sub-nodes -- so the
# variant-set scan's seqGetData(".../FunctionalAnnotation/<feature>") read failed
# -- and (b) coerced mixed numeric/string features to one type, dropping string
# coding nodes (e.g. genecode_comprehensive_category). The writers now produce the
# STAARpipeline sub-node FOLDER format (one native-typed sub-node per feature).
#
# Unlike the other favor tests (which skip without the real-data FAVOR fixtures),
# these build a tiny synthetic SeqArray GDS + a tiny FAVOR DB INLINE, so they run
# anywhere and would have caught the bug.

# ---- Skip when the GDS toolchain is unavailable ----
skip_if_no_gds_tools <- function() {
  for (pkg in c("SNPRelate", "SeqArray", "gdsfmt")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      testthat::skip(paste0("Package '", pkg, "' not available"))
    }
  }
}

# ------------------------------------------------------------------------------
# Inline fixture builders
# ------------------------------------------------------------------------------

# Build a tiny synthetic SeqArray GDS. Returns the path plus the GDS-EXPOSED
# variant fields (chr/pos/ref/alt), because seqSNP2GDS may flip REF/ALT relative
# to the snpgdsCreateGeno input; the FAVOR DB is keyed to the exposed values so
# exact matching is deterministic.
#
# maf controls genotype frequency: a small maf keeps variants "rare" so rare-MAF
# coding filters retain them.
.build_tiny_gds <- function(path, n_variants = 10L, n_samples = 12L,
                            chr = 22L, start_pos = 1000L, step = 100L,
                            maf = 0.05, seed = 1L) {
  set.seed(seed)
  positions <- as.integer(seq(start_pos, by = step, length.out = n_variants))
  geno <- matrix(stats::rbinom(n_variants * n_samples, 2, maf),
                 nrow = n_variants, ncol = n_samples)
  # Guarantee every variant has MAC >= 1 (one heterozygote) so MAC-based filters
  # do not stochastically drop variants -- keeps the per-variant retention
  # deterministic regardless of the random draw, so a coding-mask count is exact.
  geno[cbind(seq_len(n_variants), rep_len(seq_len(n_samples), n_variants))] <-
    pmax(geno[cbind(seq_len(n_variants), rep_len(seq_len(n_samples), n_variants))], 1L)

  tmp_snp <- tempfile(fileext = ".gds")
  on.exit(unlink(tmp_snp), add = TRUE)
  SNPRelate::snpgdsCreateGeno(
    gds.fn = tmp_snp, genmat = geno,
    sample.id = sprintf("s%03d", seq_len(n_samples)),
    snp.id = seq_len(n_variants),
    snp.chromosome = rep(as.integer(chr), n_variants),
    snp.position = positions,
    snp.allele = rep("A/G", n_variants),
    snpfirstdim = TRUE
  )
  suppressMessages(SeqArray::seqSNP2GDS(tmp_snp, path, verbose = FALSE))

  g <- SeqArray::seqOpen(path)
  on.exit(SeqArray::seqClose(g), add = TRUE)
  list(
    path = path,
    chr  = SeqArray::seqGetData(g, "chromosome"),
    pos  = SeqArray::seqGetData(g, "position"),
    ref  = SeqArray::seqGetData(g, "$ref"),
    alt  = SeqArray::seqGetData(g, "$alt")
  )
}

# Build a tiny one-chunk FAVOR DB keyed to a GDS's exposed variants. Half the
# variants are exonic/nonsynonymous (so the missense coding mask selects them),
# the rest intronic. Returns the FAVOR directory path and the feature data.frame
# in GDS variant order (so tests can assert exact values).
.build_tiny_favor_db <- function(dir, gds_fields, seed = 2L) {
  set.seed(seed)
  n <- length(gds_fields$pos)
  varinfo <- paste(gds_fields$chr, gds_fields$pos,
                   gds_fields$ref, gds_fields$alt, sep = "-")

  exonic_cat <- rep(c("nonsynonymous SNV", ""), length.out = n)
  gene_cat   <- ifelse(exonic_cat == "nonsynonymous SNV", "exonic", "intronic")

  favor <- data.frame(
    variant_vcf = varinfo,
    chromosome  = gds_fields$chr,
    position    = gds_fields$pos,
    ref_vcf     = gds_fields$ref,
    alt_vcf     = gds_fields$alt,
    cadd_phred  = round(stats::runif(n, 0, 40), 3),
    linsight    = round(stats::runif(n, 0, 1), 3),
    genecode_comprehensive_category        = gene_cat,
    genecode_comprehensive_exonic_category = exonic_cat,
    stringsAsFactors = FALSE
  )

  data.table::fwrite(favor, file.path(dir, "chr22_1.csv"))
  data.table::fwrite(
    data.frame(Chr = 22L, File_No = 1L, Start_Pos = 1L, End_Pos = 1e7),
    file.path(dir, "FAVORdatabase_chrsplit.csv")
  )
  list(dir = dir, favor = favor)
}

# Feature set spanning both numeric and string (coding) annotations.
.tiny_feats <- function() {
  c("cadd_phred", "linsight",
    "genecode_comprehensive_category",
    "genecode_comprehensive_exonic_category")
}

# ------------------------------------------------------------------------------
# TEST: output_agds (.create_agds_from_gds) -- the 00 -> 03 chain path
# ------------------------------------------------------------------------------

test_that("annotate_favor(output_agds=) writes a scannable sub-node folder aGDS", {
  skip_if_no_gds_tools()

  gds_path <- tempfile(fileext = ".gds")
  on.exit(unlink(gds_path), add = TRUE)
  gf <- .build_tiny_gds(gds_path)

  favor_dir <- tempfile("favordb_"); dir.create(favor_dir)
  on.exit(unlink(favor_dir, recursive = TRUE), add = TRUE)
  fb <- .build_tiny_favor_db(favor_dir, gf)
  feats <- .tiny_feats()

  out_agds <- tempfile(fileext = ".gds")
  on.exit(unlink(out_agds), add = TRUE)
  annotate_favor(variants = gds_path, favor_db_path = favor_dir,
                 features = feats, output_agds = out_agds,
                 match_method = "exact", verbose = 0)

  expect_true(file.exists(out_agds))

  g <- SeqArray::seqOpen(out_agds)
  on.exit(SeqArray::seqClose(g), add = TRUE)

  # (1) FunctionalAnnotation is a FOLDER exposing one sub-node per feature.
  fa_node <- gdsfmt::index.gdsn(g, "annotation/info/FunctionalAnnotation")
  sub_nodes <- gdsfmt::ls.gdsn(fa_node)
  expect_true(length(sub_nodes) > 0)
  expect_setequal(sub_nodes, feats)

  # (2) The scan path reads a NUMERIC feature: right length, right values,
  #     in GDS variant order.
  cadd <- SeqArray::seqGetData(g, "annotation/info/FunctionalAnnotation/cadd_phred")
  expect_length(cadd, length(gf$pos))
  expect_type(cadd, "double")
  expect_equal(as.numeric(cadd), fb$favor$cadd_phred)

  # (3) The scan path reads a STRING coding feature WITHOUT coercion -- the case
  #     a matrix node could not carry.
  gene_cat <- SeqArray::seqGetData(
    g, "annotation/info/FunctionalAnnotation/genecode_comprehensive_category")
  expect_type(gene_cat, "character")
  expect_identical(as.character(gene_cat),
                   fb$favor$genecode_comprehensive_category)
})

test_that("extract_variant_set resolves a coding mask on the annotate_favor aGDS", {
  skip_if_no_gds_tools()

  gds_path <- tempfile(fileext = ".gds")
  on.exit(unlink(gds_path), add = TRUE)
  gf <- .build_tiny_gds(gds_path)

  favor_dir <- tempfile("favordb_"); dir.create(favor_dir)
  on.exit(unlink(favor_dir, recursive = TRUE), add = TRUE)
  fb <- .build_tiny_favor_db(favor_dir, gf)

  out_agds <- tempfile(fileext = ".gds")
  on.exit(unlink(out_agds), add = TRUE)
  annotate_favor(variants = gds_path, favor_db_path = favor_dir,
                 features = .tiny_feats(), output_agds = out_agds,
                 match_method = "exact", verbose = 0)

  # The real scan entry: a missense coding mask reads the string exonic-category
  # sub-node. Thresholds are made non-masking (rare_maf_cutoff = 0.5,
  # min_variants = 1) so the ANNOTATION CLAUSE is what is exercised.
  region <- list(chr = "22", start = min(gf$pos), end = max(gf$pos),
                 label = "TINY")
  spec <- coding_filter("missense", rare_maf_cutoff = 0.5, min_variants = 1L)

  vset <- extract_variant_set(
    out_agds, region, spec,
    annotation_names = c("cadd_phred", "linsight"), verbose = 0
  )

  expect_false(is.null(vset))
  expect_s3_class(vset, "glow_variant_set")

  # The missense set must be exactly the nonsynonymous variants.
  nonsyn_pos <- gf$pos[
    fb$favor$genecode_comprehensive_exonic_category == "nonsynonymous SNV"]
  expect_equal(vset$n_variants, length(nonsyn_pos))
  expect_setequal(vset$variant_info$pos, nonsyn_pos)

  # Annotation scores for PI are carried through.
  expect_equal(colnames(vset$annotations), c("cadd_phred", "linsight"))
  expect_equal(nrow(vset$annotations), vset$n_variants)
})

# ------------------------------------------------------------------------------
# TEST: update_gds = TRUE (.update_gds_with_annotations)
# ------------------------------------------------------------------------------

test_that("annotate_favor(update_gds=TRUE) writes a scannable sub-node folder", {
  skip_if_no_gds_tools()

  gds_path <- tempfile(fileext = ".gds")
  on.exit(unlink(gds_path), add = TRUE)
  gf <- .build_tiny_gds(gds_path)

  favor_dir <- tempfile("favordb_"); dir.create(favor_dir)
  on.exit(unlink(favor_dir, recursive = TRUE), add = TRUE)
  fb <- .build_tiny_favor_db(favor_dir, gf)
  feats <- .tiny_feats()

  annotate_favor(variants = gds_path, favor_db_path = favor_dir,
                 features = feats, update_gds = TRUE,
                 match_method = "exact", verbose = 0)

  g <- SeqArray::seqOpen(gds_path)
  on.exit(SeqArray::seqClose(g), add = TRUE)

  fa_node <- gdsfmt::index.gdsn(g, "annotation/info/FunctionalAnnotation")
  expect_true(length(gdsfmt::ls.gdsn(fa_node)) > 0)
  expect_setequal(gdsfmt::ls.gdsn(fa_node), feats)

  # Numeric + string sub-nodes both read via the scan path.
  cadd <- SeqArray::seqGetData(g, "annotation/info/FunctionalAnnotation/cadd_phred")
  expect_equal(as.numeric(cadd), fb$favor$cadd_phred)
  exonic <- SeqArray::seqGetData(
    g, "annotation/info/FunctionalAnnotation/genecode_comprehensive_exonic_category")
  expect_type(exonic, "character")
  expect_identical(as.character(exonic),
                   fb$favor$genecode_comprehensive_exonic_category)
})

# ------------------------------------------------------------------------------
# TEST: .write_agds (annotation-only path, non-GDS data.frame input)
# ------------------------------------------------------------------------------

test_that("annotate_favor(output_agds=) from a data.frame writes a sub-node folder", {
  skip_if_no_gds_tools()

  # Build the GDS only to harvest deterministic variant keys for the FAVOR DB;
  # the annotate call itself takes a plain data.frame (non-GDS) input, exercising
  # the .write_agds() path.
  helper_gds <- tempfile(fileext = ".gds")
  on.exit(unlink(helper_gds), add = TRUE)
  gf <- .build_tiny_gds(helper_gds)

  favor_dir <- tempfile("favordb_"); dir.create(favor_dir)
  on.exit(unlink(favor_dir, recursive = TRUE), add = TRUE)
  fb <- .build_tiny_favor_db(favor_dir, gf)
  feats <- .tiny_feats()

  varinfo <- paste(gf$chr, gf$pos, gf$ref, gf$alt, sep = "-")
  df_in <- data.frame(VarInfo = varinfo, stringsAsFactors = FALSE)

  out_agds <- tempfile(fileext = ".gds")
  on.exit(unlink(out_agds), add = TRUE)
  annotate_favor(variants = df_in, favor_db_path = favor_dir,
                 features = feats, output_agds = out_agds,
                 match_method = "exact", verbose = 0)

  g <- gdsfmt::openfn.gds(out_agds, readonly = TRUE)
  on.exit(gdsfmt::closefn.gds(g), add = TRUE)

  fa_node <- gdsfmt::index.gdsn(g, "annotation/info/FunctionalAnnotation")
  expect_true(length(gdsfmt::ls.gdsn(fa_node)) > 0)
  expect_setequal(gdsfmt::ls.gdsn(fa_node), feats)

  # String node preserved (read via gdsfmt directly).
  gene_cat <- gdsfmt::read.gdsn(
    gdsfmt::index.gdsn(fa_node, "genecode_comprehensive_category"))
  expect_type(gene_cat, "character")
  expect_identical(gene_cat, fb$favor$genecode_comprehensive_category)
})
