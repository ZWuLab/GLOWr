# Validation tests for coding variant set extraction against STAARpipeline reference
#
# These tests compare variant ID sets from extract_variant_set() with
# coding_variant_sets.rds, which contains STAARpipeline's variant selections
# for 77 chr22 genes x 7 coding categories.

# Resolve paths relative to project root (tests run from testthat dir)
proj_root <- normalizePath(file.path(testthat::test_path(), "..", "..", "..", ".."),
                           mustWork = FALSE)
gds_path <- file.path(proj_root, "data/local/large-data/c3-test/chr22_essentialdb.gds")
ref_path <- file.path(proj_root, "data/local/large-data/c3-test/staar_reference/coding_variant_sets.rds")

# Skip if real data files are not available
skip_if_not(file.exists(gds_path), "chr22 aGDS not available")
skip_if_not(file.exists(ref_path), "coding_variant_sets.rds reference not available")

test_that("extract_variant_set matches STAARpipeline for all 98 coding variant sets", {
  gds <- SeqArray::seqOpen(gds_path)
  on.exit(SeqArray::seqClose(gds))

  reference <- readRDS(ref_path)
  regions <- define_regions_gene(chr = 22)

  categories <- c("plof", "plof_ds", "missense", "disruptive_missense",
                   "synonymous", "ptv", "ptv_ds")

  mismatches <- character(0)
  total_compared <- 0

  for (gene_name in names(reference)) {
    region_row <- regions[regions$region_id == gene_name, ]
    if (nrow(region_row) == 0) next

    for (cat_name in categories) {
      ref_ids <- reference[[gene_name]][[cat_name]]
      if (is.null(ref_ids) || length(ref_ids) == 0) next

      total_compared <- total_compared + 1

      spec <- coding_filter(cat_name, rare_maf_cutoff = 0.5, min_variants = 1L)
      vset <- extract_variant_set(gds, region_row, spec, verbose = 0)

      if (is.null(vset)) {
        mismatches <- c(mismatches,
          paste0(gene_name, "/", cat_name, ": NULL (ref has ", length(ref_ids), ")"))
        next
      }

      our_ids <- vset$variant_info$variant_id
      if (!setequal(our_ids, ref_ids)) {
        mismatches <- c(mismatches,
          paste0(gene_name, "/", cat_name,
                 ": ours=", length(our_ids), " ref=", length(ref_ids)))
      }
    }
  }

  # Report details if there are mismatches (for debugging)
  if (length(mismatches) > 0) {
    message("Mismatches found:\n  ", paste(mismatches, collapse = "\n  "))
  }

  expect_equal(length(mismatches), 0,
    info = paste("Expected 0 mismatches across", total_compared, "comparisons"))
  expect_true(total_compared >= 90,
    info = paste("Expected >= 90 comparisons, got", total_compared))
})
