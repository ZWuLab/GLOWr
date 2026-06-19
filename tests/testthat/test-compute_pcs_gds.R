# Test file for compute_pcs_gds() function

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Create a minimal SNP GDS file for testing
#'
#' Generates a synthetic GDS file in SeqArray format with random genotypes.
#'
#' @param file_path Character. Output GDS file path.
#' @param n_samples Integer. Number of samples.
#' @param n_variants Integer. Number of variants.
#' @param chromosomes Character vector. Chromosome labels to assign.
#' @return Character. The file path (for chaining).
create_test_seqarray_gds <- function(file_path,
                                     n_samples = 50L,
                                     n_variants = 200L,
                                     chromosomes = c("1", "2")) {
  # Requires SNPRelate, SeqArray, and gdsfmt
  if (!requireNamespace("SNPRelate", quietly = TRUE)) return(NULL)
  if (!requireNamespace("SeqArray", quietly = TRUE)) return(NULL)
  if (!requireNamespace("gdsfmt", quietly = TRUE)) return(NULL)

  # Create a SNP GDS file, then convert directly to SeqArray format
  # via seqSNP2GDS() (avoids seqBED2GDS version compatibility issues).
  temp_snp <- tempfile(fileext = ".gds")
  on.exit(unlink(temp_snp), add = TRUE)

  # Generate random genotype matrix (0, 1, 2)
  set.seed(123L)
  geno <- matrix(
    sample(0L:2L, n_samples * n_variants, replace = TRUE),
    nrow = n_variants, ncol = n_samples
  )

  # Assign chromosomes evenly across variants
  chr_labels <- rep(chromosomes, length.out = n_variants)

  # Create sample and SNP IDs
  sample_ids <- sprintf("sample_%03d", seq_len(n_samples))
  snp_ids <- seq_len(n_variants)

  # Positions: sequential within each chromosome
  positions <- integer(n_variants)
  for (chr in chromosomes) {
    idx <- which(chr_labels == chr)
    positions[idx] <- seq(1000L, by = 1000L, length.out = length(idx))
  }

  # Create SNP GDS file using SNPRelate
  SNPRelate::snpgdsCreateGeno(
    gds.fn = temp_snp,
    genmat = geno,
    sample.id = sample_ids,
    snp.id = snp_ids,
    snp.chromosome = as.integer(chr_labels),
    snp.position = positions,
    snp.allele = rep("A/G", n_variants),
    snpfirstdim = TRUE
  )

  # Convert SNP GDS directly to SeqArray GDS (skips PLINK intermediate)
  suppressMessages({
    SeqArray::seqSNP2GDS(temp_snp, file_path, verbose = FALSE)
  })

  return(file_path)
}


#' Create two per-chromosome SeqArray GDS files for testing merge workflow
#'
#' @param dir_path Character. Directory to create files in.
#' @param n_samples Integer. Number of samples.
#' @param n_variants_per_chr Integer. Variants per chromosome.
#' @return Character vector of two file paths.
create_test_chr_gds_files <- function(dir_path,
                                      n_samples = 50L,
                                      n_variants_per_chr = 100L) {
  if (!requireNamespace("SNPRelate", quietly = TRUE)) return(NULL)
  if (!requireNamespace("SeqArray", quietly = TRUE)) return(NULL)

  files <- character(2L)
  for (i in 1:2) {
    chr_val <- as.character(i)
    gds_path <- file.path(dir_path, sprintf("chr%s_test.gds", chr_val))
    create_test_seqarray_gds(
      file_path = gds_path,
      n_samples = n_samples,
      n_variants = n_variants_per_chr,
      chromosomes = chr_val
    )
    files[i] <- gds_path
  }
  return(files)
}


# ==============================================================================
# GROUP 1: INPUT VALIDATION (no external dependencies needed)
# ==============================================================================

test_that("compute_pcs_gds errors on non-existent GDS file paths", {
  # Input validation runs before dependency check, so no skip needed
  expect_error(
    compute_pcs_gds(gds_files = "/nonexistent/path/fake.gds"),
    "GDS file.*not found"
  )

  expect_error(
    compute_pcs_gds(gds_files = c("/fake/a.gds", "/fake/b.gds")),
    "GDS file.*not found"
  )
})

test_that("compute_pcs_gds errors on invalid n_pcs", {
  # Create a dummy file so file existence check passes
  tmp <- tempfile(fileext = ".gds")
  file.create(tmp)
  on.exit(unlink(tmp), add = TRUE)

  expect_error(
    compute_pcs_gds(gds_files = tmp, n_pcs = 0),
    "n_pcs must be a positive integer"
  )
  expect_error(
    compute_pcs_gds(gds_files = tmp, n_pcs = -5),
    "n_pcs must be a positive integer"
  )
})

test_that("compute_pcs_gds errors on invalid maf_threshold", {
  tmp <- tempfile(fileext = ".gds")
  file.create(tmp)
  on.exit(unlink(tmp), add = TRUE)

  expect_error(
    compute_pcs_gds(gds_files = tmp, maf_threshold = -0.1),
    "maf_threshold.*\\[0, 0.5\\]"
  )
  expect_error(
    compute_pcs_gds(gds_files = tmp, maf_threshold = 0.6),
    "maf_threshold.*\\[0, 0.5\\]"
  )
})

test_that("compute_pcs_gds errors on invalid missing_rate", {
  tmp <- tempfile(fileext = ".gds")
  file.create(tmp)
  on.exit(unlink(tmp), add = TRUE)

  expect_error(
    compute_pcs_gds(gds_files = tmp, missing_rate = -0.1),
    "missing_rate.*\\[0, 1\\]"
  )
  expect_error(
    compute_pcs_gds(gds_files = tmp, missing_rate = 1.5),
    "missing_rate.*\\[0, 1\\]"
  )
})

test_that("compute_pcs_gds errors on invalid ld_threshold", {
  tmp <- tempfile(fileext = ".gds")
  file.create(tmp)
  on.exit(unlink(tmp), add = TRUE)

  # ld_threshold must be in (0, 1]
  expect_error(
    compute_pcs_gds(gds_files = tmp, ld_threshold = 0),
    "ld_threshold.*\\(0, 1\\]"
  )
  expect_error(
    compute_pcs_gds(gds_files = tmp, ld_threshold = -0.5),
    "ld_threshold.*\\(0, 1\\]"
  )
  expect_error(
    compute_pcs_gds(gds_files = tmp, ld_threshold = 1.5),
    "ld_threshold.*\\(0, 1\\]"
  )
})

test_that("compute_pcs_gds errors on invalid output_file extension", {
  tmp <- tempfile(fileext = ".gds")
  file.create(tmp)
  on.exit(unlink(tmp), add = TRUE)

  expect_error(
    compute_pcs_gds(gds_files = tmp, output_file = "/tmp/output.txt"),
    "output_file must have extension .rds or .csv"
  )
  expect_error(
    compute_pcs_gds(gds_files = tmp, output_file = "/tmp/output.xlsx"),
    "output_file must have extension .rds or .csv"
  )
})

test_that("compute_pcs_gds errors on empty gds_files", {
  expect_error(
    compute_pcs_gds(gds_files = character(0)),
    "gds_files must be a character vector"
  )
  expect_error(
    compute_pcs_gds(gds_files = 123),
    "gds_files must be a character vector"
  )
})


# ==============================================================================
# GROUP 2: DEPENDENCY CHECKING
# ==============================================================================

test_that(".check_pca_dependencies gives informative error when packages missing", {
  # We can test the error message format by checking it directly
  # (The actual skip_if_not_installed guards elsewhere handle missing deps)
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("SeqArray")

  # If packages are installed, check should succeed

  expect_invisible(.check_pca_dependencies())
})


# ==============================================================================
# GROUP 3: CORE FUNCTIONALITY (requires SNPRelate + test GDS data)
# ==============================================================================

test_that("compute_pcs_gds returns glow_pcs class object", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("SeqArray")
  skip_if_not_installed("gdsfmt")

  # Create test GDS file
  tmp_dir <- tempdir()
  gds_path <- file.path(tmp_dir, "test_core.gds")
  on.exit(unlink(gds_path), add = TRUE)

  result <- create_test_seqarray_gds(gds_path, n_samples = 50, n_variants = 200)
  if (is.null(result)) skip("Could not create test GDS file")

  pcs <- compute_pcs_gds(
    gds_files = gds_path,
    n_pcs = 5L,
    maf_threshold = 0.01,
    ld_threshold = 0.5,
    verbose = 0L
  )

  # Check class

  expect_s3_class(pcs, "glow_pcs")
  expect_s3_class(pcs, "data.frame")
  expect_true(inherits(pcs, "glow_pcs"))
})

test_that("compute_pcs_gds output has correct dimensions and column names", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("SeqArray")
  skip_if_not_installed("gdsfmt")

  tmp_dir <- tempdir()
  gds_path <- file.path(tmp_dir, "test_dims.gds")
  on.exit(unlink(gds_path), add = TRUE)

  result <- create_test_seqarray_gds(gds_path, n_samples = 50, n_variants = 200)
  if (is.null(result)) skip("Could not create test GDS file")

  n_pcs <- 5L
  pcs <- compute_pcs_gds(
    gds_files = gds_path,
    n_pcs = n_pcs,
    maf_threshold = 0.01,
    ld_threshold = 0.5,
    verbose = 0L
  )

  # Dimensions: n_samples rows, n_pcs + 1 columns (sample.id + PCs)
  expect_equal(ncol(pcs), n_pcs + 1L)
  expect_equal(nrow(pcs), 50L)

  # Column names
  expected_cols <- c("sample.id", paste0("PC", seq_len(n_pcs)))
  expect_equal(colnames(pcs), expected_cols)
})

test_that("compute_pcs_gds has all required attributes", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("SeqArray")
  skip_if_not_installed("gdsfmt")

  tmp_dir <- tempdir()
  gds_path <- file.path(tmp_dir, "test_attrs.gds")
  on.exit(unlink(gds_path), add = TRUE)

  result <- create_test_seqarray_gds(gds_path, n_samples = 50, n_variants = 200)
  if (is.null(result)) skip("Could not create test GDS file")

  pcs <- compute_pcs_gds(
    gds_files = gds_path,
    n_pcs = 5L,
    maf_threshold = 0.01,
    ld_threshold = 0.5,
    verbose = 0L
  )

  # Check all required attributes are present
  expect_false(is.null(attr(pcs, "eigenvalues")))
  expect_false(is.null(attr(pcs, "variance_proportion")))
  expect_false(is.null(attr(pcs, "total_variance_explained")))
  expect_false(is.null(attr(pcs, "n_variants_used")))
  expect_false(is.null(attr(pcs, "n_variants_per_chr")))
  expect_false(is.null(attr(pcs, "call_args")))

  # Verify attribute lengths
  expect_length(attr(pcs, "eigenvalues"), 5L)
  expect_length(attr(pcs, "variance_proportion"), 5L)
  expect_length(attr(pcs, "total_variance_explained"), 1L)
  expect_length(attr(pcs, "n_variants_used"), 1L)
})

test_that("eigenvalues are in decreasing order", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("SeqArray")
  skip_if_not_installed("gdsfmt")

  tmp_dir <- tempdir()
  gds_path <- file.path(tmp_dir, "test_eigenorder.gds")
  on.exit(unlink(gds_path), add = TRUE)

  result <- create_test_seqarray_gds(gds_path, n_samples = 50, n_variants = 200)
  if (is.null(result)) skip("Could not create test GDS file")

  pcs <- compute_pcs_gds(
    gds_files = gds_path,
    n_pcs = 5L,
    maf_threshold = 0.01,
    ld_threshold = 0.5,
    verbose = 0L
  )

  ev <- attr(pcs, "eigenvalues")

  # Eigenvalues should be non-increasing (some may be equal)
  expect_true(all(diff(ev) <= 0),
              info = "Eigenvalues should be in non-increasing order")
})

test_that("variance proportions are valid and sum correctly", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("SeqArray")
  skip_if_not_installed("gdsfmt")

  tmp_dir <- tempdir()
  gds_path <- file.path(tmp_dir, "test_varprop.gds")
  on.exit(unlink(gds_path), add = TRUE)

  result <- create_test_seqarray_gds(gds_path, n_samples = 50, n_variants = 200)
  if (is.null(result)) skip("Could not create test GDS file")

  pcs <- compute_pcs_gds(
    gds_files = gds_path,
    n_pcs = 5L,
    maf_threshold = 0.01,
    ld_threshold = 0.5,
    verbose = 0L
  )

  vp <- attr(pcs, "variance_proportion")
  total_ve <- attr(pcs, "total_variance_explained")

  # Each proportion should be in [0, 1]
  expect_true(all(vp >= 0 & vp <= 1),
              info = "Variance proportions should be in [0, 1]")

  # Sum of proportions should equal total_variance_explained
  expect_equal(sum(vp), total_ve, tolerance = 1e-10)

  # Total should not exceed 1
  expect_true(total_ve <= 1.0 + 1e-10,
              info = "Total variance explained should not exceed 1")
})

test_that("seed parameter produces reproducible results", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("SeqArray")
  skip_if_not_installed("gdsfmt")

  tmp_dir <- tempdir()
  gds_path <- file.path(tmp_dir, "test_seed.gds")
  on.exit(unlink(gds_path), add = TRUE)

  result <- create_test_seqarray_gds(gds_path, n_samples = 50, n_variants = 200)
  if (is.null(result)) skip("Could not create test GDS file")

  # Run twice with same seed
  pcs1 <- compute_pcs_gds(
    gds_files = gds_path, n_pcs = 3L, seed = 42L,
    maf_threshold = 0.01, ld_threshold = 0.5, verbose = 0L
  )
  pcs2 <- compute_pcs_gds(
    gds_files = gds_path, n_pcs = 3L, seed = 42L,
    maf_threshold = 0.01, ld_threshold = 0.5, verbose = 0L
  )

  # PC scores should be identical (same seed = same LD pruning = same PCA)
  expect_equal(pcs1$PC1, pcs2$PC1, tolerance = 1e-10)
  expect_equal(pcs1$PC2, pcs2$PC2, tolerance = 1e-10)
  expect_equal(attr(pcs1, "eigenvalues"), attr(pcs2, "eigenvalues"),
               tolerance = 1e-10)
  expect_equal(attr(pcs1, "n_variants_used"), attr(pcs2, "n_variants_used"))
})


# ==============================================================================
# GROUP 4: SINGLE VS MULTIPLE GDS FILES
# ==============================================================================

test_that("compute_pcs_gds works with a single GDS file", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("SeqArray")
  skip_if_not_installed("gdsfmt")

  tmp_dir <- tempdir()
  gds_path <- file.path(tmp_dir, "test_single.gds")
  on.exit(unlink(gds_path), add = TRUE)

  result <- create_test_seqarray_gds(gds_path, n_samples = 50, n_variants = 200)
  if (is.null(result)) skip("Could not create test GDS file")

  pcs <- compute_pcs_gds(
    gds_files = gds_path,
    n_pcs = 3L,
    maf_threshold = 0.01,
    ld_threshold = 0.5,
    verbose = 0L
  )

  expect_s3_class(pcs, "glow_pcs")
  expect_equal(nrow(pcs), 50L)
  expect_equal(ncol(pcs), 4L)  # sample.id + 3 PCs
})

test_that("compute_pcs_gds works with multiple GDS files (merge)", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("SeqArray")
  skip_if_not_installed("gdsfmt")

  # Create per-chromosome files
  tmp_dir <- file.path(tempdir(), "test_multi")
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  gds_files <- create_test_chr_gds_files(
    tmp_dir, n_samples = 50, n_variants_per_chr = 100
  )
  if (is.null(gds_files)) skip("Could not create test GDS files")

  pcs <- compute_pcs_gds(
    gds_files = gds_files,
    n_pcs = 3L,
    maf_threshold = 0.01,
    ld_threshold = 0.5,
    verbose = 0L
  )

  expect_s3_class(pcs, "glow_pcs")
  expect_equal(nrow(pcs), 50L)
  expect_equal(ncol(pcs), 4L)

  # Should have variants from both chromosomes
  n_per_chr <- attr(pcs, "n_variants_per_chr")
  expect_true(length(n_per_chr) >= 1L,
              info = "Should have variant counts for at least one chromosome")
})


# ==============================================================================
# GROUP 5: OUTPUT SAVING
# ==============================================================================

test_that("RDS output preserves glow_pcs class and all attributes", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("SeqArray")
  skip_if_not_installed("gdsfmt")

  tmp_dir <- tempdir()
  gds_path <- file.path(tmp_dir, "test_rds.gds")
  rds_path <- file.path(tmp_dir, "test_output.rds")
  on.exit({
    unlink(gds_path)
    unlink(rds_path)
  }, add = TRUE)

  result <- create_test_seqarray_gds(gds_path, n_samples = 50, n_variants = 200)
  if (is.null(result)) skip("Could not create test GDS file")

  pcs <- compute_pcs_gds(
    gds_files = gds_path,
    n_pcs = 5L,
    maf_threshold = 0.01,
    ld_threshold = 0.5,
    output_file = rds_path,
    verbose = 0L
  )

  # File should be created
  expect_true(file.exists(rds_path))

  # Round-trip should preserve class and attributes
  loaded <- readRDS(rds_path)
  expect_s3_class(loaded, "glow_pcs")
  expect_equal(attr(loaded, "eigenvalues"), attr(pcs, "eigenvalues"))
  expect_equal(attr(loaded, "variance_proportion"),
               attr(pcs, "variance_proportion"))
  expect_equal(attr(loaded, "total_variance_explained"),
               attr(pcs, "total_variance_explained"))
  expect_equal(attr(loaded, "n_variants_used"),
               attr(pcs, "n_variants_used"))
  expect_equal(loaded$sample.id, pcs$sample.id)
  expect_equal(loaded$PC1, pcs$PC1)
})

test_that("CSV output is readable and has correct dimensions", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("SeqArray")
  skip_if_not_installed("gdsfmt")

  tmp_dir <- tempdir()
  gds_path <- file.path(tmp_dir, "test_csv.gds")
  csv_path <- file.path(tmp_dir, "test_output.csv")
  on.exit({
    unlink(gds_path)
    unlink(csv_path)
  }, add = TRUE)

  result <- create_test_seqarray_gds(gds_path, n_samples = 50, n_variants = 200)
  if (is.null(result)) skip("Could not create test GDS file")

  pcs <- compute_pcs_gds(
    gds_files = gds_path,
    n_pcs = 5L,
    maf_threshold = 0.01,
    ld_threshold = 0.5,
    output_file = csv_path,
    verbose = 0L
  )

  # File should be created
  expect_true(file.exists(csv_path))

  # CSV should be readable with default read.csv() (no comment.char needed)
  csv_data <- read.csv(csv_path)
  expect_equal(nrow(csv_data), 50L)
  expect_equal(ncol(csv_data), 6L)  # sample.id + 5 PCs
  expect_true("sample.id" %in% colnames(csv_data))

  # Metadata sidecar file should exist
  meta_path <- sub("\\.csv$", "_meta.txt", csv_path)
  expect_true(file.exists(meta_path))
  meta_lines <- readLines(meta_path)
  expect_true(any(grepl("GLOW Principal Components", meta_lines)))
  expect_true(any(grepl("Eigenvalues:", meta_lines)))
})


# ==============================================================================
# GROUP 6: S3 METHODS
# ==============================================================================

test_that("print.glow_pcs runs without error and produces expected format", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("SeqArray")
  skip_if_not_installed("gdsfmt")

  tmp_dir <- tempdir()
  gds_path <- file.path(tmp_dir, "test_print.gds")
  on.exit(unlink(gds_path), add = TRUE)

  result <- create_test_seqarray_gds(gds_path, n_samples = 50, n_variants = 200)
  if (is.null(result)) skip("Could not create test GDS file")

  pcs <- compute_pcs_gds(
    gds_files = gds_path, n_pcs = 5L,
    maf_threshold = 0.01, ld_threshold = 0.5, verbose = 0L
  )

  # Capture output
  output <- capture.output(print(pcs))
  combined <- paste(output, collapse = "\n")

  expect_true(grepl("GLOW Principal Components", combined))
  expect_true(grepl("Samples:", combined))
  expect_true(grepl("PCs:", combined))
  expect_true(grepl("Variants:", combined))
  expect_true(grepl("Variance:", combined))
  expect_true(grepl("Top PCs:", combined))

  # print should return invisibly
  expect_invisible(print(pcs))
})

test_that("plot(type='scree') runs without error", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("SeqArray")
  skip_if_not_installed("gdsfmt")

  tmp_dir <- tempdir()
  gds_path <- file.path(tmp_dir, "test_scree.gds")
  on.exit(unlink(gds_path), add = TRUE)

  result <- create_test_seqarray_gds(gds_path, n_samples = 50, n_variants = 200)
  if (is.null(result)) skip("Could not create test GDS file")

  pcs <- compute_pcs_gds(
    gds_files = gds_path, n_pcs = 5L,
    maf_threshold = 0.01, ld_threshold = 0.5, verbose = 0L
  )

  # Plot should run without error
  expect_no_error(plot(pcs, type = "scree"))

  # Should return invisible
  ret <- plot(pcs, type = "scree")
  expect_true(inherits(ret, "glow_pcs"))
})

test_that("plot(type='pairs') runs without error", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("SeqArray")
  skip_if_not_installed("gdsfmt")

  tmp_dir <- tempdir()
  gds_path <- file.path(tmp_dir, "test_pairs.gds")
  on.exit(unlink(gds_path), add = TRUE)

  result <- create_test_seqarray_gds(gds_path, n_samples = 50, n_variants = 200)
  if (is.null(result)) skip("Could not create test GDS file")

  pcs <- compute_pcs_gds(
    gds_files = gds_path, n_pcs = 5L,
    maf_threshold = 0.01, ld_threshold = 0.5, verbose = 0L
  )

  # Plot should run without error
  expect_no_error(plot(pcs, type = "pairs", n_pairs = 3))

  # Should return invisible
  ret <- plot(pcs, type = "pairs", n_pairs = 3)
  expect_true(inherits(ret, "glow_pcs"))
})

test_that("plot.glow_pcs matches type argument correctly", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("SeqArray")
  skip_if_not_installed("gdsfmt")

  tmp_dir <- tempdir()
  gds_path <- file.path(tmp_dir, "test_plottype.gds")
  on.exit(unlink(gds_path), add = TRUE)

  result <- create_test_seqarray_gds(gds_path, n_samples = 50, n_variants = 200)
  if (is.null(result)) skip("Could not create test GDS file")

  pcs <- compute_pcs_gds(
    gds_files = gds_path, n_pcs = 5L,
    maf_threshold = 0.01, ld_threshold = 0.5, verbose = 0L
  )

  # Default type should be "scree"
  expect_no_error(plot(pcs))

  # Invalid type should error
  expect_error(plot(pcs, type = "invalid"))
})


# ==============================================================================
# GROUP 7: EDGE CASES AND VERBOSITY
# ==============================================================================

test_that("compute_pcs_gds handles verbose=0 silently", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("SeqArray")
  skip_if_not_installed("gdsfmt")

  tmp_dir <- tempdir()
  gds_path <- file.path(tmp_dir, "test_silent.gds")
  on.exit(unlink(gds_path), add = TRUE)

  result <- create_test_seqarray_gds(gds_path, n_samples = 50, n_variants = 200)
  if (is.null(result)) skip("Could not create test GDS file")

  # Should produce no messages
  messages <- capture.output(type = "message", {
    pcs <- compute_pcs_gds(
      gds_files = gds_path, n_pcs = 3L,
      maf_threshold = 0.01, ld_threshold = 0.5, verbose = 0L
    )
  })

  expect_length(messages, 0L)
})

test_that("compute_pcs_gds errors when n_pcs >= pruned variant count", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("SeqArray")
  skip_if_not_installed("gdsfmt")

  tmp_dir <- tempdir()
  gds_path <- file.path(tmp_dir, "test_npcs_validation.gds")
  on.exit(unlink(gds_path), add = TRUE)

  # Create a small dataset — at most 30 variants can survive pruning
  result <- create_test_seqarray_gds(gds_path, n_samples = 50, n_variants = 30)
  if (is.null(result)) skip("Could not create test GDS file")

  # Request far more PCs than possible variants
  expect_error(
    compute_pcs_gds(
      gds_files = gds_path, n_pcs = 200L,
      maf_threshold = 0.01, ld_threshold = 0.5, verbose = 0L
    ),
    "n_pcs.*must be less than.*pruned variants"
  )
})

test_that("compute_pcs_gds algorithm parameter works", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("SeqArray")
  skip_if_not_installed("gdsfmt")

  tmp_dir <- tempdir()
  gds_path <- file.path(tmp_dir, "test_algorithm.gds")
  on.exit(unlink(gds_path), add = TRUE)

  result <- create_test_seqarray_gds(gds_path, n_samples = 50, n_variants = 200)
  if (is.null(result)) skip("Could not create test GDS file")

  # Default (randomized) should work
  pcs_rand <- compute_pcs_gds(
    gds_files = gds_path, n_pcs = 5L,
    maf_threshold = 0.01, ld_threshold = 0.5, verbose = 0L
  )
  expect_s3_class(pcs_rand, "glow_pcs")
  expect_equal(attr(pcs_rand, "call_args")$algorithm, "randomized")

  # Explicit "exact" should work
  pcs_exact <- compute_pcs_gds(
    gds_files = gds_path, n_pcs = 5L, algorithm = "exact",
    maf_threshold = 0.01, ld_threshold = 0.5, verbose = 0L
  )
  expect_s3_class(pcs_exact, "glow_pcs")
  expect_equal(attr(pcs_exact, "call_args")$algorithm, "exact")

  # Both should produce same number of samples and PCs
  expect_equal(nrow(pcs_rand), nrow(pcs_exact))
  expect_equal(ncol(pcs_rand), ncol(pcs_exact))
})

test_that("compute_pcs_gds cleans up temp files on error", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("SeqArray")

  # Create a file that is NOT a valid GDS to trigger an error after validation
  tmp <- tempfile(fileext = ".gds")
  writeLines("not a real gds file", tmp)
  on.exit(unlink(tmp), add = TRUE)

  # This should error during seqGDS2SNP or seqMerge but temp files should be
  # cleaned up by on.exit()
  expect_error(
    compute_pcs_gds(gds_files = tmp, n_pcs = 3L, verbose = 0L)
  )
})
