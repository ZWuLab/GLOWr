# test-train_PI_models.R
# Unit tests for PI model training functions

# ==============================================================================
# Setup: Create test data paths
# ==============================================================================

# Use sample data from data/local/large-data/test/
test_case_csv <- file.path(
  testthat::test_path(), "..", "..", "..", "..",
  "data/local/large-data/test/sample_cases_50.csv"
)
test_ctrl_csv <- file.path(
  testthat::test_path(), "..", "..", "..", "..",
  "data/local/large-data/test/sample_controls_chr21_500.csv"
)

# Skip tests if test data doesn't exist
skip_if_no_test_data <- function() {
  if (!file.exists(test_case_csv) || !file.exists(test_ctrl_csv)) {
    skip("Test data files not available")
  }
}

# ==============================================================================
# Tests for .default_PI_features()
# ==============================================================================

test_that(".default_PI_features returns correct features", {
  features <- .default_PI_features()

  expect_type(features, "character")
  expect_length(features, 16)

  # Check specific features
  expect_true("apc_conservation_v2" %in% features)
  expect_true("apc_micro_rna" %in% features)
  expect_true("apc_proximity_to_coding_v2" %in% features)
  expect_true("cadd_phred" %in% features)
  expect_true("linsight" %in% features)
  expect_true("fathmm_xf" %in% features)
})


# ==============================================================================
# Tests for impute_na_median()
# ==============================================================================

test_that("impute_na_median handles matrix with NAs", {
  # Matrix by column: col1 = c(1, 2, NA), col2 = c(4, NA, 6)
  mat <- matrix(c(1, 2, NA, 4, NA, 6), ncol = 2)
  result <- impute_na_median(mat)

  expect_equal(sum(is.na(result)), 0)
  expect_equal(result[3, 1], 1.5)  # median of c(1, 2) = 1.5
  expect_equal(result[2, 2], 5)    # median of c(4, 6) = 5
})

test_that("impute_na_median handles matrix without NAs", {
  mat <- matrix(1:6, ncol = 2)
  result <- impute_na_median(mat)

  expect_equal(result, mat)
})

test_that("impute_na_median handles data frame input", {
  df <- data.frame(a = c(1, 2, NA), b = c(NA, 5, 6))
  result <- impute_na_median(df)

  expect_true(is.matrix(result))
  expect_equal(sum(is.na(result)), 0)
})

test_that("impute_na_median handles column with all NAs", {
  mat <- matrix(c(1, 2, 3, NA, NA, NA), ncol = 2)
  result <- impute_na_median(mat)

  # Column 2 is all NA, so should remain NA
  expect_true(all(is.na(result[, 2])))
  expect_equal(sum(is.na(result[, 1])), 0)
})


# ==============================================================================
# Tests for load_case_annotations()
# ==============================================================================

test_that("load_case_annotations loads CSV correctly", {
  skip_if_no_test_data()

  result <- load_case_annotations(test_case_csv)

  expect_true(is.matrix(result))
  expect_equal(ncol(result), length(.default_PI_features()))
  expect_true(nrow(result) > 0)

  # Check column names
  expect_equal(colnames(result), .default_PI_features())

  # Check row names (VarInfo)
  expect_true(all(grepl("-", rownames(result))))  # VarInfo format: CHR-POS-REF-ALT
})

test_that("load_case_annotations errors on missing file", {
  expect_error(
    load_case_annotations("nonexistent.csv"),
    "File not found"
  )
})

test_that("load_case_annotations allows custom features", {
  skip_if_no_test_data()

  custom_features <- c("cadd_phred", "linsight")
  result <- load_case_annotations(test_case_csv, features = custom_features)

  expect_equal(ncol(result), 2)
  expect_equal(colnames(result), custom_features)
})

test_that("load_case_annotations errors on missing feature columns", {
  skip_if_no_test_data()

  expect_error(
    load_case_annotations(test_case_csv, features = c("nonexistent_col")),
    "Missing annotation columns"
  )
})


# ==============================================================================
# Tests for load_control_annotations()
# ==============================================================================

test_that("load_control_annotations loads single CSV correctly", {
  skip_if_no_test_data()

  result <- load_control_annotations(test_ctrl_csv)

  expect_true(is.matrix(result))
  expect_equal(ncol(result), length(.default_PI_features()))
  expect_true(nrow(result) > 0)
})

test_that("load_control_annotations errors on nonexistent source", {
  expect_error(
    load_control_annotations("nonexistent_path"),
    "Source not found"
  )
})

test_that("load_control_annotations respects max_controls", {
  skip_if_no_test_data()

  set.seed(123)  # For reproducible sampling
  result <- load_control_annotations(test_ctrl_csv, max_controls = 100)

  expect_equal(nrow(result), 100)
  expect_equal(ncol(result), length(.default_PI_features()))
})

test_that("load_control_annotations loads from directory", {
  skip_if_no_test_data()

  # Create a temp directory with two CSV files
  temp_dir <- tempfile("chr_dir_")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))

  # Copy test file as chr21 and chr22
  file.copy(test_ctrl_csv, file.path(temp_dir, "chr21_test.csv"))
  file.copy(test_ctrl_csv, file.path(temp_dir, "chr22_test.csv"))

  result <- load_control_annotations(temp_dir, chromosomes = c(21, 22))

  expect_true(is.matrix(result))
  # Should have 2x the rows (loaded both files)
  expect_equal(nrow(result), 2 * 500)
})

test_that("load_control_annotations filters chromosomes", {
  skip_if_no_test_data()

  temp_dir <- tempfile("chr_dir_")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))

  file.copy(test_ctrl_csv, file.path(temp_dir, "chr21_test.csv"))
  file.copy(test_ctrl_csv, file.path(temp_dir, "chr22_test.csv"))

  # Only load chr21
  result <- load_control_annotations(temp_dir, chromosomes = 21)

  expect_equal(nrow(result), 500)  # Only one file loaded
})


# ==============================================================================
# Tests for GDS Support
# ==============================================================================

# Path to GDS test file
test_ctrl_gds <- file.path(
  testthat::test_path(), "..", "..", "..", "..",
  "data/local/large-data/test/my_controls_annotated.gds"
)

skip_if_no_gds_test_data <- function() {
  if (!file.exists(test_ctrl_gds)) {
    skip("GDS test data files not available")
  }
}

# The control GDS fixture (my_controls_annotated.gds) is a git-ignored real-data
# extract built by annotate_favor_batch.R before the default PI feature set grew
# to 16 (.default_PI_features()); it carries only the original 11 FAVOR columns
# and is not regenerated when the default changes. These GDS tests exercise GDS
# I/O mechanics (matrix-format reading, variant filtering, directory globbing,
# proportional sampling), which are feature-set-agnostic, so they request
# exactly the columns the fixture provides -- read from the file so they upgrade
# automatically if it is ever re-annotated with the full default set. The
# 16-feature default itself is exercised by the CSV-based tests above
# (sample_cases_50.csv / sample_controls_chr21_500.csv carry all 16 columns).
.ctrl_gds_features <- function() {
  gds <- SeqArray::seqOpen(test_ctrl_gds, readonly = TRUE)
  on.exit(SeqArray::seqClose(gds))
  node <- gdsfmt::index.gdsn(gds, "annotation/info/FunctionalAnnotation")
  gdsfmt::get.attr.gdsn(node)$feature_names
}

test_that(".load_control_gds loads GDS file correctly", {
  skip_if_no_gds_test_data()

  features <- .ctrl_gds_features()
  result <- GLOWr:::.load_control_gds(test_ctrl_gds, features)

  expect_true(is.matrix(result))
  expect_equal(ncol(result), length(features))
  expect_equal(nrow(result), 800)  # Test file has 800 variants
  expect_equal(colnames(result), features)
  expect_true(all(grepl("-", rownames(result))))  # VarInfo format
})

test_that(".load_control_gds respects variant_indices", {
  skip_if_no_gds_test_data()

  features <- .ctrl_gds_features()
  indices <- c(1, 10, 50, 100)
  result <- GLOWr:::.load_control_gds(
    test_ctrl_gds,
    features,
    variant_indices = indices
  )

  expect_equal(nrow(result), 4)

  # Verify correct rows selected
  full_result <- GLOWr:::.load_control_gds(test_ctrl_gds, features)
  expect_equal(result[1, ], full_result[1, ])
  expect_equal(result[2, ], full_result[10, ])
})

test_that("load_control_annotations loads GDS from directory", {
  skip_if_no_gds_test_data()

  # Create temp dir with GDS files
  temp_dir <- tempfile("gds_dir_")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))

  # Symlink as chr21 and chr22
  file.symlink(normalizePath(test_ctrl_gds), file.path(temp_dir, "chr21_test.gds"))
  file.symlink(normalizePath(test_ctrl_gds), file.path(temp_dir, "chr22_test.gds"))

  # Resolve fixture features eagerly: the symlinks alias the same physical file,
  # and gdsfmt refuses to open it twice, so it must be read before the load.
  features <- .ctrl_gds_features()
  result <- load_control_annotations(
    temp_dir, chromosomes = c(21, 22), features = features
  )

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2 * 800)  # 800 variants per file × 2 files
})

test_that("load_control_annotations respects max_controls for GDS", {
  skip_if_no_gds_test_data()

  # Create temp dir with GDS files
  temp_dir <- tempfile("gds_dir_")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))

  file.symlink(normalizePath(test_ctrl_gds), file.path(temp_dir, "chr21_test.gds"))
  file.symlink(normalizePath(test_ctrl_gds), file.path(temp_dir, "chr22_test.gds"))

  set.seed(123)
  features <- .ctrl_gds_features()
  result <- load_control_annotations(
    temp_dir,
    chromosomes = c(21, 22),
    max_controls = 500,
    features = features
  )

  expect_equal(nrow(result), 500)
  expect_equal(ncol(result), length(features))
})


# ==============================================================================
# .load_control_gds: reads BOTH the sub-node folder format (written by
# annotate_favor) AND the legacy single-matrix format. These build tiny aGDS
# fixtures INLINE (no real-data dependency), so they always run.
# ==============================================================================

# Build a tiny SeqArray aGDS whose FunctionalAnnotation node is written in either
# the "folder" (one sub-node per feature) or "matrix" layout, with identical
# values, so the two .load_control_gds branches can be compared head-to-head.
.build_loader_fixture <- function(path, features, values, layout) {
  n <- nrow(values)
  tmp_snp <- tempfile(fileext = ".gds")
  on.exit(unlink(tmp_snp), add = TRUE)
  geno <- matrix(0L, nrow = n, ncol = 6L)
  SNPRelate::snpgdsCreateGeno(
    gds.fn = tmp_snp, genmat = geno,
    sample.id = sprintf("s%02d", 1:6), snp.id = seq_len(n),
    snp.chromosome = rep(22L, n),
    snp.position = as.integer(seq(1000, by = 100, length.out = n)),
    snp.allele = rep("A/G", n), snpfirstdim = TRUE
  )
  suppressMessages(SeqArray::seqSNP2GDS(tmp_snp, path, verbose = FALSE))

  g <- gdsfmt::openfn.gds(path, readonly = FALSE)
  on.exit(gdsfmt::closefn.gds(g), add = TRUE)
  info <- gdsfmt::index.gdsn(g, "annotation/info")

  if (layout == "folder") {
    fa <- gdsfmt::addfolder.gdsn(info, "FunctionalAnnotation")
    for (f in features) {
      gdsfmt::add.gdsn(fa, f, val = values[[f]],
                       compress = "LZMA_RA", closezip = TRUE)
    }
    gdsfmt::put.attr.gdsn(fa, "feature_names", features)
  } else if (layout == "matrix") {
    mat <- as.matrix(values[, features, drop = FALSE])
    fa <- gdsfmt::add.gdsn(info, "FunctionalAnnotation", val = mat,
                           compress = "LZMA_RA", closezip = TRUE)
    gdsfmt::put.attr.gdsn(fa, "feature_names", features)
  } else {
    stop("unknown layout")
  }
  invisible(path)
}

test_that(".load_control_gds returns identical matrices from folder and matrix aGDS", {
  for (pkg in c("SNPRelate", "SeqArray", "gdsfmt")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      skip(paste0("Package '", pkg, "' not available"))
    }
  }

  set.seed(7)
  features <- c("cadd_phred", "linsight", "apc_conservation_v2")
  values <- data.frame(
    cadd_phred         = round(runif(8, 0, 40), 4),
    linsight           = round(runif(8, 0, 1), 4),
    apc_conservation_v2 = round(runif(8, 0, 50), 4)
  )

  folder_gds <- tempfile(fileext = ".gds")
  matrix_gds <- tempfile(fileext = ".gds")
  on.exit(unlink(c(folder_gds, matrix_gds)), add = TRUE)
  .build_loader_fixture(folder_gds, features, values, "folder")
  .build_loader_fixture(matrix_gds, features, values, "matrix")

  res_folder <- GLOWr:::.load_control_gds(folder_gds, features)
  res_matrix <- GLOWr:::.load_control_gds(matrix_gds, features)

  expect_true(is.matrix(res_folder))
  expect_true(is.matrix(res_matrix))
  expect_equal(colnames(res_folder), features)
  expect_equal(colnames(res_matrix), features)
  # Same values from both layouts (rownames are VarInfo, identical GDS variants).
  expect_equal(res_folder, res_matrix)
  # And the loaded values equal the source.
  expect_equivalent(res_folder[, "cadd_phred"], values$cadd_phred)
})

test_that(".load_control_gds folder branch respects variant_indices", {
  for (pkg in c("SNPRelate", "SeqArray", "gdsfmt")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      skip(paste0("Package '", pkg, "' not available"))
    }
  }

  set.seed(8)
  features <- c("cadd_phred", "linsight")
  values <- data.frame(
    cadd_phred = round(runif(10, 0, 40), 4),
    linsight   = round(runif(10, 0, 1), 4)
  )
  folder_gds <- tempfile(fileext = ".gds")
  on.exit(unlink(folder_gds), add = TRUE)
  .build_loader_fixture(folder_gds, features, values, "folder")

  full <- GLOWr:::.load_control_gds(folder_gds, features)
  sub  <- GLOWr:::.load_control_gds(folder_gds, features,
                                    variant_indices = c(1L, 4L, 9L))
  expect_equal(nrow(sub), 3L)
  expect_equal(sub[1, ], full[1, ])
  expect_equal(sub[2, ], full[4, ])
  expect_equal(sub[3, ], full[9, ])
})

test_that(".load_control_gds errors on a missing feature in folder format", {
  for (pkg in c("SNPRelate", "SeqArray", "gdsfmt")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      skip(paste0("Package '", pkg, "' not available"))
    }
  }
  set.seed(9)
  features <- c("cadd_phred", "linsight")
  values <- data.frame(cadd_phred = runif(5), linsight = runif(5))
  folder_gds <- tempfile(fileext = ".gds")
  on.exit(unlink(folder_gds), add = TRUE)
  .build_loader_fixture(folder_gds, features, values, "folder")

  expect_error(
    GLOWr:::.load_control_gds(folder_gds, c("cadd_phred", "not_a_feature")),
    "Missing annotation columns"
  )
})


# ==============================================================================
# Tests for train_PI_models()
# ==============================================================================

test_that("train_PI_models produces correct number of models", {
  skip_if_no_test_data()

  output_dir <- tempfile("piModels_")

  result <- train_PI_models(
    case_csv = test_case_csv,
    control_source = test_ctrl_csv,
    output_dir = output_dir,
    n_models = 3,
    controls_per_model = 50,
    model_type = "GLM",
    random_seed = 42,
    verbose = 0
  )

  expect_length(result$models, 3)
  expect_equal(length(list.files(output_dir, "\\.rds$")), 3)

  # Clean up
  unlink(output_dir, recursive = TRUE)
})

test_that("train_PI_models returns correct metadata", {
  skip_if_no_test_data()

  output_dir <- tempfile("piModels_")

  result <- train_PI_models(
    case_csv = test_case_csv,
    control_source = test_ctrl_csv,
    output_dir = output_dir,
    n_models = 2,
    controls_per_model = 50,
    model_type = "GLM",
    random_seed = 123,
    verbose = 0
  )

  expect_equal(result$metadata$n_cases, 50)
  expect_equal(result$metadata$n_controls, 500)
  expect_equal(result$metadata$n_models, 2)
  expect_equal(result$metadata$controls_per_model, 50)
  expect_equal(result$metadata$model_type, "GLM")
  expect_equal(result$metadata$random_seed, 123)

  unlink(output_dir, recursive = TRUE)
})

test_that("train_PI_models is reproducible with seed", {
  skip_if_no_test_data()

  output_dir1 <- tempfile("piModels_")
  output_dir2 <- tempfile("piModels_")

  result1 <- train_PI_models(
    case_csv = test_case_csv,
    control_source = test_ctrl_csv,
    output_dir = output_dir1,
    n_models = 2,
    controls_per_model = 50,
    random_seed = 42,
    verbose = 0
  )

  result2 <- train_PI_models(
    case_csv = test_case_csv,
    control_source = test_ctrl_csv,
    output_dir = output_dir2,
    n_models = 2,
    controls_per_model = 50,
    random_seed = 42,
    verbose = 0
  )

  # Coefficients should be identical
  expect_equal(
    coef(result1$models[[1]]),
    coef(result2$models[[1]])
  )

  unlink(output_dir1, recursive = TRUE)
  unlink(output_dir2, recursive = TRUE)
})

test_that("train_PI_models supports LASSO model type", {
  skip_if_no_test_data()

  output_dir <- tempfile("piModels_")

  result <- train_PI_models(
    case_csv = test_case_csv,
    control_source = test_ctrl_csv,
    output_dir = output_dir,
    n_models = 2,
    controls_per_model = 50,
    model_type = "LASSO",
    random_seed = 42,
    verbose = 0
  )

  # LASSO models should be glmnet objects
  expect_true(inherits(result$models[[1]], "glmnet") ||
              inherits(result$models[[1]], "cv.glmnet"))

  unlink(output_dir, recursive = TRUE)
})

test_that("train_PI_models errors on invalid model_type", {
  skip_if_no_test_data()

  expect_error(
    train_PI_models(
      case_csv = test_case_csv,
      control_source = test_ctrl_csv,
      output_dir = tempfile(),
      model_type = "INVALID"
    ),
    "model_type must be"
  )
})

test_that("train_PI_models errors on insufficient controls", {
  skip_if_no_test_data()

  expect_error(
    train_PI_models(
      case_csv = test_case_csv,
      control_source = test_ctrl_csv,
      output_dir = tempfile(),
      controls_per_model = 10000,  # More than available
      verbose = 0
    ),
    "exceeds available controls"
  )
})

test_that("train_PI_models auto-calculates controls_per_model from n_cases", {
  skip_if_no_test_data()

  output_dir <- tempfile("piModels_")

  # Default: controls_per_model = n_cases * 1 = 50
  result <- train_PI_models(
    case_csv = test_case_csv,
    control_source = test_ctrl_csv,
    output_dir = output_dir,
    n_models = 2,
    # controls_per_model = NULL (default)
    model_type = "GLM",
    random_seed = 42,
    verbose = 0
  )

  # Should auto-calculate controls_per_model = n_cases (50)
  expect_equal(result$metadata$n_cases, 50)
  expect_equal(result$metadata$controls_per_model, 50)
  expect_equal(result$metadata$controls_multiplier, 1)

  unlink(output_dir, recursive = TRUE)
})

test_that("train_PI_models respects controls_multiplier", {
  skip_if_no_test_data()

  output_dir <- tempfile("piModels_")

  # controls_per_model = n_cases * 2 = 100
  result <- train_PI_models(
    case_csv = test_case_csv,
    control_source = test_ctrl_csv,
    output_dir = output_dir,
    n_models = 2,
    controls_multiplier = 2,  # 2x cases
    model_type = "GLM",
    random_seed = 42,
    verbose = 0
  )

  expect_equal(result$metadata$n_cases, 50)
  expect_equal(result$metadata$controls_per_model, 100)  # 50 * 2
  expect_equal(result$metadata$controls_multiplier, 2)

  unlink(output_dir, recursive = TRUE)
})

test_that("train_PI_models explicit controls_per_model overrides multiplier", {
  skip_if_no_test_data()

  output_dir <- tempfile("piModels_")

  # Explicit controls_per_model should override multiplier
  result <- train_PI_models(
    case_csv = test_case_csv,
    control_source = test_ctrl_csv,
    output_dir = output_dir,
    n_models = 2,
    controls_per_model = 75,  # Explicit value
    controls_multiplier = 2,  # Should be ignored
    model_type = "GLM",
    random_seed = 42,
    verbose = 0
  )

  expect_equal(result$metadata$controls_per_model, 75)

  unlink(output_dir, recursive = TRUE)
})

test_that("train_PI_models errors on invalid controls_multiplier", {
  skip_if_no_test_data()

  expect_error(
    train_PI_models(
      case_csv = test_case_csv,
      control_source = test_ctrl_csv,
      output_dir = tempfile(),
      controls_multiplier = -1,
      verbose = 0
    ),
    "controls_multiplier must be a positive number"
  )
})


# ==============================================================================
# Tests for save_PI_models() and load_PI_models()
# ==============================================================================

test_that("save_PI_models creates correct files", {
  skip_if_no_test_data()

  # Create dummy models
  models <- list(
    lm(mpg ~ wt, data = mtcars),
    lm(mpg ~ hp, data = mtcars)
  )

  output_dir <- tempfile("piModels_")
  paths <- save_PI_models(models, output_dir)

  expect_length(paths, 2)
  expect_true(file.exists(file.path(output_dir, "model_1.rds")))
  expect_true(file.exists(file.path(output_dir, "model_2.rds")))

  unlink(output_dir, recursive = TRUE)
})

test_that("load_PI_models loads models correctly", {
  skip_if_no_test_data()

  # Create and save GLM models (load_PI_models requires glm, not lm)
  models <- list(
    glm(vs ~ wt, data = mtcars, family = binomial()),
    glm(vs ~ hp, data = mtcars, family = binomial())
  )

  output_dir <- tempfile("piModels_")
  save_PI_models(models, output_dir)

  # Load models (returns list with $models, $model_type, etc.)
  loaded <- load_PI_models(output_dir)

  expect_equal(loaded$n_models, 2)
  expect_equal(loaded$model_type, "GLM")
  expect_equal(coef(loaded$models[[1]]), coef(models[[1]]))

  unlink(output_dir, recursive = TRUE)
})

test_that("load_PI_models errors on nonexistent directory", {
  expect_error(
    load_PI_models("nonexistent_dir"),
    "Model directory does not exist"
  )
})

test_that("load_PI_models errors on empty directory", {
  empty_dir <- tempfile("empty_")
  dir.create(empty_dir)
  on.exit(unlink(empty_dir, recursive = TRUE))

  expect_error(
    load_PI_models(empty_dir),
    "No model files matching pattern"
  )
})
