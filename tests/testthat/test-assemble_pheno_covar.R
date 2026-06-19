# Tests for assemble_pheno_covar() — the cohort-agnostic pheno/covar bundler (Phase-3 3b).

test_that("assembles a bundle in the documented shape, sample order preserved", {
  ids <- paste0("S", 1:6)
  pheno <- data.frame(iid = ids, status = c(2, 1, 2, 1, 2, -9),
                      sex = c("F", "M", "F", "M", "F", "F"),
                      age = c(50, 60, 55, 65, 58, 70), stringsAsFactors = FALSE)
  pcs <- data.frame(iid = ids, PC1 = 1:6 / 10, PC2 = 6:1 / 10)
  b <- assemble_pheno_covar(
    sample_ids = ids,
    outcome = list(values = pheno$status, map = c("1" = 0L, "2" = 1L)),
    trait = "binary", pheno = pheno, pheno_id_col = "iid",
    covariates = list(sex_numeric = list(col = "sex", map = c(F = 1, M = 0)), age = "age"),
    pcs = list(scores = pcs, id_col = "iid", cols = c("PC1", "PC2")),
    covar_set = "demo", verbose = 0)

  expect_named(b, c("sample_id", "Y", "X", "trait", "n_total",
                    "n_excluded", "covar_set", "covar_names"))
  expect_equal(colnames(b$X), c("sex_numeric", "age", "PC1", "PC2"))
  expect_equal(b$trait, "binary")
  expect_equal(b$covar_set, "demo")
  expect_equal(b$n_total, 6L)
  # S6 has status -9 (-> NA) and is dropped; the other 5 are complete
  expect_equal(b$sample_id, paste0("S", 1:5))
  expect_equal(b$n_excluded, 1L)
  expect_equal(unname(b$Y), c(1, 0, 1, 0, 1))            # 2->1, 1->0
  expect_equal(unname(b$X[, "sex_numeric"]), c(1, 0, 1, 0, 1))
  expect_equal(unname(b$X[, "PC1"]), (1:5) / 10)         # PC rows matched by id
})

test_that("recode maps send unmapped values to NA and such rows are dropped", {
  ids <- paste0("S", 1:4)
  pheno <- data.frame(iid = ids, grp = c("a", "b", "zzz", "a"),
                      y = c(1, 0, 1, 0), stringsAsFactors = FALSE)
  b <- assemble_pheno_covar(
    sample_ids = ids, outcome = list(values = pheno$y), trait = "binary",
    pheno = pheno, pheno_id_col = "iid",
    covariates = list(grp_num = list(col = "grp", map = c(a = 1, b = 2))),  # 'zzz' -> NA
    drop_incomplete = TRUE, verbose = 0)
  expect_equal(b$sample_id, c("S1", "S2", "S4"))   # S3 dropped (grp NA)
  expect_equal(unname(b$X[, "grp_num"]), c(1, 2, 1))
  expect_equal(b$n_excluded, 1L)
})

test_that("drop_incomplete = FALSE keeps NA-covariate rows (only NA-outcome dropped)", {
  ids <- paste0("S", 1:4)
  pheno <- data.frame(iid = ids, grp = c("a", "b", "zzz", "a"),  # S3 'zzz' -> NA covar
                      y = c(1, 0, 1, 0), stringsAsFactors = FALSE)  # all outcomes valid
  b_keep <- assemble_pheno_covar(
    sample_ids = ids, outcome = list(values = pheno$y), trait = "binary",
    pheno = pheno, pheno_id_col = "iid",
    covariates = list(grp_num = list(col = "grp", map = c(a = 1, b = 2))),
    drop_incomplete = FALSE, verbose = 0)
  expect_equal(b_keep$sample_id, paste0("S", 1:4))            # nothing dropped
  expect_true(is.na(b_keep$X[b_keep$sample_id == "S3", "grp_num"]))  # NA covar kept
  expect_equal(b_keep$n_excluded, 0L)

  b_drop <- assemble_pheno_covar(
    sample_ids = ids, outcome = list(values = pheno$y), trait = "binary",
    pheno = pheno, pheno_id_col = "iid",
    covariates = list(grp_num = list(col = "grp", map = c(a = 1, b = 2))),
    drop_incomplete = TRUE, verbose = 0)
  expect_equal(b_drop$sample_id, c("S1", "S2", "S4"))         # S3 dropped (NA covar)
})

test_that("PCs can be row-aligned (no id_col) and a function map works for outcome", {
  ids <- paste0("S", 1:3)
  pcmat <- cbind(PC1 = c(0.1, 0.2, 0.3))
  b <- assemble_pheno_covar(
    sample_ids = ids,
    outcome = list(values = c(10, 20, 30), map = function(v) as.integer(v > 15)),
    trait = "binary",
    pcs = list(scores = pcmat, id_col = NULL, cols = "PC1"), verbose = 0)
  expect_equal(unname(b$Y), c(0L, 1L, 1L))
  expect_equal(unname(b$X[, "PC1"]), c(0.1, 0.2, 0.3))
})

test_that("input validation errors are clear", {
  expect_error(assemble_pheno_covar(character(0), list(values = numeric(0))),
               "empty")
  expect_error(assemble_pheno_covar(c("a", "b"), list(values = 1)),
               "must equal length")
  expect_error(assemble_pheno_covar(c("a", "b"), list(map = c(x = 1))),
               "values")
})
