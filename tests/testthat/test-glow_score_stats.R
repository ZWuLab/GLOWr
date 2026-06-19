# ==============================================================================
# Tests for compute_score_stats()
# ==============================================================================

test_that("compute_score_stats dispatches SPA for binary by default", {
  set.seed(42)
  n <- 100; p <- 5
  G <- matrix(rbinom(n * p, 2, 0.1), n, p)
  X <- cbind(1, rnorm(n))
  Y <- rbinom(n, 1, 0.3)
  null_model <- fit_null_model(X, Y, trait = "binary")

  result <- compute_score_stats(G, null_model, verbose = 0)
  expect_true(is.list(result))
  expect_equal(length(result$Zscores), p)
  expect_equal(dim(result$M_Z), c(p, p))
  expect_equal(result$s0, 1)
})

test_that("compute_score_stats dispatches standard for continuous", {
  set.seed(42)
  n <- 100; p <- 5
  G <- matrix(rbinom(n * p, 2, 0.2), n, p)
  X <- cbind(1, rnorm(n))
  Y <- rnorm(n)
  null_model <- fit_null_model(X, Y, trait = "continuous")

  result <- compute_score_stats(G, null_model, verbose = 0)
  expect_true(is.list(result))
  expect_equal(length(result$Zscores), p)
})

test_that("compute_score_stats respects use_spa override", {
  set.seed(42)
  n <- 100; p <- 5
  G <- matrix(rbinom(n * p, 2, 0.1), n, p)
  X <- cbind(1, rnorm(n))
  Y <- rbinom(n, 1, 0.3)
  null_model <- fit_null_model(X, Y, trait = "binary")

  # Force non-SPA for binary
  result <- compute_score_stats(G, null_model, use_spa = FALSE, verbose = 0)
  expect_true(is.list(result))
  expect_equal(length(result$Zscores), p)
})

test_that("compute_score_stats errors on SPA + continuous", {
  set.seed(42)
  n <- 100; p <- 5
  G <- matrix(rbinom(n * p, 2, 0.2), n, p)
  X <- cbind(1, rnorm(n))
  Y <- rnorm(n)
  null_model <- fit_null_model(X, Y, trait = "continuous")

  expect_error(
    compute_score_stats(G, null_model, use_spa = TRUE, verbose = 0),
    "binary"
  )
})

test_that("compute_score_stats validates dimension mismatch", {
  set.seed(42)
  G <- matrix(rbinom(50 * 3, 2, 0.1), 50, 3)
  X <- cbind(1, rnorm(100))
  Y <- rbinom(100, 1, 0.3)
  null_model <- fit_null_model(X, Y, trait = "binary")

  expect_error(
    compute_score_stats(G, null_model, verbose = 0),
    "row|dimension|sample"
  )
})

test_that("compute_score_stats validates null_model class", {
  set.seed(42)
  G <- matrix(rbinom(100 * 3, 2, 0.1), 100, 3)
  fake_model <- list(trait = "binary", n = 100)

  expect_error(
    compute_score_stats(G, fake_model, verbose = 0),
    "glow_null_model"
  )
})

test_that("compute_score_stats validates G is numeric matrix", {
  set.seed(42)
  n <- 100; p <- 3
  X <- cbind(1, rnorm(n))
  Y <- rbinom(n, 1, 0.3)
  null_model <- fit_null_model(X, Y, trait = "binary")

  # Not a matrix
  expect_error(
    compute_score_stats(rnorm(n), null_model, verbose = 0),
    "numeric matrix"
  )

  # Matrix with NAs
  G_na <- matrix(rbinom(n * p, 2, 0.1), n, p)
  G_na[1, 1] <- NA
  expect_error(
    compute_score_stats(G_na, null_model, verbose = 0),
    "NA"
  )
})

test_that("compute_score_stats emits message when verbose >= 1", {
  set.seed(42)
  n <- 100; p <- 3
  G <- matrix(rbinom(n * p, 2, 0.2), n, p)
  X <- cbind(1, rnorm(n))
  Y <- rnorm(n)
  null_model <- fit_null_model(X, Y, trait = "continuous")

  expect_message(
    compute_score_stats(G, null_model, verbose = 1),
    "Score stats"
  )
})

test_that("compute_score_stats suppresses message when verbose = 0", {
  set.seed(42)
  n <- 100; p <- 3
  G <- matrix(rbinom(n * p, 2, 0.2), n, p)
  X <- cbind(1, rnorm(n))
  Y <- rnorm(n)
  null_model <- fit_null_model(X, Y, trait = "continuous")

  expect_silent(
    compute_score_stats(G, null_model, verbose = 0)
  )
})

test_that("compute_score_stats SPA and non-SPA give same correlation structure", {
  # For binary trait, M_Z should be identical regardless of SPA
  set.seed(42)
  n <- 200; p <- 4
  G <- matrix(rbinom(n * p, 2, 0.15), n, p)
  X <- cbind(1, rnorm(n))
  Y <- rbinom(n, 1, 0.4)
  null_model <- fit_null_model(X, Y, trait = "binary")

  result_spa <- compute_score_stats(G, null_model, use_spa = TRUE, verbose = 0)
  result_std <- compute_score_stats(G, null_model, use_spa = FALSE, verbose = 0)

  # Correlation matrices should be identical (both based on variance structure)
  expect_equal(result_spa$M_Z, result_std$M_Z, tolerance = 1e-10)
  expect_equal(result_spa$M_s, result_std$M_s, tolerance = 1e-10)
})

test_that("compute_score_stats returns scores for both SPA and non-SPA", {
  set.seed(42)
  n <- 200; p <- 4
  G <- matrix(rbinom(n * p, 2, 0.15), n, p)
  X <- cbind(1, rnorm(n))
  Y <- rbinom(n, 1, 0.4)
  null_model <- fit_null_model(X, Y, trait = "binary")

  result_spa <- compute_score_stats(G, null_model, use_spa = TRUE, verbose = 0)
  result_std <- compute_score_stats(G, null_model, use_spa = FALSE, verbose = 0)

  # Both should have scores field
  expect_true("scores" %in% names(result_spa))
  expect_true("scores" %in% names(result_std))
  expect_equal(length(result_spa$scores), p)
  expect_equal(length(result_std$scores), p)

  # Raw scores should be identical regardless of SPA
  expect_equal(as.vector(result_spa$scores), as.vector(result_std$scores),
               tolerance = 1e-10)
})
