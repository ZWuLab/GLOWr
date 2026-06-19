# Test file for optimal weight helper functions
#
# This file tests all helper functions in helpers_optimalWeights.R

library(testthat)
library(GLOWr)

# Temporary g_GFisher for testing (will be properly ported in Phase 3)
g_GFisher_two <- function(x, df = 2) {
  qchisq(log(2) + pnorm(abs(x), lower.tail = FALSE, log.p = TRUE),
         df = df, lower.tail = FALSE, log.p = TRUE)
}

# Identity transformation
g_identity <- function(x) x

#################### Tests for select_best_model ####################

test_that("select_best_model works with linear relationship", {
  set.seed(123)
  X <- runif(100, 0.01, 0.5)
  Y <- 2 * X + 0.5 + rnorm(100, 0, 0.1)

  model <- select_best_model(X, Y)

  # Model should be a linear model
  expect_s3_class(model, "lm")

  # R-squared should be reasonably high for this linear relationship
  expect_gt(summary(model)$r.squared, 0.7)
})

test_that("select_best_model handles f(X) relationship", {
  set.seed(123)
  X <- runif(100, 0.01, 0.5)
  Y <- X * (1 - X) + rnorm(100, 0, 0.01)

  model <- select_best_model(X, Y)

  expect_s3_class(model, "lm")
  expect_gt(summary(model)$r.squared, 0.8)
})

test_that("select_best_model validates input", {
  X <- c(0.1, 0.2, 0.3)
  Y <- c(1, 2)

  expect_error(select_best_model(X, Y), "same length")

  X_bad <- c(-0.1, 0.2, 0.3)
  Y_ok <- c(1, 2, 3)
  expect_error(select_best_model(X_bad, Y_ok), "non-negative|interval")  # New validation catches this earlier

  X_ok <- c(0.1, 0.2, 0.3)
  Y_bad <- c(-1, 2, 3)
  expect_error(select_best_model(X_ok, Y_bad), "non-negative|positive")  # New validation catches this earlier
})


#################### Tests for getTrainData ####################

test_that("getTrainData creates proper training data", {
  set.seed(123)
  case_anno <- matrix(runif(50 * 5), ncol = 5)
  control_anno <- matrix(runif(1000 * 5), ncol = 5)

  train_data <- getTrainData(case_anno, control_anno, control_need_N = 100)

  # Check structure
  expect_type(train_data, "list")
  expect_named(train_data, c("x", "y"))

  # Check dimensions
  expect_equal(nrow(train_data$x), 150)  # 50 cases + 100 controls
  expect_equal(ncol(train_data$x), 5)    # 5 annotation features
  expect_equal(nrow(train_data$y), 150)
  expect_equal(ncol(train_data$y), 1)

  # Check labels: first 50 should be 1 (cases), rest should be 0 (controls)
  expect_equal(sum(train_data$y == 1), 50)
  expect_equal(sum(train_data$y == 0), 100)
})

test_that("getTrainData validates control_need_N", {
  case_anno <- matrix(runif(10 * 3), ncol = 3)
  control_anno <- matrix(runif(50 * 3), ncol = 3)

  expect_error(
    getTrainData(case_anno, control_anno, control_need_N = 100),
    "exceeds"
  )
})


#################### Tests for Hermite Polynomials ####################

test_that("hermite computes correct polynomial values", {
  x <- 0
  mu <- 0

  # Order 1: He_1(x) = x
  expect_equal(hermite(x, mu, 1), 0)

  # Order 2: He_2(x) = x^2 - 1
  expect_equal(hermite(x, mu, 2), -1)

  # Order 3: He_3(x) = x^3 - 3x
  expect_equal(hermite(x, mu, 3), 0)

  # Order 4: He_4(x) = x^4 - 6x^2 + 3
  expect_equal(hermite(x, mu, 4), 3)
})

test_that("hermite handles shift parameter correctly", {
  x <- 2
  mu <- 1

  # Order 1: He_1(x-mu) = (x-mu) = (2-1) = 1
  expect_equal(hermite(x, mu, 1), 1)

  # Order 2: He_2(x-mu) = (x-mu)^2 - 1 = 1 - 1 = 0
  expect_equal(hermite(x, mu, 2), 0)

  # Order 3: He_3(x-mu) = (x-mu)^3 - 3(x-mu) = 1 - 3 = -2
  expect_equal(hermite(x, mu, 3), -2)
})

test_that("hermite is vectorized", {
  x <- c(0, 1, 2)
  mu <- 0

  result <- hermite(x, mu, 1)

  expect_length(result, 3)
  expect_equal(result, c(0, 1, 2))
})

test_that("hermite validates order", {
  expect_error(hermite(0, 0, 0), "between 1 and 8")
  expect_error(hermite(0, 0, 9), "between 1 and 8")
})


#################### Tests for E_gX_p ####################

test_that("E_gX_p computes correct expectation for identity", {
  g <- g_identity

  # E[X] for X ~ N(0, 1) should be 0
  result <- E_gX_p(g, mu = 0, p = 1, sigma = 1)
  expect_equal(result, 0, tolerance = 1e-6)

  # E[X^2] for X ~ N(0, 1) should be 1
  result <- E_gX_p(g, mu = 0, p = 2, sigma = 1)
  expect_equal(result, 1, tolerance = 1e-6)

  # E[X] for X ~ N(2, 1) should be 2
  result <- E_gX_p(g, mu = 2, p = 1, sigma = 1)
  expect_equal(result, 2, tolerance = 1e-6)
})

test_that("E_gX_p handles non-identity transformations", {
  g_square <- function(x) x^2

  # E[X^2] for X ~ N(0, 1) should be 1
  result <- E_gX_p(g_square, mu = 0, p = 1, sigma = 1)
  expect_equal(result, 1, tolerance = 1e-6)

  # E[(X^2)^2] = E[X^4] for X ~ N(0, 1) should be 3
  result <- E_gX_p(g_square, mu = 0, p = 2, sigma = 1)
  expect_equal(result, 3, tolerance = 1e-5)
})


#################### Tests for Mixture Distribution Functions ####################

test_that("E_T_mix computes correct expectation", {
  g <- g_identity

  # When pi = 0, should equal E[g(Z_0)]
  result <- E_T_mix(g, mu = 2, pi = 0)
  expect_equal(result, 0, tolerance = 1e-6)

  # When pi = 1, should equal E[g(Z_0 + mu)]
  result <- E_T_mix(g, mu = 2, pi = 1)
  expect_equal(result, 2, tolerance = 1e-6)

  # When pi = 0.5, should be midpoint
  result <- E_T_mix(g, mu = 2, pi = 0.5)
  expect_equal(result, 1, tolerance = 1e-6)
})

test_that("Var_T_mix computes correct variance", {
  g <- g_identity

  # When pi = 0, should equal Var[Z_0] = 1
  result <- Var_T_mix(g, mu = 2, pi = 0)
  expect_equal(result, 1, tolerance = 1e-6)

  # When pi = 1, should equal Var[Z_0 + mu] = 1 (shift doesn't change variance)
  result <- Var_T_mix(g, mu = 2, pi = 1)
  expect_equal(result, 1, tolerance = 1e-6)

  # When 0 < pi < 1, variance includes mixture component
  result <- Var_T_mix(g, mu = 2, pi = 0.5)
  expect_gt(result, 1)  # Should be greater than 1 due to mixture
})


#################### Tests for Weight Calculation Functions ####################

test_that("get_wts computes optimal weights correctly", {
  # Simple independent case
  Sigma <- diag(3)
  r <- c(1, 2, 0.5)

  wts <- get_wts(Sigma, r, is.posi.wts = FALSE)

  # For independent statistics with identity covariance, w = r
  expect_equal(as.vector(wts$w), r, tolerance = 1e-10)

  # Normalized weights should be w / mean(|w|)
  expected_norm <- r / mean(abs(r))
  expect_equal(as.vector(wts$w_normalized), expected_norm, tolerance = 1e-10)
})

test_that("get_wts handles positive weights constraint", {
  Sigma <- diag(3)
  r <- c(1, -0.5, 2)

  wts <- get_wts(Sigma, r, is.posi.wts = TRUE)

  # Negative weight should be set to 0
  expect_equal(as.vector(wts$w)[2], 0)
  expect_equal(as.vector(wts$w)[1], 1, tolerance = 1e-10)
  expect_equal(as.vector(wts$w)[3], 2, tolerance = 1e-10)
})

test_that("get_r_tilde computes mean differences", {
  g <- g_identity
  MU <- c(0.5, 1.0, 0.2)
  SD <- c(1.1, 1.2, 1.05)

  r_tilde <- get_r_tilde(g, MU, SD)

  # For identity function: E[Z_i] - E[Z_0] = MU_i - 0 = MU_i
  expect_equal(r_tilde, MU, tolerance = 1e-6)
})

test_that("get_r computes mean differences for mixture", {
  g <- g_identity
  MU <- c(0.5, 1.0, 0.2)
  PI <- c(0.01, 0.05, 0.001)

  r <- get_r(g, MU, PI)

  # For identity: E_mix - E_null = pi * mu - 0 = pi * mu
  expected <- MU * PI
  expect_equal(r, expected, tolerance = 1e-6)
})


#################### Tests for Covariance Matrix Functions ####################

test_that("CovM_gXgY computes identity covariance correctly", {
  skip_if_not(interactive(), "Hermite expansion tests are slow")

  g <- g_identity
  M <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  MU1 <- c(0, 0)
  MU2 <- c(0, 0)

  Cov <- CovM_gXgY(g, MU1, MU2, M)

  # For identity transformation with zero means, should equal M
  expect_equal(Cov, M, tolerance = 1e-4)
})


#################### Integration Test for Optimal_Weights_M ####################

test_that("Optimal_Weights_M works for Burden test", {
  g <- g_identity
  Bstar <- c(0.5, 1.0, 0.3)
  PI <- c(0.01, 0.05, 0.02)
  M <- diag(3)

  weights <- Optimal_Weights_M(g, Bstar, PI, M)

  # Should return only BE and APE weights for Burden
  expect_named(weights, c("wts_BE", "wts_APE"))

  # Weights should be vectors of length 3
  expect_length(weights$wts_BE, 3)
  expect_length(weights$wts_APE, 3)

  # All weights should be non-negative (default is.posi.wts = TRUE)
  expect_true(all(weights$wts_BE >= 0))
  expect_true(all(weights$wts_APE >= 0))

  # Weights should be normalized (mean of abs = 1)
  expect_equal(mean(abs(weights$wts_BE)), 1, tolerance = 1e-10)
  expect_equal(mean(abs(weights$wts_APE)), 1, tolerance = 1e-10)
})

test_that("Optimal_Weights_M validates inputs", {
  g <- g_identity
  Bstar <- c(0.5, 1.0, 0.3)
  PI <- c(0.01, 0.05)  # Wrong length
  M <- diag(3)

  expect_error(Optimal_Weights_M(g, Bstar, PI, M), "same length")

  PI <- c(0.01, 0.05, 1.5)  # Invalid PI value
  expect_error(Optimal_Weights_M(g, Bstar, PI, M), "\\[0, 1\\]")

  PI <- c(0.01, 0.05, 0.02)
  M_wrong <- diag(2)  # Wrong dimension
  expect_error(Optimal_Weights_M(g, Bstar, PI, M_wrong), "compatible dimensions")
})
