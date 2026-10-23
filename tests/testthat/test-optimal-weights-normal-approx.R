# ==============================================================================
# File: test-optimal-weights-normal-approx.R
# Purpose: The normal approximation of Optimal_Weights_M() (the _N schemes) on the
#          standardized scale: CovM_gXgY() with scales, get_r_tilde() with standard
#          deviations, and the resulting weights, each against a closed form that holds
#          for g(x) = x and g(x) = x^2 (the Hermite expansion of a polynomial of degree
#          two is exact at order two, so only the quadrature error remains).

g_identity <- function(x) x
g_square <- function(x) x^2

# A covariance with a non-unit diagonal, as M + MVM is, and its standardization
.make_sigma <- function(seed = 1, n = 4) {
  set.seed(seed)
  L <- matrix(rnorm(n * 2), n, 2)
  Sigma <- L %*% t(L) + diag(runif(n, 0.5, 2))
  list(Sigma = Sigma, R = stats::cov2cor(Sigma), s = sqrt(diag(Sigma)))
}

test_that("CovM_gXgY with scales returns the covariance itself for the identity transformation", {
  z <- .make_sigma(1)
  mu <- c(0.3, -0.5, 1.2, 0)
  Cov <- GLOWr:::CovM_gXgY(g_identity, mu, mu, z$R, SD1 = z$s, SD2 = z$s)
  expect_equal(Cov, z$Sigma, tolerance = 1e-6)
})

test_that("CovM_gXgY with scales matches the closed form for g(x) = x^2", {
  z <- .make_sigma(2)
  m1 <- c(0.3, -0.5, 1.2, 0)
  m2 <- c(1.0, 0.2, -0.4, 0.7)
  # X_i = m1_i + s_i Z_i, Y_j = m2_j + s_j Z_j, Cov(Z_i, Z_j) = R_ij:
  # Cov(X_i^2, Y_j^2) = 2 Sigma_ij^2 + 4 m1_i m2_j Sigma_ij with Sigma_ij = s_i s_j R_ij
  expected <- 2 * z$Sigma^2 + 4 * outer(m1, m2) * z$Sigma
  Cov <- GLOWr:::CovM_gXgY(g_square, m1, m2, z$R, SD1 = z$s, SD2 = z$s)
  expect_equal(Cov, expected, tolerance = 1e-6)
})

test_that("CovM_gXgY defaults reproduce the unit-scale call", {
  M <- matrix(c(1, 0.4, 0.4, 1), 2, 2)
  mu <- c(0.5, -0.2)
  a <- GLOWr:::CovM_gXgY(g_square, mu, mu, M)
  b <- GLOWr:::CovM_gXgY(g_square, mu, mu, M, SD1 = c(1, 1), SD2 = c(1, 1))
  expect_identical(a, b)
  expect_equal(a, 2 * M^2 + 4 * outer(mu, mu) * M, tolerance = 1e-6)
})

test_that("the identity fast path uses both mean vectors and both scales", {
  # M = I: only the diagonal is computed; with X_i = m1_i + s1_i Z, Y_i = m2_i + s2_i Z,
  # Cov(X_i^2, Y_i^2) = 2 s1_i^2 s2_i^2 + 4 m1_i m2_i s1_i s2_i
  m1 <- c(0.5, 1.5); m2 <- c(0, -1); s1 <- c(1, 2); s2 <- c(1.5, 0.5)
  Cov <- GLOWr:::CovM_gXgY(g_square, m1, m2, diag(2), SD1 = s1, SD2 = s2)
  expected <- diag(2 * s1^2 * s2^2 + 4 * m1 * m2 * s1 * s2)
  expect_equal(Cov, expected, tolerance = 1e-6)
  # the unit-scale, equal-mean case is the variance 2 + 4 m^2 as before
  Cov0 <- GLOWr:::CovM_gXgY(g_square, m1, m1, diag(2))
  expect_equal(Cov0, diag(2 + 4 * m1^2), tolerance = 1e-6)
})

test_that("CovM_gXgY refuses a covariance in place of a correlation, and non-positive scales", {
  z <- .make_sigma(3)
  mu <- rep(0, 4)
  expect_error(GLOWr:::CovM_gXgY(g_square, mu, mu, z$Sigma), "correlation matrix")
  expect_error(GLOWr:::CovM_gXgY(g_square, mu, mu, z$R, SD1 = c(1, 0, 1, 1), SD2 = rep(1, 4)), "positive")
  expect_error(GLOWr:::CovM_gXgY(g_square, mu, mu[1:3], z$R), "same length")
})

test_that("get_r_tilde takes standard deviations: E[(mu + s Z)^2] - 1 = mu^2 + s^2 - 1", {
  mu <- c(0.5, 1.0, 0.2)
  s <- c(1.1, 2.0, 1.5)
  expect_equal(GLOWr:::get_r_tilde(g_square, mu, s), mu^2 + s^2 - 1, tolerance = 1e-6)
})

test_that("the _N weights for g(x) = x^2 match the closed form of the standardized normal approximation", {
  set.seed(11)
  n <- 5
  L <- matrix(rnorm(n * 2), n, 2)
  M <- stats::cov2cor(L %*% t(L) + diag(runif(n, 0.3, 1)))
  Bstar <- runif(n, 0.5, 2)
  PI <- runif(n, 0.2, 0.8)

  MU <- Bstar * PI
  V <- Bstar^2 * PI * (1 - PI)
  MMU <- as.vector(M %*% MU)
  Sigma <- M + t(M * V) %*% M                 # Var(Z~) under the random-effects model

  # Closed forms for T_i = Z_i^2: r_i = E[Z~_i^2] - E[Z_0^2] = MMU_i^2 + Sigma_ii - 1,
  # Sigma_BE = 2 M o M (under H0), Sigma_APE = 2 Sigma o Sigma + 4 (MMU MMU') o Sigma
  r <- MMU^2 + diag(Sigma) - 1
  Sigma_BE <- 2 * M^2
  Sigma_APE <- 2 * Sigma^2 + 4 * outer(MMU, MMU) * Sigma
  expected_BE <- GLOWr:::get_wts(Sigma_BE, r, is.posi.wts = TRUE)$w_normalized
  expected_APE <- GLOWr:::get_wts(Sigma_APE, r, is.posi.wts = TRUE)$w_normalized

  wts <- Optimal_Weights_M(g_square, Bstar, PI, M, is.posi.wts = TRUE)
  expect_equal(as.vector(wts$wts_BE_N), as.vector(expected_BE), tolerance = 1e-5)
  expect_equal(as.vector(wts$wts_APE_N), as.vector(expected_APE), tolerance = 1e-5)

  # The pre-0.3.3 form (the covariance in the expansion, the variance as the SD) is not this
  r_old <- GLOWr:::get_r_tilde(g_square, MMU, diag(Sigma))
  expect_false(isTRUE(all.equal(r_old, r, tolerance = 1e-3)))
})
