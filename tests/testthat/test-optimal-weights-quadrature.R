# ==============================================================================
# File: test-optimal-weights-quadrature.R
# Purpose: The quadrature and the Hermite expansion behind Optimal_Weights_M() (GLOWr
#          0.3.4): the kink-split composite Gauss-Legendre rule against closed forms and
#          against integrate() at a tight tolerance, the order-64 expansion of the Fisher
#          transform against a direct integration of every pairwise covariance (the
#          synthetic fixture), the order-8 expansion against its full-precision reference,
#          the structure of the covariances, the boundary cases, the reuse of coefficients,
#          the cache key, and hermite() by recurrence.

g_fisher <- function(x) g_GFisher_two(x, 2)        # -2 log(2 Phi(-|x|)), a kink at 0
g_square <- function(x) x^2

.fixture <- function() {
  path <- testthat::test_path("fixtures", "optimal-weights-synthetic.rds")
  skip_if_not(file.exists(path), "fixture optimal-weights-synthetic.rds not found")
  readRDS(path)
}
.relmax <- function(a, b) max(abs(a - b)) / max(abs(b))
# the standardized quantities of the normal approximation
.standardize <- function(Bstar, PI, M) {
  MU <- Bstar * PI; V <- Bstar^2 * PI * (1 - PI)
  Sigma <- M + t(M * V) %*% M
  list(mu = as.vector(M %*% MU), s = sqrt(diag(Sigma)), R = stats::cov2cor(Sigma), Sigma = Sigma)
}

# ---- G1: the rule ----

test_that("the composite rule integrates polynomial moments exactly, at every scale", {
  for (mu in c(-2, 0.3, 4)) for (s in c(0.5, 1, 7.5, 30)) {
    m <- GLOWr:::.hermite_moments(g_square, mu, s, 6)
    # b_0 = E[(mu + sU)^2] = mu^2 + s^2; b_1 = 2 mu s; b_2 = sqrt(2) s^2; b_3.. = 0
    b <- c(mu^2 + s^2, 2 * mu * s, sqrt(2) * s^2, 0, 0, 0, 0)
    expect_lt(max(abs(m$b - b)), 1e-12 * max(1, max(abs(b))))
    # E[(mu + sU)^4] = mu^4 + 6 mu^2 s^2 + 3 s^4, and the centered variance 4 mu^2 s^2 + 2 s^4
    Eg2 <- mu^4 + 6 * mu^2 * s^2 + 3 * s^4
    expect_lt(abs(m$Eg2 - Eg2), 1e-12 * max(1, Eg2))
    expect_lt(abs(m$Var - (4 * mu^2 * s^2 + 2 * s^4)), 1e-12 * max(1, Eg2))
  }
  expect_equal(GLOWr:::E_gX_p(g_square, 1.3, 1, 2.5), 1.3^2 + 2.5^2, tolerance = 1e-12)
  expect_equal(GLOWr:::E_gX_p(function(x) x, 1.3, 2, 2.5), 1.3^2 + 2.5^2, tolerance = 1e-12)
  expect_equal(GLOWr:::E_gX_p(function(x) x, 2, 3), 14, tolerance = 1e-12)   # E[(2 + U)^3]
  expect_equal(GLOWr:::get_r_tilde(g_square, c(0.5, 1, 0.2), c(1.1, 2, 1.5)),
               c(0.5, 1, 0.2)^2 + c(1.1, 2, 1.5)^2 - 1, tolerance = 1e-12)
})

test_that("the null Fisher transform has E[g] = 2 and Var[g] = 4 (a chi-square with 2 df)", {
  m <- GLOWr:::.hermite_moments(g_fisher, 0, 1, 64)
  expect_equal(m$b[1], 2, tolerance = 1e-10)
  expect_equal(m$Var, 4, tolerance = 1e-10)
  expect_equal(m$Eg2 - m$b[1]^2, 4, tolerance = 1e-10)
  # g is even: the odd coefficients vanish; the truncated variance rises toward 4
  expect_true(all(abs(m$b[seq(2, 65, by = 2)]) < 1e-10))
  partial <- cumsum(m$b[-1]^2)
  expect_lt(partial[8], partial[64])
  expect_lt(partial[64], 4)
  expect_gt(partial[64], 3.999)
})

test_that("the rule agrees with integrate() at a tight tolerance split at the kink", {
  mu <- 2.75; s <- 7.54                                      # a stretched variant
  m <- GLOWr:::.hermite_moments(g_fisher, mu, s, 8)
  x0 <- mu - mu / s                                          # where mu + s (x - mu) = 0
  # the same domain as the rule, mu +/- 12 on the standardized scale
  for (k in 1:8) {
    f <- function(x) g_fisher(mu + s * (x - mu)) * dnorm(x - mu) * GLOWr:::hermite(x, mu, k)
    a_k <- integrate(f, mu - 12, x0, rel.tol = 1e-10, subdivisions = 2000L)$value +
           integrate(f, x0, mu + 12, rel.tol = 1e-10, subdivisions = 2000L)$value
    b_k <- m$b[k + 1] * sqrt(factorial(k))
    expect_lt(abs(b_k - a_k), 1e-8 * abs(a_k) + 1e-8)
  }
})

test_that("a transform returning a non-finite value, or a non-finite moment, is refused", {
  expect_error(GLOWr:::.hermite_moments(function(x) ifelse(x > 50, Inf, x), 1, 10, 4), "non-finite")
  expect_error(suppressWarnings(GLOWr:::E_gX_p(function(x) log(x), 0, 1)), "non-finite")
  # g finite everywhere on the domain (exp(360) at x = 12) but its square overflows
  expect_error(GLOWr:::.hermite_moments(function(x) exp(30 * abs(x)), 0, 1, 2), "not finite")
  expect_error(GLOWr:::.hermite_moments(function(x) x[1], 0, 1, 2), "one value per input")
  expect_error(GLOWr:::.hermite_moments(g_square, 0, -1, 2), "positive")
})

# ---- G2: the weights against the direct integration ----

test_that("the order-64 weights of the fixture are within the truncation tolerance of the direct integration", {
  fx <- .fixture()
  w <- Optimal_Weights_M(g_fisher, fx$inputs$Bstar, fx$inputs$PI, fx$inputs$M, is.posi.wts = TRUE)
  expect_lt(.relmax(as.vector(w$wts_APE_N), fx$direct$APE_N), 1e-3)
  expect_lt(.relmax(as.vector(w$wts_BE_N), fx$direct$BE_N), 1e-3)
  expect_lt(.relmax(as.vector(w$wts_APE_sparse), fx$direct$APE_sparse), 1e-4)
  expect_lt(.relmax(as.vector(w$wts_BE_sparse), fx$direct$BE_sparse), 5e-3)
})

test_that("ORD = 1:8 reproduces the order-8 expansion with accurate coefficients", {
  fx <- .fixture()
  n <- length(fx$inputs$Bstar); M <- fx$inputs$M
  z <- .standardize(fx$inputs$Bstar, fx$inputs$PI, M)
  r <- GLOWr:::get_r_tilde(g_fisher, z$mu, z$s)
  expect_equal(r, fx$r_tilde, tolerance = 1e-10)
  S_APE <- GLOWr:::CovM_gXgY(g_fisher, z$mu, z$mu, z$R, ORD = 1:8, SD1 = z$s, SD2 = z$s)
  S_BE <- GLOWr:::CovM_gXgY(g_fisher, rep(0, n), rep(0, n), M, ORD = 1:8)
  w_APE <- as.vector(GLOWr:::get_wts(S_APE, r, TRUE)$w_normalized)
  w_BE <- as.vector(GLOWr:::get_wts(S_BE, r, TRUE)$w_normalized)
  expect_lt(.relmax(w_APE, fx$order8_accurate$APE_N), 1e-9)
  expect_lt(.relmax(w_BE, fx$order8_accurate$BE_N), 1e-9)
})

# ---- G3: the structure of the covariances ----

test_that("the _N covariance's diagonal is the truncated sum and the matrix is positive semidefinite", {
  fx <- .fixture()
  n <- 5L; M <- fx$inputs$M[1:n, 1:n]
  mu <- c(0.2, -1, 0.5, 2, 0); s <- c(1, 1.5, 3, 0.8, 1)
  S <- GLOWr:::CovM_gXgY(g_fisher, mu, mu, M, SD1 = s, SD2 = s)
  tab <- GLOWr:::.hermite_table(g_fisher, mu, s, 64)
  expect_equal(diag(S), colSums(tab$b[-1, , drop = FALSE]^2), tolerance = 1e-12)
  expect_true(isSymmetric(S, tol = 1e-12))
  expect_gt(min(eigen(S, symmetric = TRUE, only.values = TRUE)$values), -1e-8 * max(diag(S)))
})

test_that("the sparse covariance carries the mixture variance and matches the SKAT closed form", {
  fx <- .fixture()
  n <- 5L; M <- fx$inputs$M[1:n, 1:n]; Bstar <- fx$inputs$Bstar[1:n]; PI <- fx$inputs$PI[1:n]
  Smix <- GLOWr:::CovMT_mix(g_fisher, Bstar, PI, M)
  expect_equal(diag(Smix), sapply(1:n, function(i) GLOWr:::Var_T_mix(g_fisher, Bstar[i], PI[i])),
               tolerance = 1e-10)
  # g(x) = x^2: Z_i = U_i + Bstar_i C_i, Cov(Z_i^2, Z_j^2) = 2 M_ij^2 + 4 pi_i pi_j Bstar_i Bstar_j M_ij
  # off the diagonal, and the total variance Var(Z_i^2) on it
  S2 <- GLOWr:::CovMT_mix(g_square, Bstar, PI, M)
  expected <- 2 * M^2 + 4 * outer(PI * Bstar, PI * Bstar) * M
  # Var(Z^2) with Z = U + B C: E[Z^4] - E[Z^2]^2 from the moments of the mixture
  EZ2 <- 1 + PI * Bstar^2
  EZ4 <- PI * (Bstar^4 + 6 * Bstar^2 + 3) + (1 - PI) * 3
  diag(expected) <- EZ4 - EZ2^2
  expect_equal(S2, expected, tolerance = 1e-10)
  # and get_Sigma()/get_r() agree with the one-table moments of Optimal_Weights_M()
  mom <- GLOWr:::.optimal_weight_moments(g_fisher, Bstar, PI, M)
  expect_equal(mom$Sigma_APE_sparse, Smix, tolerance = 1e-12)
  expect_equal(mom$Sigma_BE_sparse, GLOWr:::get_Sigma(g_fisher, Bstar, PI, M, "H0"), tolerance = 1e-12)
  expect_equal(mom$r_sparse, GLOWr:::get_r(g_fisher, Bstar, PI), tolerance = 1e-12)
  z <- .standardize(Bstar, PI, M)
  expect_equal(mom$r_N, GLOWr:::get_r_tilde(g_fisher, z$mu, z$s), tolerance = 1e-12)
  expect_equal(mom$Sigma_APE_N, GLOWr:::CovM_gXgY(g_fisher, z$mu, z$mu, z$R, SD1 = z$s, SD2 = z$s), tolerance = 1e-12)
  expect_equal(mom$Sigma_BE_N, GLOWr:::CovM_gXgY(g_fisher, rep(0, n), rep(0, n), M), tolerance = 1e-12)
})

test_that("the SKAT transform's weights are the same at any order (its expansion is exact at order 2)", {
  fx <- .fixture()
  n <- 6L; M <- fx$inputs$M[1:n, 1:n]; Bstar <- fx$inputs$Bstar[1:n]; PI <- fx$inputs$PI[1:n]
  z <- .standardize(Bstar, PI, M)
  S2 <- GLOWr:::CovM_gXgY(g_square, z$mu, z$mu, z$R, ORD = 1:2, SD1 = z$s, SD2 = z$s)
  S64 <- GLOWr:::CovM_gXgY(g_square, z$mu, z$mu, z$R, SD1 = z$s, SD2 = z$s)
  expect_equal(S2, S64, tolerance = 1e-10)
  expect_equal(S64, 2 * z$Sigma^2 + 4 * outer(z$mu, z$mu) * z$Sigma, tolerance = 1e-10)
  w2 <- GLOWr:::.optimal_weight_moments(g_square, Bstar, PI, M, ORD = 1:2)
  w64 <- GLOWr:::.optimal_weight_moments(g_square, Bstar, PI, M)
  expect_equal(w2, w64, tolerance = 1e-10)
})

# ---- boundaries ----

test_that("K = 0 and K = 1, and ORD subsets, keep their dimensions and values", {
  m0 <- GLOWr:::.hermite_moments(g_fisher, 0.5, 2, 0)
  m1 <- GLOWr:::.hermite_moments(g_fisher, 0.5, 2, 1)
  expect_length(m0$b, 1L); expect_length(m1$b, 2L)
  expect_equal(m0$b[1], m1$b[1])
  tab <- GLOWr:::.hermite_table(g_fisher, c(0.5, 1), c(2, 1), 0)
  expect_equal(dim(tab$b), c(1L, 2L))
  M <- matrix(c(1, 0.4, 0.4, 1), 2, 2); mu <- c(0.5, -0.2)
  S2 <- GLOWr:::CovM_gXgY(g_fisher, mu, mu, M, ORD = 2L)
  S24 <- GLOWr:::CovM_gXgY(g_fisher, mu, mu, M, ORD = c(2L, 4L))
  S4 <- GLOWr:::CovM_gXgY(g_fisher, mu, mu, M, ORD = 4L)
  expect_equal(dim(S2), c(2L, 2L)); expect_equal(dim(S24), c(2L, 2L))
  expect_equal(S24, S2 + S4, tolerance = 1e-12)
  tab <- GLOWr:::.hermite_table(g_fisher, mu, c(1, 1), 4)
  expect_equal(S2, tcrossprod(tab$b[3, ]) * M^2, tolerance = 1e-12)
  expect_error(GLOWr:::CovM_gXgY(g_fisher, mu, mu, M, ORD = 0L), "positive integers")
  expect_error(GLOWr:::CovM_gXgY(g_fisher, mu, mu, M, ORD = 1.5), "positive integers")
})

test_that("one variant, repeated pairs, zero and negative Bstar, and endpoint PI all work", {
  fx <- .fixture()
  w1 <- Optimal_Weights_M(g_fisher, 3, 0.5, matrix(1, 1, 1))
  expect_true(all(is.finite(unlist(w1))))
  n <- 4L; M <- fx$inputs$M[1:n, 1:n]
  w <- Optimal_Weights_M(g_fisher, c(2, 2, 2, 2), c(0.1, 0.1, 0.1, 0.1), M)   # repeated pairs
  expect_true(all(is.finite(unlist(w))))
  w <- Optimal_Weights_M(g_fisher, c(0, -3, 5, 0), c(0.2, 0.3, 0.1, 0.4), M)  # zero and negative
  expect_true(all(is.finite(unlist(w))))
  for (PI in list(c(0, 0, 0, 0), c(1, 1, 1, 1), c(1, 0, 0, 1))) {
    z <- .standardize(fx$inputs$Bstar[1:n], PI, M)
    expect_equal(z$s, rep(1, n))                                   # MVM = 0
    w <- suppressWarnings(Optimal_Weights_M(g_fisher, fx$inputs$Bstar[1:n], PI, M))
    expect_true(all(is.finite(unlist(w))))
  }
})

# ---- reuse and the cache ----

test_that("each distinct (mu, s) pair is integrated once, in the table and in Optimal_Weights_M()", {
  calls <- 0L
  g_count <- function(x) { calls <<- calls + 1L; x^2 }
  tab <- GLOWr:::.hermite_table(g_count, c(0, 1, 0, 1, 0), c(1, 2, 1, 2, 1), 4)
  expect_equal(calls, 2L)
  expect_equal(dim(tab$b), c(5L, 5L))
  expect_equal(tab$b[, 1], tab$b[, 3])
  expect_equal(tab$Eg2[2], tab$Eg2[4])
  fx <- .fixture()
  n <- 5L; M <- fx$inputs$M[1:n, 1:n]; Bstar <- fx$inputs$Bstar[1:n]; PI <- fx$inputs$PI[1:n]
  calls <- 0L
  GLOWr:::clear_glow_cache()
  invisible(Optimal_Weights_M(g_count, Bstar, PI, M))
  expect_equal(calls, 2L * n + 1L)   # the null pair, n normal-approximation pairs, n sparse pairs
})

test_that("the cached null covariance is keyed by the null coefficients, M and ORD", {
  fx <- .fixture()
  n <- 4L; M <- fx$inputs$M[1:n, 1:n]; Bstar <- fx$inputs$Bstar[1:n]; PI <- fx$inputs$PI[1:n]
  # closures with the same body and different captured values
  make_g <- function(df) function(x) g_GFisher_two(x, df)
  g2 <- make_g(2); g3 <- make_g(3)
  GLOWr:::clear_glow_cache()
  warm2 <- GLOWr:::CovMT_mix(g2, Bstar, PI, M)
  warm3 <- GLOWr:::CovMT_mix(g3, Bstar, PI, M)   # must not reuse g2's null covariance
  GLOWr:::clear_glow_cache()
  cold3 <- GLOWr:::CovMT_mix(g3, Bstar, PI, M)
  expect_equal(warm3, cold3, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(warm2, warm3, tolerance = 1e-3)))
  # two transforms that agree at a probe set and differ elsewhere
  g_a <- function(x) x^2
  g_b <- function(x) x^2 + 0.1 * (x^2 * (x^2 - 0.25) * (x^2 - 2.25) * (x^2 - 9) * exp(-x^2))^2
  probe <- c(-3, -1.5, -0.5, 0, 0.5, 1.5, 3)
  expect_equal(g_a(probe), g_b(probe))
  GLOWr:::clear_glow_cache()
  invisible(GLOWr:::CovMT_mix(g_a, Bstar, PI, M))
  warm_b <- GLOWr:::CovMT_mix(g_b, Bstar, PI, M)
  GLOWr:::clear_glow_cache()
  cold_b <- GLOWr:::CovMT_mix(g_b, Bstar, PI, M)
  expect_equal(warm_b, cold_b, tolerance = 1e-12)
  # different orders do not share an entry, and the cache can be switched off
  GLOWr:::clear_glow_cache()
  s8 <- GLOWr:::CovMT_mix(g_fisher, Bstar, PI, M, ORD = 1:8)
  s64 <- GLOWr:::CovMT_mix(g_fisher, Bstar, PI, M)
  expect_false(isTRUE(all.equal(s8, s64, tolerance = 1e-8)))
  old <- options(GLOWr.use_cache = FALSE); on.exit(options(old), add = TRUE)
  expect_equal(GLOWr:::CovMT_mix(g_fisher, Bstar, PI, M), s64, tolerance = 1e-12)
})

test_that("hermite() by recurrence equals the closed forms and accepts any positive order", {
  x <- seq(-3, 3, by = 0.5); mu <- 0.5; y <- x - mu
  closed <- list(y, y^2 - 1, y^3 - 3 * y, y^4 - 6 * y^2 + 3, y^5 - 10 * y^3 + 15 * y,
                 y^6 - 15 * y^4 + 45 * y^2 - 15, y^7 - 21 * y^5 + 105 * y^3 - 105 * y,
                 y^8 - 28 * y^6 + 210 * y^4 - 420 * y^2 + 105)
  for (k in 1:8) expect_equal(GLOWr:::hermite(x, mu, k), closed[[k]], tolerance = 1e-12)
  # He_12(1.5) by the three-term recurrence, computed here independently
  h0 <- 1; h1 <- 1.5; for (k in 1:11) { h2 <- 1.5 * h1 - k * h0; h0 <- h1; h1 <- h2 }
  expect_equal(GLOWr:::hermite(1.5, 0, 12), h1, tolerance = 1e-12)
  expect_error(GLOWr:::hermite(0, 0, 0), "positive integer")
  expect_error(GLOWr:::hermite(0, 0, 2.5), "positive integer")
})
