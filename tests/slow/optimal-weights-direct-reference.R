# ==============================================================================
# File: optimal-weights-direct-reference.R
# Purpose: Regenerate the reference fixture of the optimal-weight tests: the optimal
#          weights of one 20-variant window under the normal approximation (BE_N, APE_N)
#          and the sparse mixture (BE_sparse, APE_sparse) for the Fisher transform
#          g(x) = -2 log(2 Phi(-|x|)), computed by a DIRECT nested one-dimensional
#          integration of every pairwise covariance (no Hermite expansion), plus the
#          order-8 expansion with accurate coefficients (the reference of Victoria Pan's
#          report of 2026-10-29) at full precision. The package fixture is synthetic
#          (tests/testthat/fixtures/optimal-weights-synthetic-inputs.R); the script takes
#          any inputs file that defines Bstar, PI and M.
#
#          Not run by R CMD check (tests/slow/ is not a test directory). About six
#          minutes. Run from the package root:
#            Rscript tests/slow/optimal-weights-direct-reference.R \
#              [tests/testthat/fixtures/optimal-weights-synthetic-inputs.R]
#          The fixture is written beside the inputs file, its stem without "-inputs".
args <- commandArgs(trailingOnly = TRUE)
inputs_path <- if (length(args) >= 1L) args[[1L]] else
  file.path("tests", "testthat", "fixtures", "optimal-weights-synthetic-inputs.R")
if (!file.exists(inputs_path)) stop("inputs file not found: ", inputs_path, " (run from the package root)")
out_path <- sub("-inputs\\.R$", ".rds", inputs_path)
if (identical(out_path, inputs_path)) out_path <- sub("\\.R$", ".rds", inputs_path)
t0 <- Sys.time()
source(inputs_path)   # defines Bstar, PI, M
n <- length(Bstar)
g <- function(x) -2 * (log(2) + pnorm(abs(x), lower.tail = FALSE, log.p = TRUE))   # g_GFisher_two(x, 2)

# ---- composite Gauss-Legendre with a cut at a kink ------------------------------------
gl_rule <- function(m) { i <- seq_len(m - 1); b <- i / sqrt(4 * i^2 - 1); J <- matrix(0, m, m)
  J[cbind(i, i + 1)] <- b; J[cbind(i + 1, i)] <- b; e <- eigen(J, symmetric = TRUE); o <- order(e$values)
  list(x = e$values[o], w = 2 * e$vectors[1, o]^2) }
RULE <- gl_rule(10)
panel_nodes <- function(a, b, panels) { edges <- seq(a, b, length.out = panels + 1); h <- diff(edges) / 2
  c0 <- (edges[-1] + edges[-length(edges)]) / 2
  list(x = c(outer(RULE$x, h) + rep(c0, each = length(RULE$x))), w = c(outer(RULE$w, h))) }
nodes_split <- function(a, b, kink, panels) { if (is.finite(kink) && kink > a && kink < b) {
  L <- panel_nodes(a, kink, panels); U <- panel_nodes(kink, b, panels); list(x = c(L$x, U$x), w = c(L$w, U$w)) } else panel_nodes(a, b, 2 * panels) }
LIM <- 10
E_g <- function(mu, s, p = 1, panels = 200) { nd <- nodes_split(-LIM, LIM, -mu / s, panels); sum(g(mu + s * nd$x)^p * dnorm(nd$x) * nd$w) }
# Cov{g(mi + si U), g(mj + sj V)} with (U, V) standard bivariate normal, correlation r:
# V = r U + c W, the inner integral over W for each outer node, both cut at the kink of g
cov_pair <- function(mi, si, mj, sj, r, p_out = 100, p_in = 50) {
  if (abs(r) >= 1) stop("|r| must be below 1 for a pair")
  cc <- sqrt(1 - r^2); out <- nodes_split(-LIM, LIM, -mi / si, p_out)
  h <- vapply(out$x, function(u) { nd <- nodes_split(-LIM, LIM, (-mj / sj - r * u) / cc, p_in)
    sum(g(mj + sj * (r * u + cc * nd$x)) * dnorm(nd$x) * nd$w) }, 0)
  sum(g(mi + si * out$x) * dnorm(out$x) * out$w * h) - E_g(mi, si) * E_g(mj, sj)
}
var_one <- function(m, s) E_g(m, s, 2) - E_g(m, s)^2

# ---- the normal approximation ----------------------------------------------------------
MU <- Bstar * PI; V <- Bstar^2 * PI * (1 - PI); mu <- c(M %*% MU)
Sigma_tilde <- M + t(M * V) %*% M; s <- sqrt(diag(Sigma_tilde)); R <- stats::cov2cor(Sigma_tilde)
r_tilde <- vapply(1:n, function(i) E_g(mu[i], s[i]), 0) - E_g(0, 1)
Cov_APE <- matrix(0, n, n); Cov_BE <- matrix(0, n, n)
for (i in 1:n) for (j in i:n) {
  Cov_APE[i, j] <- Cov_APE[j, i] <- if (i == j) var_one(mu[i], s[i]) else cov_pair(mu[i], s[i], mu[j], s[j], R[i, j])
  Cov_BE[i, j] <- Cov_BE[j, i] <- if (i == j) var_one(0, 1) else cov_pair(0, 1, 0, 1, M[i, j])
}
cat(sprintf("normal approximation done (%.0f s)\n", as.numeric(Sys.time() - t0, units = "secs")))

# ---- the sparse mixture -------------------------------------------------------------------
E1 <- vapply(1:n, function(i) E_g(Bstar[i], 1), 0); E0 <- E_g(0, 1)
ET <- PI * E1 + (1 - PI) * E0
VarT <- PI * vapply(1:n, function(i) E_g(Bstar[i], 1, 2), 0) + (1 - PI) * E_g(0, 1, 2) - ET^2
r_sp <- ET - E0
Cov_sp <- matrix(0, n, n)
for (i in 1:n) for (j in 1:n) if (i < j) {
  rr <- M[i, j]
  e11 <- cov_pair(Bstar[i], 1, Bstar[j], 1, rr) + E1[i] * E1[j]
  e10 <- cov_pair(Bstar[i], 1, 0, 1, rr) + E1[i] * E0
  e01 <- cov_pair(0, 1, Bstar[j], 1, rr) + E0 * E1[j]
  e00 <- Cov_BE[i, j] + E0^2
  ETT <- PI[i] * PI[j] * e11 + PI[i] * (1 - PI[j]) * e10 + (1 - PI[i]) * PI[j] * e01 + (1 - PI[i]) * (1 - PI[j]) * e00
  Cov_sp[i, j] <- Cov_sp[j, i] <- ETT - ET[i] * ET[j]
}
diag(Cov_sp) <- VarT
cat(sprintf("sparse mixture done (%.0f s)\n", as.numeric(Sys.time() - t0, units = "secs")))

# ---- convergence check of the pair quadrature on the most correlated pair --------------
ij <- which(abs(R) == max(abs(R[upper.tri(R)])), arr.ind = TRUE)[1, ]
conv <- abs(cov_pair(mu[ij[1]], s[ij[1]], mu[ij[2]], s[ij[2]], R[ij[1], ij[2]]) -
            cov_pair(mu[ij[1]], s[ij[1]], mu[ij[2]], s[ij[2]], R[ij[1], ij[2]], 200, 100))

# ---- the order-8 expansion with accurate coefficients (Victoria Pan's reference) ---------
hermite_norm <- function(u, K) { H <- matrix(0, length(u), K + 1); H[, 1] <- 1; if (K >= 1) H[, 2] <- u
  for (k in 1:(K - 1)) H[, k + 2] <- (u * H[, k + 1] - sqrt(k) * H[, k]) / sqrt(k + 1); H }
coefs <- function(m, sc, K = 8, panels = 200) { nd <- nodes_split(-LIM, LIM, -m / sc, panels)
  c(crossprod(hermite_norm(nd$x, K), g(m + sc * nd$x) * dnorm(nd$x) * nd$w)) }
B1 <- sapply(1:n, function(i) coefs(mu[i], s[i])); b0 <- coefs(0, 1)
Cov8_APE <- matrix(0, n, n); Cov8_BE <- matrix(0, n, n)
for (k in 1:8) { Cov8_APE <- Cov8_APE + outer(B1[k + 1, ], B1[k + 1, ]) * R^k; Cov8_BE <- Cov8_BE + b0[k + 1]^2 * M^k }

# ---- weights, as get_wts() computes them (positive part, normalized by the mean |w|) ----
wts <- function(Sigma, r) { w <- solve(Sigma, r); w <- w * (w > 0); w / mean(abs(w)) }
fixture <- list(
  source = inputs_path,
  transform = "g(x) = -2 log(2 Phi(-|x|)), the two-sided Fisher transform with df = 2 (g_GFisher_two(x, 2))",
  inputs = list(Bstar = Bstar, PI = PI, M = M),
  direct = list(BE_N = wts(Cov_BE, r_tilde), APE_N = wts(Cov_APE, r_tilde),
                BE_sparse = wts(Cov_BE_sp <- { S <- Cov_BE; diag(S) <- var_one(0, 1); S }, r_sp), APE_sparse = wts(Cov_sp, r_sp)),
  direct_cov = list(APE_N = Cov_APE, BE_N = Cov_BE, APE_sparse = Cov_sp),
  order8_accurate = list(BE_N = wts(Cov8_BE, r_tilde), APE_N = wts(Cov8_APE, r_tilde)),
  r_tilde = r_tilde, r_sparse = r_sp,
  pair_quadrature_convergence = conv,
  method = "nested kink-split composite Gauss-Legendre (10-point, 100/50 panels per side, [-10, 10]); the sparse diagonal is the mixture variance",
  generated = format(Sys.time(), "%Y-%m-%d %H:%M %Z"), R_version = R.version.string)
saveRDS(fixture, out_path)
cat(sprintf("fixture written: %s (pair quadrature convergence %.1e; %.0f s)\n", out_path, conv, as.numeric(Sys.time() - t0, units = "secs")))
