# ==============================================================================
# File: optimal_weights_quadrature.R
# Purpose: The numerical integration behind the optimal weights of Optimal_Weights_M():
#          a composite Gauss-Legendre rule on [-12, 12] against the standard normal
#          density, cut where the argument of the transformation g is zero (the kink of
#          the two-sided transforms), the normalized probabilists' Hermite functions by
#          their three-term recurrence, and the Hermite moments of g(mu + s U), U ~ N(0,1),
#          from one evaluation of g per distinct (mu, s) pair. Mehler's expansion
#          Cov{g(mu_i + s_i U_i), g(mu_j + s_j U_j)} = sum_k b_k(i) b_k(j) r_ij^k is then
#          one matrix sum (.mehler_sum). Design and verification:
#
# INTERNAL HELPERS:
#   - .gl_rule()           Gauss-Legendre nodes and weights on [-1, 1], cached
#   - .normal_nodes()      the composite rule's nodes and weights, cut at a point
#   - .hermite_functions() he_k(u) = He_k(u) / sqrt(k!), k = 0..K, by recurrence
#   - .hermite_moments()   b_0..b_K, E[g^2] and Var[g] of g(mu + s U) for one (mu, s)
#   - .hermite_table()     the moments for vectors mu, s, each distinct pair once
#   - .mehler_sum()        sum over the orders of c_k(i) d_k(j) M_ij^k

#################### INTERNAL HELPER FUNCTIONS ####################

# Package-local environment for the cached Gauss-Legendre rules
.ow_env <- new.env(parent = emptyenv())

#' Gauss-Legendre nodes and weights on the interval \eqn{[-1, 1]}
#'
#' Golub-Welsch: the nodes are the eigenvalues of the symmetric tridiagonal Jacobi matrix
#' of the Legendre recurrence, the weights twice the squared first components of its
#' eigenvectors. Computed once per \code{m} and kept in a package environment.
#'
#' @param m Integer, the number of nodes (default 10).
#' @return A list with the nodes \code{x} (increasing) and the weights \code{w}.
#' @keywords internal
#' @noRd
.gl_rule <- function(m = 10L) {
  key <- paste0("gl", m)
  rule <- .ow_env[[key]]
  if (!is.null(rule)) return(rule)
  i <- seq_len(m - 1L)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, m, m)
  J[cbind(i, i + 1L)] <- b
  J[cbind(i + 1L, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  rule <- list(x = e$values[o], w = 2 * e$vectors[1L, o]^2)
  assign(key, rule, envir = .ow_env)
  rule
}

#' Nodes and weights of the composite rule for integrals against the normal density
#'
#' The domain \code{[-lim, lim]} (the standard normal mass outside \code{[-12, 12]} is
#' \eqn{3.6 \times 10^{-33}}) is cut into panels of width at most \code{width}, with a
#' cut at \code{cut} when it lies inside, so that a kink of the integrand at \code{cut}
#' falls on a panel boundary. Each panel carries the \code{m}-point Gauss-Legendre rule.
#' The density is not included: the caller multiplies by \code{dnorm(x)}.
#'
#' @param cut Numeric, the point to cut at, or \code{NA} for no cut.
#' @param lim Numeric, the half-width of the domain (default 12).
#' @param width Numeric, the largest panel width (default 0.1).
#' @param m Integer, nodes per panel (default 10).
#' @return A list with the nodes \code{x} and the weights \code{w}.
#' @keywords internal
#' @noRd
.normal_nodes <- function(cut = NA_real_, lim = 12, width = 0.1, m = 10L) {
  rule <- .gl_rule(m)
  panel <- function(a, b) {
    k <- max(1L, ceiling((b - a) / width))
    edges <- seq(a, b, length.out = k + 1L)
    h <- diff(edges) / 2
    c0 <- (edges[-1L] + edges[-length(edges)]) / 2
    list(x = c(outer(rule$x, h) + rep(c0, each = m)), w = c(outer(rule$w, h)))
  }
  if (is.finite(cut) && cut > -lim && cut < lim) {
    lo <- panel(-lim, cut)
    hi <- panel(cut, lim)
    list(x = c(lo$x, hi$x), w = c(lo$w, hi$w))
  } else {
    panel(-lim, lim)
  }
}

#' Normalized probabilists' Hermite functions at the nodes
#'
#' \eqn{he_k(u) = He_k(u) / \sqrt{k!}} for \eqn{k = 0, \dots, K}, orthonormal under the
#' standard normal density, by the forward-stable recurrence
#' \eqn{he_{k+1}(u) = \{u\,he_k(u) - \sqrt{k}\,he_{k-1}(u)\} / \sqrt{k + 1}}.
#'
#' @param u Numeric vector, the evaluation points.
#' @param K Integer, the largest order (0 is valid).
#' @return A \code{length(u)} by \code{K + 1} matrix; column \code{k + 1} holds order
#'   \code{k}.
#' @keywords internal
#' @noRd
.hermite_functions <- function(u, K) {
  H <- matrix(0, length(u), K + 1L)
  H[, 1L] <- 1
  if (K >= 1L) H[, 2L] <- u
  if (K >= 2L) {
    for (k in seq_len(K - 1L)) {
      H[, k + 2L] <- (u * H[, k + 1L] - sqrt(k) * H[, k]) / sqrt(k + 1)
    }
  }
  H
}

#' Hermite moments of g(mu + s U), U ~ N(0, 1), for one (mu, s)
#'
#' From one evaluation of \code{g} on the nodes of \code{.normal_nodes()}, cut at
#' \eqn{u_0 = -\mu / s} where the argument of \code{g} is zero (the kink of the two-sided
#' transforms, which depend on \eqn{|x|}; a smooth transform is unaffected by the cut):
#' the coefficients \eqn{b_k = E[g(\mu + sU)\,he_k(U)]} for \eqn{k = 0, \dots, K} (so
#' \eqn{b_0 = E[g]}), the second moment \eqn{E[g^2]} and the centered variance
#' \eqn{E[(g - E g)^2]}, which does not cancel a large mean against a large second moment.
#'
#' The rule is validated for the Fisher and polynomial transforms over the scales GLOW
#' meets (design section 4.3). It assumes \eqn{g(\mu + sU)} square-integrable and smooth
#' away from the cut: a transform with another kink, a cusp, a discontinuity or a growth
#' like \eqn{e^{8u}} needs its own convergence check.
#'
#' @param g Function, the transformation, vectorized over its argument.
#' @param mu,s Numeric scalars, the mean and the scale (\code{s > 0}).
#' @param K Integer, the largest Hermite order (0 is valid).
#' @return A list: \code{b} (numeric, length \code{K + 1}), \code{Eg2}, \code{Var}.
#' @keywords internal
#' @noRd
.hermite_moments <- function(g, mu, s, K) {
  if (!is.finite(mu) || !is.finite(s) || s <= 0) {
    stop("The mean must be finite and the scale finite and positive (got mu = ",
         format(mu), ", s = ", format(s), ").", call. = FALSE)
  }
  nd <- .normal_nodes(cut = -mu / s)
  x <- mu + s * nd$x
  gx <- g(x)
  if (length(gx) != length(x)) {
    stop("The transformation g must return one value per input (", length(x),
         " inputs, ", length(gx), " values).", call. = FALSE)
  }
  bad <- which(!is.finite(gx))
  if (length(bad)) {
    stop("The transformation g returned a non-finite value at x = ", format(x[bad[1L]]),
         " (mu = ", format(mu), ", s = ", format(s), ").", call. = FALSE)
  }
  wphi <- stats::dnorm(nd$x) * nd$w
  H <- .hermite_functions(nd$x, K)
  b <- as.vector(crossprod(H, gx * wphi))
  Eg2 <- sum(gx^2 * wphi)
  Var <- sum((gx - b[1L])^2 * wphi)
  if (!all(is.finite(b)) || !is.finite(Eg2) || !is.finite(Var)) {
    stop("A moment of g(mu + s U) is not finite (mu = ", format(mu), ", s = ", format(s),
         "): the transformation grows too fast for the fixed quadrature domain.",
         call. = FALSE)
  }
  list(b = b, Eg2 = Eg2, Var = Var)
}

#' Hermite moments for vectors of means and scales, each distinct pair once
#'
#' Pairs are matched exactly (both doubles at 17 significant digits), so the null pair
#' (0, 1) repeated over every variant costs one evaluation of \code{g}.
#'
#' @param g Function, the transformation.
#' @param mu,s Numeric vectors of the same length.
#' @param K Integer, the largest Hermite order.
#' @return A list: \code{b} (a \code{K + 1} by \code{length(mu)} matrix, row \code{k + 1}
#'   holding order \code{k}), \code{Eg2} and \code{Var} (vectors), indexed like \code{mu}.
#' @keywords internal
#' @noRd
.hermite_table <- function(g, mu, s, K) {
  n <- length(mu)
  if (length(s) != n) stop("mu and s must have the same length", call. = FALSE)
  key <- paste(format(mu, digits = 17L), format(s, digits = 17L))
  uk <- unique(key)
  idx <- match(key, uk)
  first <- match(uk, key)
  mom <- lapply(first, function(i) .hermite_moments(g, mu[i], s[i], K))
  B <- matrix(unlist(lapply(mom, `[[`, "b")), nrow = K + 1L)
  list(b = B[, idx, drop = FALSE],
       Eg2 = vapply(mom, `[[`, 0, "Eg2")[idx],
       Var = vapply(mom, `[[`, 0, "Var")[idx])
}

#' Mehler's sum over the orders
#'
#' \eqn{S_{ij} = \sum_{k \in ORD} C_{k+1,i} D_{k+1,j} M_{ij}^k}, the powers of \code{M}
#' accumulated elementwise. With \code{C = D} the coefficient matrices of
#' \code{.hermite_table()} and \code{M} a correlation matrix, this is the truncated
#' covariance of the transformed variables; with \code{M = I} only the diagonal survives.
#'
#' @param C,D Numeric matrices, \code{max(ORD) + 1} or more rows, \code{n} columns.
#' @param M Numeric \code{n} by \code{n} matrix of correlations.
#' @param ORD Integer vector, the orders summed (positive integers).
#' @return An \code{n} by \code{n} matrix.
#' @keywords internal
#' @noRd
.mehler_sum <- function(C, D, M, ORD) {
  n <- ncol(C)
  out <- matrix(0, n, n)
  Mk <- matrix(1, n, n)
  K <- max(ORD)
  for (k in seq_len(K)) {
    Mk <- Mk * M
    if (k %in% ORD) out <- out + tcrossprod(C[k + 1L, ], D[k + 1L, ]) * Mk
  }
  out
}
