# This file contains g-transformation functions for GFisher integration

#' Two-sided p-value transformation for GFisher
#'
#' Transforms Z-scores to chi-square statistics using two-sided p-values.
#' This transformation is essential for combining Z-scores using the Generalized
#' Fisher method.
#'
#' @param x Numeric vector of Z-scores
#' @param df Degrees of freedom for the chi-square distribution (typically 2 for
#'   two-sided tests)
#'
#' @return Numeric vector of transformed chi-square statistics
#'
#' @details
#' The transformation implements:
#' \deqn{g(x) = F_{\chi^2_{df}}^{-1}(\log(2) + \log(\bar{\Phi}(|x|)))}
#' where \eqn{\bar{\Phi}(|x|)} is the upper tail probability of the standard normal
#' distribution evaluated at \eqn{|x|}.
#'
#' The factor \eqn{\log(2)} accounts for the two-sided p-value calculation:
#' \eqn{p = 2\bar{\Phi}(|x|)}.
#'
#' All calculations are performed in log space for numerical stability.
#'
#' Computational complexity: O(n) where n is the length of x
#'
#' @note This function must maintain exact numerical fidelity with the legacy
#'   implementation for correct Fisher test results. Do not simplify or modify
#'   the transformation.
#'
#' @references
#' Van Lishout, F., et al. (2013). An efficient algorithm to perform multiple
#' testing in epistasis screening. BMC Bioinformatics, 14, 138.
#'
#' @noRd
#' @keywords internal
g_GFisher_two <- function(x, df) {
  # Transform Z-scores to chi-square statistics using two-sided p-values
  # p = 2 * Phi_bar(|x|), where Phi_bar is the upper tail of standard normal
  # Use log.p = TRUE throughout for numerical stability
  qchisq(
    log(2) + pnorm(abs(x), lower.tail = FALSE, log.p = TRUE),
    df = df,
    lower.tail = FALSE,
    log.p = TRUE
  )
}


#' One-sided p-value transformation for GFisher
#'
#' Transforms Z-scores to chi-square statistics using one-sided p-values.
#' This transformation is essential for combining Z-scores using the Generalized
#' Fisher method when directional information is important.
#'
#' @param x Numeric vector of Z-scores
#' @param df Degrees of freedom for the chi-square distribution (typically 2 for
#'   one-sided tests)
#'
#' @return Numeric vector of transformed chi-square statistics
#'
#' @details
#' The transformation implements:
#' \deqn{g(x) = F_{\chi^2_{df}}^{-1}(\log(\bar{\Phi}(x)))}
#' where \eqn{\bar{\Phi}(x)} is the upper tail probability of the standard normal
#' distribution evaluated at \eqn{x}.
#'
#' Unlike the two-sided version, this does not use the absolute value of x,
#' preserving directional information. The p-value is simply \eqn{p = \bar{\Phi}(x)}.
#'
#' All calculations are performed in log space for numerical stability.
#'
#' Computational complexity: O(n) where n is the length of x
#'
#' @note This function must maintain exact numerical fidelity with the legacy
#'   implementation for correct Fisher test results. Do not simplify or modify
#'   the transformation.
#'
#' @references
#' Van Lishout, F., et al. (2013). An efficient algorithm to perform multiple
#' testing in epistasis screening. BMC Bioinformatics, 14, 138.
#'
#' @noRd
#' @keywords internal
g_GFisher_one <- function(x, df) {
  # Transform Z-scores to chi-square statistics using one-sided p-values
  # p = Phi_bar(x), where Phi_bar is the upper tail of standard normal
  # Use log.p = TRUE throughout for numerical stability
  qchisq(
    pnorm(x, lower.tail = FALSE, log.p = TRUE),
    df = df,
    lower.tail = FALSE,
    log.p = TRUE
  )
}


#' Identity transformation for Burden test
#'
#' Returns the input Z-scores unchanged. This is used for Burden/Laptik tests
#' where no transformation is needed.
#'
#' @param x Numeric vector of Z-scores
#'
#' @return Numeric vector identical to input x
#'
#' @details
#' The identity transformation: \eqn{g(x) = x}
#'
#' This is used when combining statistics that are already on the appropriate
#' scale for the test being performed (e.g., Burden test).
#'
#' Computational complexity: O(1)
#'
#' @noRd
#' @keywords internal
g_Burden <- function(x) {
  # Identity transformation: g(x) = x
  x
}


#' Wrapper for GFisher transformation functions
#'
#' Convenient wrapper that calls either \code{g_GFisher_two} or \code{g_GFisher_one}
#' based on the specified p-value type.
#'
#' @param x Numeric vector of Z-scores
#' @param df Degrees of freedom for the chi-square distribution
#' @param p.type Character string specifying p-value type: "two" for two-sided
#'   (default) or "one" for one-sided
#'
#' @return Numeric vector of transformed chi-square statistics
#'
#' @details
#' This wrapper provides a unified interface for both one-sided and two-sided
#' transformations. It routes to the appropriate transformation function based
#' on \code{p.type}.
#'
#' For \code{p.type = "two"}: calls \code{g_GFisher_two(x, df)}
#' For \code{p.type = "one"}: calls \code{g_GFisher_one(x, df)}
#'
#' Using this wrapper can be more convenient than directly calling the specific
#' transformation functions when the p-value type is determined programmatically.
#'
#' Computational complexity: O(n) where n is the length of x
#'
#' @note If \code{p.type} is neither "one" nor "two", the function will return
#'   NULL (no error checking implemented to match legacy behavior).
#'
#' @noRd
#' @keywords internal
g_GFisher <- function(x, df, p.type = "two") {
  # Wrapper function that routes to appropriate transformation
  if (p.type == "two") {
    g_GFisher_two(x, df)
  } else if (p.type == "one") {
    g_GFisher_one(x, df)
  }
}
