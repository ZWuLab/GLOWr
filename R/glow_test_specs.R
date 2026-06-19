# ==============================================================================
# GLOW Test Specification Functions
# ==============================================================================
#
# Defines the test specification data structure for the GLOW extensible
# composition layer. Test specs describe statistical test families (Burden,
# SKAT, Fisher) with their transformation functions, degrees of freedom,
# weight schemes, and p-value types.
#
# EXPORTED FUNCTIONS:
#   - default_test_specs()       Return the default 3-family test spec list
#
# INTERNAL HELPER FUNCTIONS:
#   - .default_equal_weights()   Generate equal-weight vector
#   - validate_test_spec()       Validate and apply defaults to a test spec
#   - make_test_spec()           Convenience constructor for test specs


#################### EXPORTED MAIN FUNCTIONS ####################

#' Default Test Specifications for GLOW Omnibus Test
#'
#' Returns the default three-family test specification list used by the GLOW
#' omnibus test. The families are SKAT (df=1), Burden (df=Inf), and Fisher
#' (df=2), matching the stacking order used by \code{run_bsf_tests}.
#'
#' @details
#' Each test spec is a list with:
#' \describe{
#'   \item{family}{Character string identifying the test family}
#'   \item{g}{Transformation function applied to Z-scores}
#'   \item{df}{Degrees of freedom for the test}
#'   \item{weight_fns}{Named list of weight-generating functions}
#'   \item{p.type}{P-value type: "two" for two-sided}
#' }
#'
#' The default weight functions include \code{optimal} (via
#' \code{\link{Optimal_Weights_M}}) and \code{equal} (unit weights).
#'
#' Family order is SKAT, Burden, Fisher to match the stacking order
#' used by \code{run_bsf_tests}.
#'
#' **Computational Complexity**: O(1) -- returns a static list.
#'
#' @return A list of 3 test spec lists, one per family (SKAT, Burden, Fisher).
#'
#' @note
#' Each spec's \code{g} function is used for **weight computation** via
#' \code{Optimal_Weights_M}. The \code{df} value determines the **test
#' statistic transformation**, applied via \code{g_GFisher(x, df=...)}.
#' For the built-in families, these are consistent: SKAT (\code{g(x)=x^2},
#' \code{df=1}), Burden (\code{g(x)=x}, \code{df=Inf}), and Fisher
#' (\code{g(x)=g_GFisher_two(x,2)}, \code{df=2}). Custom specs must
#' maintain this consistency; see \code{make_test_spec} for details.
#'
#' @references
#' Zhang, H., Liu, M., Landers, J. E., and Wu, Z. Integrated Weighted Association
#' Test with Application to Genetic Association Studies. Annals of Applied
#' Statistics (in revision).
#'
#' @examples
#' specs <- default_test_specs()
#' length(specs)  # 3
#' specs[[1]]$family  # "SKAT"
#' specs[[2]]$family  # "Burden"
#' specs[[3]]$family  # "Fisher"
#'
#' @seealso \code{\link{Optimal_Weights_M}} for the optimal weight function
#'
#' @export
default_test_specs <- function() {
  list(
    # --- SKAT: df=1, g(x) = x^2 ---
    list(
      family = "SKAT",
      g = function(x) x^2,
      df = 1,
      weight_fns = list(
        optimal = Optimal_Weights_M,
        equal = .default_equal_weights
      ),
      p.type = "two"
    ),
    # --- Burden: df=Inf, g(x) = x ---
    list(
      family = "Burden",
      g = function(x) x,
      df = Inf,
      weight_fns = list(
        optimal = Optimal_Weights_M,
        equal = .default_equal_weights
      ),
      p.type = "two"
    ),
    # --- Fisher: df=2, g(x) = g_GFisher_two(x, df) ---
    list(
      family = "Fisher",
      g = function(x, df = 2) g_GFisher_two(x, df),
      df = 2,
      weight_fns = list(
        optimal = Optimal_Weights_M,
        equal = .default_equal_weights
      ),
      p.type = "two"
    )
  )
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' Generate equal (unit) weights
#'
#' Returns a named list with a single element \code{equ} containing a vector
#' of ones of length \code{p}. Used as the default equal-weight scheme in
#' test specs.
#'
#' @param p Integer, number of variants.
#'
#' @return A named list: \code{list(equ = rep(1, p))}.
#'
#' @keywords internal
#' @noRd
.default_equal_weights <- function(p) {
  list(equ = rep(1, p))
}


#' Validate a test specification
#'
#' Checks that a test spec has all required fields (\code{family}, \code{g},
#' \code{df}) and applies defaults for optional fields (\code{weight_fns},
#' \code{p.type}). Returns the validated and defaulted spec.
#'
#' @param spec A list representing a test specification.
#'
#' @return The validated spec list with defaults applied.
#'
#' @details
#' Validation rules:
#' \itemize{
#'   \item \code{family}: must be a non-empty character string
#'   \item \code{g}: must be a function
#'   \item \code{df}: must be a numeric scalar > 0 (Inf is valid)
#'   \item \code{weight_fns} (optional): must be a named list of functions;
#'     defaults to optimal + equal
#'   \item \code{p.type} (optional): defaults to "two"
#' }
#'
#' @keywords internal
#' @noRd
validate_test_spec <- function(spec) {
  # --- Required: family ---
  if (is.null(spec$family) || !is.character(spec$family) ||
      length(spec$family) != 1 || nchar(spec$family) == 0) {
    stop("Test spec 'family' must be a non-empty character string")
  }

  # --- Required: g ---
  if (is.null(spec$g) || !is.function(spec$g)) {
    stop("Test spec 'g' must be a function (family: ", spec$family, ")")
  }

  # --- Required: df ---
  if (is.null(spec$df) || !is.numeric(spec$df) || length(spec$df) != 1 ||
      is.na(spec$df) || spec$df <= 0) {
    stop("Test spec 'df' must be a numeric scalar > 0 (family: ",
         spec$family, ")")
  }

  # --- Optional: weight_fns (default: optimal + equal) ---
  if (is.null(spec$weight_fns)) {
    spec$weight_fns <- list(
      optimal = Optimal_Weights_M,
      equal = .default_equal_weights
    )
  } else {
    # Must be a named list of functions
    if (!is.list(spec$weight_fns)) {
      stop("Test spec 'weight_fns' must be a list (family: ",
           spec$family, ")")
    }
    wfn_names <- names(spec$weight_fns)
    if (is.null(wfn_names) || any(wfn_names == "")) {
      stop("Test spec 'weight_fns' must be a named list (family: ",
           spec$family, ")")
    }
    for (nm in wfn_names) {
      if (!is.function(spec$weight_fns[[nm]])) {
        stop("weight_fns[['", nm, "']] must be a function (family: ",
             spec$family, ")")
      }
    }
  }

  # --- Optional: p.type (default: "two") ---
  if (is.null(spec$p.type)) {
    spec$p.type <- "two"
  } else if (!spec$p.type %in% c("one", "two")) {
    stop("Test spec 'p.type' must be \"one\" or \"two\" (family: ",
         spec$family, ")")
  }

  spec
}


#' Convenience constructor for a test specification
#'
#' Creates a test spec list and validates it via \code{validate_test_spec()}.
#'
#' @param family Character string identifying the test family.
#' @param g Transformation function applied to Z-scores. Used for weight
#'   computation via \code{Optimal_Weights_M}.
#' @param df Numeric scalar, degrees of freedom (use Inf for Burden). Determines
#'   the test statistic transformation, applied via
#'   \code{g_GFisher(x, df=..., p.type=...)}.
#' @param weight_fns Named list of weight-generating functions. Default: NULL
#'   (uses optimal + equal).
#' @param p.type Character string for p-value type. Default: NULL (uses "two").
#'
#' @note
#' Custom families must ensure \code{g} and \code{df} are consistent:
#' \code{g(x)} should equal \code{g_GFisher(x, df, p.type)} for correct
#' results. For the built-in families (SKAT: \code{g(x) = x^2, df = 1};
#' Burden: \code{g(x) = x, df = Inf}; Fisher: \code{g(x) = g_GFisher_two(x, 2),
#' df = 2}), this consistency holds. See \code{run_bsf_tests} for details.
#'
#' @return A validated test spec list.
#'
#' @keywords internal
#' @noRd
make_test_spec <- function(family, g, df, weight_fns = NULL, p.type = NULL) {
  spec <- list(family = family, g = g, df = df)
  if (!is.null(weight_fns)) {
    spec$weight_fns <- weight_fns
  }
  if (!is.null(p.type)) {
    spec$p.type <- p.type
  }
  validate_test_spec(spec)
}
