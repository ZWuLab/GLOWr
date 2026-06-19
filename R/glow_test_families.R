# ==============================================================================
# GLOW Test Family Classification
# ==============================================================================
#
# Maps GLOW (and, via overrides, non-GLOW) test column names to test families
# (SKAT, Burden, Fisher, Omni, ...) and groups a vector of test names into a
# panel-ready named list of families. These classifiers sit beside the test-spec
# machinery (glow_test_specs.R / default_test_specs()) and drive per-family QQ
# panels and family-aware reporting in the summary layer.
#
# EXPORTED FUNCTIONS:
#   - glow_family_of()                Vectorised test name -> family label
#   - group_glow_tests_by_family()    Split test names into a named family list


#################### EXPORTED MAIN FUNCTIONS ####################

#' Classify GLOW Test Names into Test Families
#'
#' Vectorised classifier mapping test column names to test families. By
#' convention a name of the form \code{GLOW_<Family>_<scheme>...} resolves to
#' \code{<Family>} (the first token after the \code{GLOW_} prefix), e.g.
#' \code{GLOW_SKAT_BE_N} -> \code{"SKAT"}, \code{GLOW_Burden_equ} ->
#' \code{"Burden"}, \code{GLOW_Fisher_APR_sparse} -> \code{"Fisher"}. The three
#' omnibus columns (\code{GLOW_BSF_Omni}, \code{GLOW_SNV_CCT}, \code{GLOW_Omni})
#' are bucketed into the synthetic \code{"Omni"} family.
#'
#' @details
#' Non-\code{GLOW_*} names resolve to \code{NA} unless mapped by
#' \code{family_overrides}, so callers can detect typos (an \code{NA} return for
#' an expected test name flags a name that matched no rule). The function is
#' fully vectorised; complexity is \eqn{O(n)} in \code{length(names)}.
#'
#' @param names Character vector of test column names.
#' @param family_overrides Optional named character vector mapping a non-GLOW
#'   test name (names of the vector) to a family label (values). Example:
#'   \code{c(STAAR_O = "Omni")} routes \code{STAAR_O} into the \code{Omni}
#'   family. Default \code{NULL} preserves the base contract (non-GLOW names
#'   resolve to \code{NA}).
#'
#' @return Character vector of family labels, the same length as \code{names}.
#'
#' @examples
#' glow_family_of(c("GLOW_SKAT_BE_N", "GLOW_Burden_equ", "GLOW_Omni"))
#' glow_family_of("STAAR_O", family_overrides = c(STAAR_O = "Omni"))
#'
#' @seealso \code{\link{group_glow_tests_by_family}}
#' @export
glow_family_of <- function(names, family_overrides = NULL) {
  omni_set <- c("GLOW_BSF_Omni", "GLOW_SNV_CCT", "GLOW_Omni")
  out <- rep(NA_character_, length(names))
  is_omni <- names %in% omni_set
  out[is_omni] <- "Omni"
  is_other_glow <- grepl("^GLOW_", names) & !is_omni
  if (any(is_other_glow)) {
    parts <- strsplit(sub("^GLOW_", "", names[is_other_glow]), "_",
                      fixed = TRUE)
    out[is_other_glow] <- vapply(parts, `[`, character(1L), 1L)
  }
  if (!is.null(family_overrides)) {
    if (is.null(names(family_overrides)) ||
        any(!nzchar(names(family_overrides)))) {
      stop("`family_overrides` must be a named character vector ",
           "(names = test names, values = family labels)")
    }
    hit <- match(names, names(family_overrides))
    out[!is.na(hit)] <- unname(family_overrides[hit[!is.na(hit)]])
  }
  out
}


#' Group GLOW Test Names into a Named List of Families
#'
#' Groups a vector of test column names into a named list, one entry per test
#' family, suitable for driving per-family QQ panels.
#'
#' @details
#' Family ordering: families in \code{priority} appear first in the listed order
#' (only those actually present), then any remaining non-Omni families in
#' alphabetical order, then the synthetic \code{"Omni"} family last when
#' \code{omni_last = TRUE} (so it lands in a sensible position on a square panel
#' grid). Within each family, columns retain their input order. The function
#' errors if any input column does not resolve to a family (via the
#' \code{GLOW_*} convention or \code{family_overrides}). Complexity is
#' \eqn{O(n)} in \code{length(tests)}.
#'
#' @param tests Character vector of test column names (typically
#'   \code{grep("^GLOW_", names(df), value = TRUE)}, optionally augmented with
#'   non-GLOW tests such as \code{"STAAR_O"}).
#' @param priority Character vector giving the family-name priority for the
#'   leading panels. Default \code{c("SKAT", "Burden", "Fisher")}.
#' @param omni_last Logical. If \code{TRUE} (default), the \code{"Omni"} family
#'   is forced last.
#' @param family_overrides Optional named character vector forwarded to
#'   \code{\link{glow_family_of}}. Example: \code{c(STAAR_O = "Omni")}. Default
#'   \code{NULL} requires every input to be a recognised \code{GLOW_*} name.
#'
#' @return Named list of character vectors. Names are family labels; each
#'   element holds the column names of that family.
#'
#' @examples
#' group_glow_tests_by_family(
#'   c("GLOW_SKAT_BE_N", "GLOW_Burden_equ", "GLOW_Fisher_equ", "GLOW_Omni"))
#'
#' @seealso \code{\link{glow_family_of}}
#' @export
group_glow_tests_by_family <- function(
    tests,
    priority = c("SKAT", "Burden", "Fisher"),
    omni_last = TRUE,
    family_overrides = NULL) {
  fam <- glow_family_of(tests, family_overrides = family_overrides)
  bad <- is.na(fam)
  if (any(bad)) {
    stop("Unrecognised test column name(s) (not GLOW_* and not in ",
         "family_overrides): ",
         paste(tests[bad], collapse = ", "))
  }
  present <- unique(fam)
  leading  <- intersect(priority, present)
  trailing <- setdiff(present, leading)
  if (omni_last) trailing <- setdiff(trailing, "Omni")
  trailing <- sort(trailing)
  ordered <- c(leading, trailing)
  if (omni_last && "Omni" %in% present) ordered <- c(ordered, "Omni")

  out <- lapply(ordered, function(f) tests[fam == f])
  names(out) <- ordered
  out
}
