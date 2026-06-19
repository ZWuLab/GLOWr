# ==============================================================================
# Per-Region Variant Drill-down Data Frame
# ==============================================================================
#
# Flattens one region's evidence packet into a per-variant drill-down
# data.frame: one row per post-Component-3-filter raw variant, with unit-level
# fields repeated across the variants composing the same post-collapse test
# unit and NA in the unit-level / weight columns for LD-removed variants. The
# single source of truth for the per-region drill-down schema (the chr-level
# developer viewer is migrated onto this).
#
# EXPORTED FUNCTIONS:
#   - glow_drilldown_df()   Region evidence packet -> flat per-variant df
#
# INTERNAL HELPER FUNCTIONS:
#   - `%||%`                File-local NULL-coalescing operator (not exported)


#################### EXPORTED MAIN FUNCTIONS ####################

#' Build a Per-Variant Drill-down Data Frame from a Region Evidence Packet
#'
#' Builds the flat per-region drill-down data.frame from one region's evidence
#' packet: one row per post-Component-3-filter raw variant, with unit-level
#' fields repeated across the variants composing the same post-collapse test
#' unit and \code{NA} in the unit-level + weight columns for LD-removed
#' variants.
#'
#' @details
#' Per-variant identity/metadata (\code{rsid}, \code{chr}, \code{pos},
#' \code{ref}, \code{alt}, \code{MAF}, \code{MAC}, \code{ld_kept}) is reported at
#' the per-variant level; the post-collapse aggregate, where one exists, is
#' explicitly prefixed \code{unit_*} (only \code{MAF}/\code{MAC} have one, hence
#' \code{unit_MAF} / \code{unit_MAC}; \code{unit_idx} / \code{unit_n_components}
#' / \code{unit_is_collapsed} describe the unit-membership mapping itself). Test
#' results and unit parameters (\code{Z_standard}, \code{pvalue_standard},
#' \code{Z_SPA}, \code{pvalue_SPA}, \code{B}, \code{PI}, \code{weight_<scheme>})
#' are unit-level -- one value per post-collapse test unit, repeated across that
#' unit's components, \code{NA} for LD-removed rows. The number of
#' \code{weight_<scheme>} columns follows \code{rownames(e$test_units$weights)}.
#'
#' @param e One region's evidence-packet entry, i.e. the value of
#'   \code{evidence_chr\{N\}.rds[[label]]} as produced by the per-chr/-chunk
#'   GLOW runner. Required sub-objects:
#'   \describe{
#'     \item{\code{e$post_filter}}{data.frame, one row per post-Component-3 raw
#'       variant. Required columns: \code{variant_id} (integer), \code{chr},
#'       \code{pos} (integer), \code{ref}, \code{alt} (character), \code{MAF},
#'       \code{MAC} (numeric), \code{ld_kept} (logical), \code{collapse_group}
#'       (integer index into \code{e$test_units}, \code{NA} when
#'       \code{ld_kept = FALSE}). Optional: \code{rsid} (character; NA-filled if
#'       absent).}
#'     \item{\code{e$test_units}}{list of length \code{n_after_collapse} with
#'       parallel per-unit fields. Required: \code{component_post_filter_idx}
#'       (list of integer vectors, one per unit), \code{is_collapsed} (logical),
#'       \code{MAF}, \code{B}, \code{PI}, \code{Z_standard},
#'       \code{pvalue_standard} (numeric), \code{weights}
#'       (\code{n_schemes x n_after_collapse} matrix; rownames are the scheme
#'       labels and become the suffixes of the output \code{weight_<scheme>}
#'       columns). Optional: \code{MAC} (NA-filled if absent), \code{Z_SPA},
#'       \code{pvalue_SPA} (if \code{NULL}, the output \code{Z_SPA} /
#'       \code{pvalue_SPA} columns are all \code{NA} -- these are
#'       binary-trait-only fields).}
#'   }
#'
#' @return A data.frame with \code{nrow == n_after_filter} rows. Column count is
#'   ~31 (exact width depends on the scheme count via
#'   \code{rownames(e$test_units$weights)}). Row order matches the order of
#'   \code{e$post_filter} (genomic order within the region).
#'
#' @export
glow_drilldown_df <- function(e) {
  pf <- e$post_filter
  tu <- e$test_units
  unit_idx <- pf$collapse_group   # NA for ld_kept = FALSE

  # Map a vector indexed by post-collapse position to the post-filter rows.
  # default_NA controls the fill value for LD-removed rows.
  .map_unit <- function(vec, unit_idx, default_NA) {
    out <- rep(default_NA, length(unit_idx))
    ok <- !is.na(unit_idx)
    out[ok] <- vec[unit_idx[ok]]
    out
  }

  unit_n_components <- .map_unit(lengths(tu$component_post_filter_idx),
                                 unit_idx, NA_integer_)
  unit_is_collapsed <- .map_unit(tu$is_collapsed, unit_idx, NA)
  unit_MAF <- .map_unit(tu$MAF, unit_idx, NA_real_)
  unit_MAC <- .map_unit(tu$MAC %||% rep(NA_real_, length(tu$MAF)),
                        unit_idx, NA_real_)
  B_col  <- .map_unit(tu$B,  unit_idx, NA_real_)
  PI_col <- .map_unit(tu$PI, unit_idx, NA_real_)
  Z_standard      <- .map_unit(tu$Z_standard,      unit_idx, NA_real_)
  pvalue_standard <- .map_unit(tu$pvalue_standard, unit_idx, NA_real_)
  if (!is.null(tu$Z_SPA)) {
    Z_SPA      <- .map_unit(tu$Z_SPA,      unit_idx, NA_real_)
    pvalue_SPA <- .map_unit(tu$pvalue_SPA, unit_idx, NA_real_)
  } else {
    Z_SPA      <- rep(NA_real_, length(unit_idx))
    pvalue_SPA <- rep(NA_real_, length(unit_idx))
  }

  # 13 weight columns, ordered as in the evidence packet's rownames
  schemes <- rownames(tu$weights)
  weights_df <- as.data.frame(
    lapply(seq_along(schemes), function(s)
      .map_unit(tu$weights[s, ], unit_idx, NA_real_)),
    col.names = paste0("weight_", schemes),
    stringsAsFactors = FALSE
  )

  cbind(
    data.frame(
      rsid              = if (!is.null(pf$rsid)) pf$rsid else NA_character_,
      chr               = pf$chr,
      pos               = pf$pos,
      ref               = pf$ref,
      alt               = pf$alt,
      MAF               = pf$MAF,
      MAC               = pf$MAC,
      ld_kept           = pf$ld_kept,
      unit_idx          = unit_idx,
      unit_n_components = unit_n_components,
      unit_is_collapsed = unit_is_collapsed,
      unit_MAF          = unit_MAF,
      unit_MAC          = unit_MAC,
      B                 = B_col,
      PI                = PI_col,
      Z_standard        = Z_standard,
      pvalue_standard   = pvalue_standard,
      Z_SPA             = Z_SPA,
      pvalue_SPA        = pvalue_SPA,
      stringsAsFactors  = FALSE
    ),
    weights_df
  )
}


#################### INTERNAL HELPER FUNCTIONS ####################

# Internal: NULL-coalescing operator (used by glow_drilldown_df()). Copied
# not otherwise provide a `%||%`. Not exported.
`%||%` <- function(a, b) if (is.null(a)) b else a
