################################################################################
# assemble_pheno_covar.R
#
# Assemble a phenotype/covariate analysis bundle from a GDS sample order, an
# external phenotype/covariate table, and a PC matrix — the cohort-agnostic core
# of the single-variant / SNV-set "prepare pheno-covar" step.

#' Assemble a phenotype/covariate analysis bundle
#'
#' @description
#' Build the \code{(sample_id, Y, X, trait, ...)} bundle that
#' \code{\link{fit_null_model}} / \code{\link{marginal_scan}} and the GLOW SNV-set
#' pipeline consume, by merging three sources onto a single GDS sample order:
#' a raw outcome vector, an external phenotype/covariate table, and (optionally) a
#' principal-component matrix. Cohort-specific column names and recoding rules are
#' passed as arguments; the file I/O (opening the GDS, reading the CSV / PCs) is the
#' caller's responsibility, so this function is pure and testable.
#'
#' This is the packaged form of the per-cohort "prepare pheno-covariates" step.
#' Sample order is anchored to \code{sample_ids} (typically the GDS sample order);
#' the returned bundle preserves that order minus any dropped rows.
#'
#' @param sample_ids Character vector: the sample-order anchor (e.g. a GDS's
#'   \code{sample.id}). All other sources are aligned to this order.
#' @param outcome A list \code{list(values=, map=)} giving the raw outcome
#'   \code{values} (a vector aligned to \code{sample_ids}) and an optional
#'   \code{map}: a named vector (e.g. \code{c("1"=0L, "2"=1L)}; unmapped values
#'   become \code{NA}) or a function \code{function(values) -> Y}. If \code{map} is
#'   \code{NULL}, \code{values} are used as-is (coerced with \code{as.numeric}).
#' @param trait Either \code{"binary"} or \code{"continuous"} — recorded in the bundle.
#' @param pheno Optional data frame: the external phenotype/covariate table from
#'   which covariate columns are drawn. Merged onto \code{sample_ids} by \code{pheno_id_col}.
#' @param pheno_id_col Name of the column in \code{pheno} matching \code{sample_ids}
#'   (required when \code{covariates} draws from \code{pheno}).
#' @param covariates Named list defining the covariate (design) columns of \code{X},
#'   in order. Each element is either a string (a \code{pheno} column name, coerced
#'   with \code{as.numeric}) or a list \code{list(col=, map=)} where \code{map} is a
#'   named recode vector (unmapped -> \code{NA}) or a function. The list \emph{names}
#'   become the \code{X} column names.
#' @param pcs Optional list \code{list(scores=, id_col=, cols=)} appending PC columns
#'   to \code{X}: \code{scores} is a data frame / matrix of PC scores, \code{cols} the
#'   PC column names to take, and \code{id_col} (optional) a column of \code{scores}
#'   matched to \code{sample_ids} (if \code{NULL}, \code{scores} rows are assumed
#'   already in \code{sample_ids} order). The taken \code{cols} keep their names in \code{X}.
#' @param drop_incomplete If \code{TRUE} (default), after dropping \code{NA}-outcome
#'   samples, drop any remaining sample with an \code{NA} in \code{X}
#'   (\code{complete.cases}).
#' @param covar_set Optional label recorded in the bundle (\code{covar_set}).
#' @param verbose 0 = silent, >=1 = progress messages.
#'
#' @return A list with:
#' \describe{
#'   \item{sample_id}{retained sample IDs (subset of \code{sample_ids}, original order)}
#'   \item{Y}{the outcome for retained samples}
#'   \item{X}{numeric covariate matrix (covariate columns then PC columns)}
#'   \item{trait}{\code{"binary"} or \code{"continuous"}}
#'   \item{n_total}{\code{length(sample_ids)}}
#'   \item{n_excluded}{\code{n_total} minus the number retained}
#'   \item{covar_set}{the \code{covar_set} label}
#'   \item{covar_names}{\code{colnames(X)}}
#' }
#'
#' @seealso \code{\link{extract_pheno_covar_gds}} (reads pheno/covar from a GDS node),
#'   \code{\link{fit_null_model}}, \code{\link{marginal_scan}}
#'
#' @examples
#' set.seed(1)
#' ids <- paste0("S", 1:6)
#' pheno <- data.frame(
#'   iid = ids, status = c(2, 1, 2, 1, 2, -9),
#'   sex = c("F", "M", "F", "M", "F", "F"), age = c(50, 60, 55, 65, 58, 70),
#'   stringsAsFactors = FALSE)
#' pcs <- data.frame(iid = ids, PC1 = rnorm(6), PC2 = rnorm(6))
#' b <- assemble_pheno_covar(
#'   sample_ids = ids,
#'   outcome = list(values = pheno$status, map = c("1" = 0L, "2" = 1L)),
#'   trait = "binary",
#'   pheno = pheno, pheno_id_col = "iid",
#'   covariates = list(sex_numeric = list(col = "sex", map = c(F = 1, M = 0)),
#'                     age = "age"),
#'   pcs = list(scores = pcs, id_col = "iid", cols = c("PC1", "PC2")),
#'   covar_set = "demo")
#' str(b)
#'
#' @export
assemble_pheno_covar <- function(sample_ids,
                                 outcome,
                                 trait = c("binary", "continuous"),
                                 pheno = NULL,
                                 pheno_id_col = NULL,
                                 covariates = NULL,
                                 pcs = NULL,
                                 drop_incomplete = TRUE,
                                 covar_set = NULL,
                                 verbose = 1) {
  trait <- match.arg(trait)
  sample_ids <- as.character(sample_ids)
  n_total <- length(sample_ids)
  if (n_total < 1L) stop("`sample_ids` is empty.")

  # ---- Outcome -> Y (aligned to sample_ids) ----
  if (!is.list(outcome) || is.null(outcome$values))
    stop("`outcome` must be a list with at least `values` (aligned to sample_ids).")
  if (length(outcome$values) != n_total)
    stop("length(outcome$values) (", length(outcome$values),
         ") must equal length(sample_ids) (", n_total, ").")
  Y <- .recode_vec(outcome$values, outcome$map)

  # ---- Align the external pheno table onto sample_ids (left join, keep order) ----
  pheno_aligned <- NULL
  if (!is.null(pheno)) {
    pheno <- as.data.frame(pheno, stringsAsFactors = FALSE)  # robust to data.table
    if (is.null(pheno_id_col) || !pheno_id_col %in% names(pheno))
      stop("`pheno_id_col` must name a column of `pheno`.")
    idx <- match(sample_ids, as.character(pheno[[pheno_id_col]]))
    if (verbose >= 1)
      message(sprintf("  pheno match: %d / %d samples", sum(!is.na(idx)), n_total))
    pheno_aligned <- pheno[idx, , drop = FALSE]
  }

  # ---- Build covariate (design) columns, in order ----
  covar_cols <- list()
  if (!is.null(covariates)) {
    if (is.null(names(covariates)) || any(!nzchar(names(covariates))))
      stop("`covariates` must be a NAMED list (names become X column names).")
    for (nm in names(covariates)) {
      spec <- covariates[[nm]]
      col  <- if (is.character(spec)) spec else spec$col
      map  <- if (is.list(spec)) spec$map else NULL
      if (is.null(pheno_aligned) || !col %in% names(pheno_aligned))
        stop("Covariate '", nm, "' references column '", col,
             "' not found in `pheno`.")
      covar_cols[[nm]] <- .recode_vec(pheno_aligned[[col]], map)
    }
  }

  # ---- PC columns (matched to sample_ids) ----
  pc_mat <- NULL
  if (!is.null(pcs)) {
    if (is.null(pcs$scores) || is.null(pcs$cols))
      stop("`pcs` must be a list with `scores` and `cols`.")
    sc <- as.data.frame(pcs$scores, stringsAsFactors = FALSE)  # robust to data.table
    if (!is.null(pcs$id_col)) {
      pidx <- match(sample_ids, as.character(sc[[pcs$id_col]]))
      if (verbose >= 1)
        message(sprintf("  PC match: %d / %d samples", sum(!is.na(pidx)), n_total))
      sc <- sc[pidx, , drop = FALSE]
    } else if (nrow(sc) != n_total) {
      stop("`pcs$scores` has ", nrow(sc), " rows but no `id_col`; expected ",
           n_total, " rows aligned to sample_ids.")
    }
    missing_cols <- setdiff(pcs$cols, colnames(sc))
    if (length(missing_cols))
      stop("`pcs$cols` not found in scores: ", paste(missing_cols, collapse = ", "))
    pc_mat <- as.matrix(sc[, pcs$cols, drop = FALSE])
  }

  # ---- Assemble X (covariates then PCs) ----
  X_parts <- list()
  if (length(covar_cols))
    X_parts$covar <- do.call(cbind, c(covar_cols, list(deparse.level = 0)))
  if (!is.null(pc_mat)) X_parts$pc <- pc_mat
  if (!length(X_parts)) stop("No covariates and no PCs given; X would be empty.")
  X <- do.call(cbind, unname(X_parts))
  if (length(covar_cols)) colnames(X)[seq_along(covar_cols)] <- names(covar_cols)
  rownames(X) <- NULL   # sample identity lives in `sample_id`; keep X free of stray rownames

  # ---- Subset: drop NA-outcome, then (optionally) NA-covariate samples ----
  keep <- !is.na(Y)
  if (verbose >= 1)
    message(sprintf("  valid outcome: %d, missing: %d", sum(keep), sum(!keep)))
  Y <- Y[keep]; X <- X[keep, , drop = FALSE]; ids <- sample_ids[keep]
  if (drop_incomplete) {
    cc <- stats::complete.cases(X)
    if (any(!cc) && verbose >= 1)
      message(sprintf("  dropping %d samples with NA covariates", sum(!cc)))
    Y <- Y[cc]; X <- X[cc, , drop = FALSE]; ids <- ids[cc]
  }
  n_final <- length(Y)
  if (verbose >= 1)
    message(sprintf("  final: %d samples, %d covariates [%s]",
                    n_final, ncol(X), paste(colnames(X), collapse = ", ")))

  list(sample_id = ids, Y = Y, X = X, trait = trait,
       n_total = n_total, n_excluded = n_total - n_final,
       covar_set = covar_set, covar_names = colnames(X))
}

# Apply a recode `map` to a raw vector:
#   - named vector: out = unname(map[as.character(x)]) (unmapped -> NA), type from map
#   - function:     out = map(x)
#   - NULL:         out = as.numeric(x)
.recode_vec <- function(x, map = NULL) {
  if (is.null(map)) return(as.numeric(x))
  if (is.function(map)) return(map(x))
  if (is.null(names(map)) || any(!nzchar(names(map))))
    stop("A recode `map` must be a NAMED vector or a function.")
  unname(map[as.character(x)])
}
