# ==============================================================================
# Marginal Single-Variant Analysis Functions
# ==============================================================================
#
# Functions for genome-wide marginal association scanning from GDS files.
#
# EXPORTED FUNCTIONS:
#   - extract_pheno_covar_gds()  Phenotype/covariate extraction from GDS
#   - marginal_scan()            Chunked marginal scan with SPA support
#   - annotate_gds_marginal()    Write results to GDS annotation nodes
#
# INTERNAL HELPERS:
#   - .precompute_scan_components()  Extract components from null model
#   - .process_chunk_standard()      Standard asymptotic Z per chunk
#   - .process_chunk_SPA()           SPA-adjusted Z per chunk
#   - .impute_geno()                 Missing genotype imputation + MAF/MAC


#################### EXPORTED MAIN FUNCTIONS ####################

#' Marginal Single-Variant Association Scan
#'
#' Performs genome-wide marginal association analysis from a GDS file with
#' chunked processing for memory efficiency. Supports standard asymptotic
#' and saddlepoint approximation (SPA) methods.
#'
#' @param gds_file Character path to GDS file.
#' @param null_model A \code{glow_null_model} object from \code{fit_null_model()}.
#' @param sample_id Optional character vector of sample IDs to override
#'   \code{null_model$sample_id} for GDS sample alignment. Default: NULL
#'   (uses null_model$sample_id).
#' @param use_SPA SPA control: NULL (auto: TRUE for binary, FALSE for
#'   continuous), TRUE (force SPA, binary only), FALSE (standard only), or
#'   "both" (equivalent to TRUE). Default: NULL.
#' @param chunk_size Integer number of variants per processing chunk.
#'   Default: 2000L.
#' @param mac_cutoff Integer minimum minor allele count. Variants with
#'   MAC < mac_cutoff are excluded from results. Default: 1L.
#' @param missing_imputation Character method for missing genotypes: "mean"
#'   or "zero". Default: "mean".
#' @param output_csv Optional path to write results as CSV. Default: NULL.
#' @param verbose Integer verbosity level (0=silent, 1=progress, 2=detailed).
#'
#' @return A data.frame with class "glow_marginal_scan" containing columns:
#'   chr, pos, ref, alt, variant_id, MAF, MAC, score, var_score, Z, pvalue.
#'   When SPA is active, additional columns Z_SPA and pvalue_SPA are included.
#'
#' @details
#' For each chunk of variants, genotype dosages are read from the GDS file
#' using SeqArray, missing values imputed (mean or zero), and marginal score
#' statistics computed from the pre-fitted null model.
#'
#' Score statistic: \eqn{S_j = G_j^T (Y - \mu_0)} where \eqn{\mu_0} are
#' fitted values under the null. Variance is computed using the hat-matrix
#' projection to account for covariates. The Z-statistic is
#' \eqn{Z_j = S_j / \sqrt{Var(S_j)}}.
#'
#' When SPA is enabled (recommended for binary traits with case-control
#' imbalance), \code{SPAtest::ScoreTest_SPA()} is called per chunk and
#' SPA-corrected p-values are stored in the pvalue_SPA column.
#'
#' Complexity: O(n * m) per chunk where n is sample count and m is
#' chunk_size. Memory per chunk is O(n * chunk_size). Total runtime is
#' linear in the number of variants.
#'
#' @examples
#' \dontrun{
#' pheno <- extract_pheno_covar_gds("chr22.gds", covar_names = "sex")
#' nm    <- fit_null_model(pheno$X, pheno$Y, trait = "binary",
#'                         sample_id = pheno$sample_id)
#' results <- marginal_scan("chr22.gds", nm, use_SPA = FALSE, verbose = 1)
#' head(results)
#' }
#'
#' @seealso \code{\link{extract_pheno_covar_gds}}, \code{\link{fit_null_model}},
#'   \code{\link{getZ_marg_score}}, \code{\link{annotate_gds_marginal}}
#'
#' @references
#' Dey, R., Schmidt, E. M., Abecasis, G. R., and Lee, S. (2017). A fast and accurate
#' algorithm to test for binary phenotypes and its application to PheWAS. American
#' Journal of Human Genetics, 101(1), 37-49. doi:10.1016/j.ajhg.2017.05.014
#'
#' @export
marginal_scan <- function(gds_file, null_model, sample_id = NULL,
                          use_SPA = NULL, chunk_size = 2000L,
                          mac_cutoff = 1L, missing_imputation = "mean",
                          output_csv = NULL, verbose = 1) {
  if (!requireNamespace("SeqArray", quietly = TRUE)) {
    stop("SeqArray package required. Install with: BiocManager::install('SeqArray')")
  }

  # Validate inputs
  stopifnot(inherits(null_model, "glow_null_model"))
  stopifnot(file.exists(gds_file))
  stopifnot(is.numeric(chunk_size), length(chunk_size) == 1, chunk_size >= 1)

  missing_imputation <- match.arg(missing_imputation, c("mean", "zero"))

  trait <- null_model$trait

  # Resolve SPA setting
  if (is.null(use_SPA)) {
    use_SPA <- if (trait == "binary") TRUE else FALSE
  }
  if (identical(use_SPA, "both")) {
    use_SPA <- TRUE
  }
  if (isTRUE(use_SPA) && trait != "binary") {
    stop("SPA is only available for binary traits")
  }

  # Open GDS
  gds <- SeqArray::seqOpen(gds_file, readonly = TRUE)
  on.exit(SeqArray::seqClose(gds), add = TRUE)

  # Resolve sample IDs
  if (is.null(sample_id)) {
    sample_id <- null_model$sample_id
  }
  if (is.null(sample_id)) {
    stop("sample_id must be provided via null_model$sample_id or the sample_id parameter")
  }

  # Sample alignment
  gds_sample_ids <- SeqArray::seqGetData(gds, "sample.id")
  match_idx <- match(sample_id, gds_sample_ids)
  if (any(is.na(match_idx))) {
    n_missing <- sum(is.na(match_idx))
    stop(sprintf("%d sample IDs from null_model not found in GDS file", n_missing))
  }

  # Compute reorder index: GDS returns samples in their file order after
  # seqSetFilter. We need to map that back to null_model order.
  gds_order <- order(match_idx)
  reorder_idx <- order(gds_order)

  # Get all variant metadata (before any filter)
  all_variant_ids <- SeqArray::seqGetData(gds, "variant.id")
  n_variants <- length(all_variant_ids)

  if (verbose >= 1) {
    message(sprintf("  Scanning %d variants across %d samples (chunk_size=%d)",
                    n_variants, length(sample_id), chunk_size))
  }

  # Precompute scan components from null model
  components <- .precompute_scan_components(null_model)

  # Set up chunks
  n_chunks <- ceiling(n_variants / chunk_size)
  results_list <- vector("list", n_chunks)

  for (i in seq_len(n_chunks)) {
    start_idx <- (i - 1L) * chunk_size + 1L
    end_idx <- min(i * chunk_size, n_variants)
    chunk_variant_ids <- all_variant_ids[start_idx:end_idx]

    # Set filter for this chunk
    SeqArray::seqSetFilter(gds, sample.id = sample_id,
                           variant.id = chunk_variant_ids, verbose = FALSE)

    # Read variant metadata for this chunk
    chr <- SeqArray::seqGetData(gds, "chromosome")
    pos <- SeqArray::seqGetData(gds, "position")
    allele_str <- SeqArray::seqGetData(gds, "allele")
    vid <- SeqArray::seqGetData(gds, "variant.id")

    # Parse alleles
    allele_split <- strsplit(allele_str, ",")
    ref <- vapply(allele_split, `[`, character(1), 1)
    alt <- vapply(allele_split, `[`, character(1), 2)

    # Read genotype dosage (n_filtered_samples x n_chunk_variants)
    G <- SeqArray::seqGetData(gds, "$dosage")

    # Reorder rows to match null model sample order
    G <- G[reorder_idx, , drop = FALSE]

    # Process chunk
    if (isTRUE(use_SPA)) {
      chunk_result <- .process_chunk_SPA(G, components, null_model,
                                          missing_imputation)
    } else {
      chunk_result <- .process_chunk_standard(G, components,
                                               missing_imputation)
    }

    # Combine with variant metadata
    chunk_df <- data.frame(
      chr = chr,
      pos = pos,
      ref = ref,
      alt = alt,
      variant_id = vid,
      chunk_result,
      stringsAsFactors = FALSE
    )

    results_list[[i]] <- chunk_df

    if (verbose >= 1 && n_chunks > 1) {
      message(sprintf("  Chunk %d/%d (%d variants)", i, n_chunks,
                      end_idx - start_idx + 1L))
    }
  }

  # Combine all chunks
  results <- do.call(rbind, results_list)
  rownames(results) <- NULL

  # Apply MAC cutoff
  if (mac_cutoff > 0) {
    keep <- results$MAC >= mac_cutoff
    n_filtered <- sum(!keep)
    results <- results[keep, , drop = FALSE]
    if (verbose >= 1 && n_filtered > 0) {
      message(sprintf("  Filtered %d variants with MAC < %d", n_filtered,
                      mac_cutoff))
    }
  }

  # Set class
  class(results) <- c("glow_marginal_scan", "data.frame")

  # Write CSV if requested
  if (!is.null(output_csv)) {
    utils::write.csv(results, output_csv, row.names = FALSE)
    if (verbose >= 1) {
      message(sprintf("  Results written to %s", output_csv))
    }
  }

  results
}


#' Extract Phenotype and Covariates from GDS File
#'
#' Reads sample annotation fields from a GDS file and prepares them
#' for null model fitting. Handles PLINK phenotype encoding and
#' covariate conversions.
#'
#' @param gds_file Character path to GDS file.
#' @param pheno_name Character name of the sample.annotation field
#'   containing the phenotype. Default: "phenotype".
#' @param covar_names Character vector of sample.annotation field names
#'   to use as covariates. Default: NULL (no covariates beyond intercept).
#' @param pheno_coding Character encoding scheme for phenotype values.
#'   "plink" (1=control, 2=case, -9=missing) or "01" (already 0/1).
#'   Default: "plink".
#' @param verbose Integer verbosity level (0=silent, 1=progress, 2=detailed).
#'
#' @return A list with components:
#'   \describe{
#'     \item{sample_id}{Character vector of sample IDs (missing-phenotype
#'       samples excluded)}
#'     \item{Y}{Numeric response vector (0/1 for binary)}
#'     \item{X}{Numeric covariate matrix (columns named by covar_names),
#'       or NULL if no covariates}
#'     \item{trait}{Character: "binary" or "continuous" (auto-detected)}
#'     \item{n_total}{Integer: total samples in GDS}
#'     \item{n_excluded}{Integer: samples excluded due to missing phenotype}
#'   }
#'
#' @details
#' Opens the GDS file read-only using SeqArray and reads sample annotation
#' nodes. Phenotype encoding is applied before missing-value filtering, so
#' PLINK sentinel values (-9) are converted to NA and those samples are
#' excluded. Character covariates are converted: sex (M/F) becomes 1/0;
#' other character covariates become integer factor codes minus 1.
#'
#' Trait type is auto-detected after encoding: if all non-missing values are
#' in \{0, 1\} the trait is "binary", otherwise "continuous".
#'
#' Complexity: O(n) in sample count. The GDS file is opened and closed
#' within this call.
#'
#' @examples
#' \dontrun{
#' pheno <- extract_pheno_covar_gds(
#'   "chr22.gds",
#'   pheno_name   = "phenotype",
#'   covar_names  = c("sex", "age"),
#'   pheno_coding = "plink",
#'   verbose      = 1
#' )
#' head(pheno$Y)
#' dim(pheno$X)
#' }
#'
#' @seealso \code{\link{fit_null_model}}, \code{\link{marginal_scan}}
#'
#' @export
extract_pheno_covar_gds <- function(gds_file, pheno_name = "phenotype",
                                    covar_names = NULL, pheno_coding = "plink",
                                    verbose = 1) {
  if (!requireNamespace("SeqArray", quietly = TRUE)) {
    stop("SeqArray package required. Install with: BiocManager::install('SeqArray')")
  }
  if (!requireNamespace("gdsfmt", quietly = TRUE)) {
    stop("gdsfmt package required. Install with: BiocManager::install('gdsfmt')")
  }

  stopifnot(file.exists(gds_file))
  stopifnot(pheno_coding %in% c("plink", "01"))

  # Open GDS
  gds <- SeqArray::seqOpen(gds_file, readonly = TRUE)
  on.exit(SeqArray::seqClose(gds), add = TRUE)

  # Read sample IDs
  sample_id <- SeqArray::seqGetData(gds, "sample.id")
  n_total <- length(sample_id)

  # Read phenotype
  pheno_path <- paste0("sample.annotation/", pheno_name)
  pheno_node <- tryCatch(
    gdsfmt::index.gdsn(gds, pheno_path),
    error = function(e) {
      stop("Phenotype field '", pheno_name, "' not found in sample.annotation/")
    }
  )
  Y_raw <- gdsfmt::read.gdsn(pheno_node)

  # Convert phenotype encoding
  if (pheno_coding == "plink") {
    Y <- ifelse(Y_raw == 1, 0, ifelse(Y_raw == 2, 1, NA))
  } else {
    Y <- as.numeric(Y_raw)
    Y[Y < 0] <- NA
  }

  # Identify valid samples (non-NA phenotype)
  valid <- !is.na(Y)
  n_excluded <- sum(!valid)

  if (all(!valid)) {
    stop("All phenotype values are missing after encoding conversion")
  }

  if (verbose >= 1) {
    message(sprintf("  Phenotype '%s': %d valid, %d excluded (of %d total)",
                    pheno_name, sum(valid), n_excluded, n_total))
  }

  # Read covariates (if specified)
  X <- NULL
  if (!is.null(covar_names) && length(covar_names) > 0) {
    X_list <- list()
    for (name in covar_names) {
      covar_path <- paste0("sample.annotation/", name)
      node <- tryCatch(
        gdsfmt::index.gdsn(gds, covar_path),
        error = function(e) {
          stop("Covariate field '", name, "' not found in sample.annotation/")
        }
      )
      vals <- gdsfmt::read.gdsn(node)

      # Convert character covariates to numeric; report mapping at verbose >= 2
      if (is.character(vals)) {
        unique_vals <- unique(stats::na.omit(vals))
        if (all(unique_vals %in% c("M", "F"))) {
          if (verbose >= 2) {
            message(sprintf("  Covariate '%s': M -> 1, F -> 0", name))
          }
          vals <- as.numeric(vals == "M")
        } else {
          lev <- sort(unique_vals)
          coded <- as.integer(factor(vals, levels = lev)) - 1L
          if (verbose >= 2) {
            mapping_str <- paste(
              sprintf("'%s'->%d", lev, seq_along(lev) - 1L),
              collapse = ", "
            )
            message(sprintf("  Covariate '%s': %s", name, mapping_str))
          }
          vals <- coded
        }
      }
      X_list[[name]] <- vals
    }
    X <- do.call(cbind, X_list)
    colnames(X) <- covar_names
  }

  # Subset to valid samples
  sample_id <- sample_id[valid]
  Y <- Y[valid]
  if (!is.null(X)) {
    X <- X[valid, , drop = FALSE]
  }

  # Auto-detect trait type
  unique_Y <- unique(Y)
  trait <- if (all(unique_Y %in% c(0, 1))) "binary" else "continuous"

  if (verbose >= 1) {
    message(sprintf("  Trait type: %s", trait))
    if (trait == "binary") {
      message(sprintf("  Cases: %d, Controls: %d", sum(Y == 1), sum(Y == 0)))
    }
  }

  list(
    sample_id = sample_id,
    Y = Y,
    X = X,
    trait = trait,
    n_total = n_total,
    n_excluded = n_excluded
  )
}


#' Write Marginal Analysis Results to GDS Annotation Node
#'
#' Stores marginal analysis results as an annotation node in a GDS file
#' for downstream reuse. Follows the same pattern as FunctionalAnnotation
#' storage.
#'
#' @param gds_path Character path to GDS file to annotate.
#' @param results A data.frame (class "glow_marginal_scan" or plain
#'   data.frame) from \code{marginal_scan()}.
#' @param node_name Character name for the annotation node.
#'   Default: "MarginalAnalysis".
#' @param overwrite Logical. If TRUE, replace existing node. Default: FALSE.
#' @param verbose Integer verbosity level (0=silent, 1=progress).
#'
#' @return Invisible gds_path.
#'
#' @details
#' Opens the GDS file for writing and navigates to (or creates)
#' annotation/info/<node_name>. The results data.frame is aligned to the
#' GDS variant order via \code{match(gds_variant_ids, results$variant_id)};
#' any GDS variant not present in results will have NA values in the stored
#' matrix. A verbose message is emitted when unmatched variants are detected.
#'
#' Columns stored: MAF, MAC, score, var_score, Z, pvalue. When SPA results
#' are present in \code{results}, Z_SPA and pvalue_SPA are also stored.
#' The matrix is compressed with LZMA_RA. Metadata attributes
#' "column_names" and "analysis_date" are attached to the node.
#'
#' Complexity: O(V) in the total number of GDS variants. The file is
#' modified in-place; no temporary copy is created.
#'
#' @examples
#' \dontrun{
#' pheno   <- extract_pheno_covar_gds("chr22.gds", covar_names = "sex")
#' nm      <- fit_null_model(pheno$X, pheno$Y, trait = "binary",
#'                           sample_id = pheno$sample_id)
#' results <- marginal_scan("chr22.gds", nm, use_SPA = FALSE)
#' annotate_gds_marginal("chr22.gds", results, node_name = "MarginalAnalysis")
#' }
#'
#' @seealso \code{\link{marginal_scan}}, \code{\link{fit_null_model}}
#'
#' @export
annotate_gds_marginal <- function(gds_path, results,
                                   node_name = "MarginalAnalysis",
                                   overwrite = FALSE, verbose = 1) {
  if (!requireNamespace("gdsfmt", quietly = TRUE)) {
    stop("gdsfmt package required. Install with: BiocManager::install('gdsfmt')")
  }

  stopifnot(file.exists(gds_path))
  stopifnot(is.data.frame(results))
  stopifnot("variant_id" %in% names(results))

  # Select numeric columns for storage
  result_cols <- c("MAF", "MAC", "score", "var_score", "Z", "pvalue")
  if ("pvalue_SPA" %in% names(results)) {
    result_cols <- c(result_cols, "Z_SPA", "pvalue_SPA")
  }
  missing_cols <- setdiff(result_cols, names(results))
  if (length(missing_cols) > 0) {
    stop("Missing expected columns in results: ",
         paste(missing_cols, collapse = ", "))
  }

  # Open GDS for writing
  gds <- gdsfmt::openfn.gds(gds_path, readonly = FALSE)
  on.exit(gdsfmt::closefn.gds(gds), add = TRUE)

  # Navigate to annotation/info (create if needed)
  annotation_exists <- "annotation" %in% gdsfmt::ls.gdsn(gds)
  if (!annotation_exists) {
    annot_folder <- gdsfmt::addfolder.gdsn(gds, "annotation")
  } else {
    annot_folder <- gdsfmt::index.gdsn(gds, "annotation")
  }

  info_exists <- "info" %in% gdsfmt::ls.gdsn(annot_folder)
  if (!info_exists) {
    info_folder <- gdsfmt::addfolder.gdsn(annot_folder, "info")
  } else {
    info_folder <- gdsfmt::index.gdsn(annot_folder, "info")
  }

  # Check for existing node
  if (node_name %in% gdsfmt::ls.gdsn(info_folder)) {
    if (!overwrite) {
      stop("Node '", node_name, "' already exists. Use overwrite = TRUE.")
    }
    gdsfmt::delete.gdsn(gdsfmt::index.gdsn(info_folder, node_name))
  }

  # Match variant ordering: GDS variant order -> results order
  gds_variant_ids <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "variant.id"))
  match_idx <- match(gds_variant_ids, results$variant_id)

  # Report unmatched GDS variants (will be NA in the stored matrix)
  n_unmatched <- sum(is.na(match_idx))
  if (verbose >= 1 && n_unmatched > 0) {
    message(sprintf("  %d GDS variants not in results (filled with NA)",
                    n_unmatched))
  }

  # Build result matrix aligned to GDS variant order
  result_matrix <- as.matrix(results[match_idx, result_cols])
  rownames(result_matrix) <- NULL

  # Write matrix
  fa_node <- gdsfmt::add.gdsn(info_folder, node_name, result_matrix,
                                compress = "LZMA_RA", closezip = TRUE)
  gdsfmt::put.attr.gdsn(fa_node, "column_names", result_cols)
  gdsfmt::put.attr.gdsn(fa_node, "analysis_date", as.character(Sys.Date()))

  if (verbose >= 1) {
    message(sprintf("  Added %s node: %d variants x %d columns",
                    node_name, nrow(result_matrix), length(result_cols)))
  }

  invisible(gds_path)
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' Precompute Scan Components from Null Model
#'
#' @param null_model A glow_null_model object
#' @return List with H_half, Y_resid, weights, trait, s0
#' @keywords internal
#' @noRd
.precompute_scan_components <- function(null_model) {
  # Get or compute H_half
  H_half <- null_model$H_half
  if (is.null(H_half)) {
    # Backward compat: compute H_half from stored X
    X <- null_model$X
    if (null_model$trait == "binary") {
      w <- sqrt(null_model$fitted_probs * (1 - null_model$fitted_probs))
      X_tilde <- X * w
      XtX_chol <- chol(crossprod(X_tilde))
      H_half <- X_tilde %*% backsolve(XtX_chol, diag(ncol(XtX_chol)))
    } else {
      XtX_chol <- chol(crossprod(X))
      H_half <- X %*% backsolve(XtX_chol, diag(ncol(XtX_chol)))
    }
  }

  list(
    H_half = H_half,
    Y_resid = if (null_model$trait == "binary") {
      null_model$Y - null_model$fitted_probs
    } else {
      null_model$residuals / null_model$s0
    },
    weights = null_model$weights,
    trait = null_model$trait,
    s0 = null_model$s0
  )
}


#' Process One Genotype Chunk Using Standard Asymptotic Method
#'
#' @param G Genotype dosage matrix (n x m)
#' @param components List from .precompute_scan_components()
#' @param missing_imputation Character: "mean" or "zero"
#' @return data.frame with MAF, MAC, score, var_score, Z, pvalue
#' @keywords internal
#' @noRd
.process_chunk_standard <- function(G, components, missing_imputation = "mean") {
  geno_info <- .impute_geno(G, method = missing_imputation)
  G <- geno_info$G

  # Compute scores: S = G'(Y - mu_0) for binary, G'r/s0 for continuous
  scores <- as.numeric(crossprod(G, components$Y_resid))

  # Compute diagonal of score variance
  if (components$trait == "binary") {
    Gw <- G * components$weights
    diag_var <- colSums(Gw^2) - colSums(crossprod(components$H_half, Gw)^2)
    var_scale <- colSums(Gw^2)
  } else {
    diag_var <- colSums(G^2) - colSums(crossprod(components$H_half, G)^2)
    var_scale <- colSums(G^2)
  }

  # Handle zero-variance (monomorphic / degenerate variants). `diag_var` is a
  # difference of two positive sums and suffers catastrophic cancellation for
  # near-constant genotypes, so an absolute `.Machine$double.eps` test classifies
  # the resulting float noise INCONSISTENTLY across `chunk_size` (colSums/crossprod
  # round differently for a 1- vs many-column matrix, flipping a variant in/out of
  # "zero variance"). Detect it deterministically instead: MAC == 0 is monomorphic
  # by definition (chunk-invariant), plus a scale-relative guard for near-degenerate
  # variance. This makes Z (= 0 here) chunk-invariant. See test "chunk_size = 1".
  zero_var <- geno_info$MAC == 0 | !is.finite(diag_var) |
    diag_var < sqrt(.Machine$double.eps) * pmax(var_scale, 1)
  zero_var[is.na(zero_var)] <- TRUE   # all-missing / undefined MAC -> degenerate
  diag_var[zero_var] <- 1

  Z <- scores / sqrt(diag_var)
  Z[zero_var] <- 0
  pvalue <- 2 * stats::pnorm(-abs(Z))
  pvalue[zero_var] <- 1

  data.frame(MAF = geno_info$MAF, MAC = geno_info$MAC, score = scores,
             var_score = diag_var, Z = Z, pvalue = pvalue)
}


#' Process One Genotype Chunk Using SPA
#'
#' @param G Genotype dosage matrix (n x m)
#' @param components List from .precompute_scan_components()
#' @param null_model glow_null_model for SPAtest inputs
#' @param missing_imputation Character: "mean" or "zero"
#' @return data.frame with MAF, MAC, score, var_score, Z, pvalue, Z_SPA,
#'   pvalue_SPA
#' @keywords internal
#' @noRd
.process_chunk_SPA <- function(G, components, null_model,
                                missing_imputation = "mean") {
  geno_info <- .impute_geno(G, method = missing_imputation)
  G <- geno_info$G

  # Standard Z-scores (always computed)
  scores <- as.numeric(crossprod(G, components$Y_resid))

  Gw <- G * components$weights
  diag_var <- colSums(Gw^2) - colSums(crossprod(components$H_half, Gw)^2)

  # Zero-variance detection: deterministic + chunk-invariant (catastrophic
  # cancellation makes an absolute-eps test chunk-dependent; see
  # .process_chunk_standard). MAC == 0 is monomorphic; the relative guard covers
  # near-degenerate variance.
  zero_var <- geno_info$MAC == 0 | !is.finite(diag_var) |
    diag_var < sqrt(.Machine$double.eps) * pmax(colSums(Gw^2), 1)
  zero_var[is.na(zero_var)] <- TRUE   # all-missing / undefined MAC -> degenerate
  diag_var_safe <- diag_var
  diag_var_safe[zero_var] <- 1

  Z <- scores / sqrt(diag_var_safe)
  Z[zero_var] <- 0
  pvalue <- 2 * stats::pnorm(-abs(Z))
  pvalue[zero_var] <- 1

  # SPA p-values via SPAtest
  # SPAtest expects: genos = SNPs in rows, samples in columns
  # cov should NOT include intercept (SPAtest adds its own)
  # Verify first column of X is the intercept before dropping it
  stopifnot(all(null_model$X[, 1] == 1))
  X_no_intercept <- null_model$X[, -1, drop = FALSE]

  spa_result <- SPAtest::ScoreTest_SPA(
    genos = t(G),
    pheno = null_model$Y,
    cov = X_no_intercept,
    method = "fastSPA",
    minmac = 1L,
    Cutoff = "BE"
  )

  pvalue_SPA <- spa_result$p.value
  Z_SPA <- stats::qnorm(pvalue_SPA / 2, lower.tail = FALSE) * sign(scores)
  # Handle NA from SPAtest (e.g., for monomorphic variants)
  Z_SPA[is.na(pvalue_SPA)] <- NA

  data.frame(MAF = geno_info$MAF, MAC = geno_info$MAC, score = scores,
             var_score = diag_var, Z = Z, pvalue = pvalue,
             Z_SPA = Z_SPA, pvalue_SPA = pvalue_SPA)
}


#' Impute Missing Genotypes and Compute Allele Frequencies
#'
#' @param G Genotype dosage matrix (n x m), may contain NAs
#' @param method Character: "mean" or "zero"
#' @return List with G (imputed), MAF, MAC
#' @keywords internal
#' @noRd
.impute_geno <- function(G, method = "mean") {
  m <- ncol(G)
  n <- nrow(G)

  # Compute allele frequency before imputation
  AF <- colMeans(G, na.rm = TRUE) / 2
  MAF <- pmin(AF, 1 - AF)

  # Impute missing
  if (method == "mean") {
    for (j in seq_len(m)) {
      na_idx <- is.na(G[, j])
      if (any(na_idx)) {
        G[na_idx, j] <- 2 * AF[j]
      }
    }
  } else if (method == "zero") {
    G[is.na(G)] <- 0
  }

  # Compute MAC (minor allele count)
  MAC <- as.integer(round(colSums(G)))
  MAC <- pmin(MAC, 2L * as.integer(n) - MAC)

  list(G = G, MAF = MAF, MAC = MAC)
}
