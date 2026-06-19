# ==============================================================================
# Region Generator Functions for Variant-Set Analysis
# ==============================================================================
#
# Functions for defining genomic regions (gene-based, sliding window, custom)
# for use with extract_variant_set().
#
# EXPORTED FUNCTIONS:
#   - define_regions_gene()    Gene-based regions from gene boundary table
#   - define_regions_window()  Sliding window regions
#   - define_regions_custom()  User-provided custom regions
#
# INTERNAL HELPERS:
#   - .normalize_gene_table()  Normalize gene table column formats
#   - .empty_regions_df()      Empty regions data frame template

#################### EXPORTED MAIN FUNCTIONS ####################

#' Define Gene-Based Genomic Regions
#'
#' Creates a standard region table from a gene boundary table. Each gene
#' becomes one region for variant-set analysis.
#'
#' @param genes_info Data frame with gene boundaries. Default NULL uses
#'   the bundled STAARpipeline genes_info (18,445 protein-coding genes).
#'   Accepts STAARpipeline format (hgnc_symbol, chromosome_name, start_position,
#'   end_position) or simplified format (gene, chr, start, end).
#' @param chr Integer or character chromosome(s) to include. NULL for all.
#' @param extend Integer bp to extend gene boundaries on each side. Default 0.
#'
#' @return Data frame with columns: region_id, chr, start, end, label.
#'
#' @examples
#' # With a user-provided gene table
#' my_genes <- data.frame(
#'   gene = c("BRCA1", "TP53"),
#'   chr = c("17", "17"),
#'   start = c(43044295L, 7668402L),
#'   end = c(43125364L, 7687550L)
#' )
#' regions <- define_regions_gene(my_genes)
#' regions
#'
#' # Extend boundaries by 1kb on each side
#' regions_ext <- define_regions_gene(my_genes, extend = 1000)
#'
#' @export
define_regions_gene <- function(genes_info = NULL, chr = NULL, extend = 0L) {

  # 1. Use bundled data if NULL
  if (is.null(genes_info)) {
    genes_info <- get("genes_info", envir = asNamespace("GLOWr"))
  }

  # 2. Detect and normalize column names
  genes_info <- .normalize_gene_table(genes_info)

  # 3. Filter by chromosome
  if (!is.null(chr)) {
    chr <- as.character(chr)
    genes_info <- genes_info[as.character(genes_info$chr) %in% chr, , drop = FALSE]
  }

  if (nrow(genes_info) == 0) {
    warning("No genes found for the specified chromosome(s)")
    return(.empty_regions_df())
  }

  # 4. Apply extension
  if (extend > 0) {
    genes_info$start <- pmax(1L, genes_info$start - as.integer(extend))
    genes_info$end <- genes_info$end + as.integer(extend)
  }

  # 5. Format output
  data.frame(
    region_id = genes_info$gene,
    chr = genes_info$chr,
    start = genes_info$start,
    end = genes_info$end,
    label = genes_info$gene,
    stringsAsFactors = FALSE
  )
}

#' Define Sliding Window Regions
#'
#' Generates overlapping fixed-size genomic windows for variant-set analysis.
#'
#' @param chr Chromosome (integer or character).
#' @param start Integer start position. Required.
#' @param end Integer end position. Required.
#' @param window_size Integer window width in bp. Default 2000.
#' @param step_size Integer step between windows. Default window_size/2 (50% overlap).
#'
#' @return Data frame with columns: region_id, chr, start, end, label.
#'
#' @examples
#' # 2kb windows with 50% overlap (default)
#' wins <- define_regions_window(chr = 22, start = 16570000, end = 16580000)
#' wins
#'
#' # 5kb non-overlapping windows
#' wins2 <- define_regions_window(chr = 1, start = 1000000, end = 1050000,
#'                                window_size = 5000, step_size = 5000)
#'
#' @export
define_regions_window <- function(chr, start, end,
                                   window_size = 2000L, step_size = NULL) {
  chr <- as.character(chr)
  start <- as.integer(start)
  end <- as.integer(end)
  window_size <- as.integer(window_size)
  if (is.null(step_size)) step_size <- as.integer(window_size / 2)
  step_size <- as.integer(step_size)

  stopifnot(end > start, window_size > 0, step_size > 0)

  # Generate window starts
  n_windows <- max(1L, floor((end - start + 1L - window_size) / step_size) + 1L)
  win_starts <- start + (seq_len(n_windows) - 1L) * step_size
  win_ends <- pmin(win_starts + window_size - 1L, end)

  data.frame(
    region_id = paste0("win_", seq_len(n_windows)),
    chr = rep(chr, n_windows),
    start = win_starts,
    end = win_ends,
    label = paste0("chr", chr, ":", win_starts, "-", win_ends),
    stringsAsFactors = FALSE
  )
}

#' Define Custom Genomic Regions
#'
#' Validates and standardizes a user-provided region table.
#'
#' @param regions_table Data frame with at least: region_id, chr, start, end.
#'   Optional: label (defaults to region_id).
#'
#' @return Standardized data frame with columns: region_id, chr, start, end, label.
#'
#' @examples
#' my_regions <- data.frame(
#'   region_id = c("enhancer_1", "enhancer_2"),
#'   chr = c("1", "1"),
#'   start = c(1000000L, 2000000L),
#'   end = c(1005000L, 2005000L)
#' )
#' regions <- define_regions_custom(my_regions)
#' regions
#'
#' @export
define_regions_custom <- function(regions_table) {
  required <- c("region_id", "chr", "start", "end")
  missing_cols <- setdiff(required, names(regions_table))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  out <- data.frame(
    region_id = as.character(regions_table$region_id),
    chr = as.character(regions_table$chr),
    start = as.integer(regions_table$start),
    end = as.integer(regions_table$end),
    stringsAsFactors = FALSE
  )
  out$label <- if ("label" %in% names(regions_table)) {
    as.character(regions_table$label)
  } else {
    out$region_id
  }
  out
}

#################### INTERNAL HELPER FUNCTIONS ####################

#' Normalize gene table column formats
#'
#' Detects the format of a gene boundary table and normalizes it to a
#' standard format with columns: gene, chr, start, end.
#'
#' @param df Data frame with gene boundaries in any supported format.
#'
#' @return Data frame with columns: gene, chr, start, end.
#'
#' @keywords internal
#' @noRd
.normalize_gene_table <- function(df) {
  # STAARpipeline format
  if (all(c("hgnc_symbol", "chromosome_name", "start_position", "end_position") %in% names(df))) {
    return(data.frame(
      gene = df$hgnc_symbol,
      chr = as.character(df$chromosome_name),
      start = as.integer(df$start_position),
      end = as.integer(df$end_position),
      stringsAsFactors = FALSE
    ))
  }
  # Simple gene-list format: 4-column df where 4th col is character (gene names)
  if (ncol(df) == 4 && is.character(df[[4]])) {
    return(data.frame(
      gene = df[[4]],
      chr = as.character(df[[1]]),
      start = as.integer(df[[2]]),
      end = as.integer(df[[3]]),
      stringsAsFactors = FALSE
    ))
  }
  # Named format: gene, chr, start, end
  if (all(c("gene", "chr", "start", "end") %in% names(df))) {
    return(data.frame(
      gene = df$gene,
      chr = as.character(df$chr),
      start = as.integer(df$start),
      end = as.integer(df$end),
      stringsAsFactors = FALSE
    ))
  }
  stop("Unrecognized gene table format. Expected columns: ",
       "(hgnc_symbol, chromosome_name, start_position, end_position) or ",
       "(gene, chr, start, end)")
}

#' Empty regions data frame template
#'
#' Returns an empty data frame with the standard region table columns.
#'
#' @return Empty data frame with columns: region_id, chr, start, end, label.
#'
#' @keywords internal
#' @noRd
.empty_regions_df <- function() {
  data.frame(
    region_id = character(0),
    chr = character(0),
    start = integer(0),
    end = integer(0),
    label = character(0),
    stringsAsFactors = FALSE
  )
}
