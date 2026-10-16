# ==============================================================================
# File: noncoding_regions.R
# Purpose: the grouping pass of the gene-centric non-coding masks (design
#   define_regions_noncoding() runs once per chromosome and mask: it reads the
#   chromosome's QC label, variant type, GENCODE category and the presence of
#   the CAGE, DHS and GeneHancer annotations once (the context, shared by the
#   masks of one unit), selects the mask's records, intersects a promoter mask
#   with the TSS windows, reads the linkage field for the selected records
#   only, links them to genes (noncoding_linkage.R), restricts the genes to the
#   mask's gene table, and returns one region per gene whose candidates are the
#   linked records' row numbers in the chromosome's variant index, unique and
#   increasing, as a list column that extract_variant_set() takes in place of
#   its interval lookup.
#
# EXPORTED FUNCTIONS:
#   - define_regions_noncoding()   One mask's regions on one chromosome
#
# INTERNAL HELPERS:
#   - .noncoding_chr_context()     The chromosome-wide vectors, read once
#   - .noncoding_gene_list()       The mask's gene table on one chromosome
#   - .promoter_windows_hg38()     The TSS +/- 3 kb windows (TxDb, memoised)
#   - .in_promoter_windows()       Positions inside the windows of a chromosome
#   - .empty_noncoding_regions()   The 0-row regions table

# Session cache of the promoter windows (and a hook for tests, which place a
# synthetic GRanges here).
.noncoding_cache <- new.env(parent = emptyenv())


#################### EXPORTED MAIN FUNCTIONS ####################

#' Define Gene-Centric Non-Coding Regions for One Chromosome and Mask
#'
#' Builds the regions of one gene-centric non-coding mask on one chromosome:
#' one region per gene, whose candidate variants are the records that the
#' mask's annotation links to that gene, wherever they lie on the chromosome.
#' This is how STAARpipeline's gene-centric non-coding analysis forms its
#' sets, computed here once per chromosome and mask instead of once per gene.
#' The regions table carries the candidates as a list column
#' \code{variant_rows} that \code{\link{extract_variant_set}} takes in place
#' of its interval lookup; everything downstream (the filter's QC, type, MAF,
#' MAC and minimum-size tests, LD pruning, collapsing, the tests) is unchanged.
#'
#' @section The eight masks:
#' Every mask starts from the records that pass the filter's QC label and
#' variant type. The masks read five annotation fields of the aGDS, named by
#' their catalog names: \code{GENCODE.Category}, \code{GENCODE.Info},
#' \code{GeneHancer}, \code{CAGE} and \code{DHS}; \code{\link{annotate_favor}}
#' writes all five by default.
#' \tabular{lllll}{
#'   \strong{Mask} \tab \strong{Gene table} \tab \strong{Selection} \tab \strong{Linkage field} \tab \strong{TSS window} \cr
#'   \code{upstream} \tab protein-coding \tab category has \code{upstream} \tab \code{GENCODE.Info} \tab no \cr
#'   \code{downstream} \tab protein-coding \tab category has \code{downstream} \tab \code{GENCODE.Info} \tab no \cr
#'   \code{UTR} \tab protein-coding \tab category has \code{UTR3} or \code{UTR5} \tab \code{GENCODE.Info} \tab no \cr
#'   \code{promoter_CAGE} \tab protein-coding \tab \code{CAGE} non-empty \tab \code{GENCODE.Info} \tab yes \cr
#'   \code{promoter_DHS} \tab protein-coding \tab \code{DHS} non-empty \tab \code{GENCODE.Info} \tab yes \cr
#'   \code{enhancer_CAGE} \tab protein-coding \tab \code{CAGE} and \code{GeneHancer} non-empty \tab \code{GeneHancer} \tab no \cr
#'   \code{enhancer_DHS} \tab protein-coding \tab \code{DHS} and \code{GeneHancer} non-empty \tab \code{GeneHancer} \tab no \cr
#'   \code{ncRNA} \tab ncRNA genes \tab category has \code{ncRNA_exonic} or \code{ncRNA_splicing} \tab \code{GENCODE.Info} \tab no
#' }
#' "Has" means that the ";"-separated token list of the category contains the
#' value. The protein-coding table is the bundled \code{genes_info}, the ncRNA
#' table the bundled \code{ncRNA_gene} (21,104 names). The two promoter masks
#' keep, among the selected records, those inside the TSS \eqn{\pm} 3 kb
#' windows of \code{GenomicFeatures::promoters(genes(TxDb.Hsapiens.UCSC.hg38.knownGene),
#' 3000, 3000)}, genome-wide, as STAARpipeline does: the window is not tied to
#' the linked gene, so a variant can lie in gene A's window and be linked to
#' gene B by its annotation. The packages \code{TxDb.Hsapiens.UCSC.hg38.knownGene},
#' \code{GenomicFeatures}, \code{GenomicRanges} and \code{IRanges} are needed
#' for those two masks only; the aGDS must be on hg38, which the code does not
#' verify.
#'
#' @section The linkage rules:
#' Under the default \code{rule = "glow"} a record is linked by the meaning of
#' its annotation field, so that FAVOR v1 strings and the FAVOR 2.0 layout
#' written by \code{\link{annotate_favor}} give the same links wherever the
#' information is present in both. Every balanced parenthesized detail is
#' removed first; ";" separates the gene groups of a multi-category record
#' (one group per category token, in order) and "," the genes of a group;
#' \code{NONE} is a sentinel; hyphenated names stay whole.
#' \itemize{
#'   \item \code{upstream}, \code{downstream}: the genes of the group at the
#'     position of the mask's category token, so an \code{upstream;downstream}
#'     record joins both masks with its own group in each.
#'   \item \code{UTR}: all genes of the UTR groups (a record annotated to two
#'     genes' UTRs joins both).
#'   \item \code{promoter_CAGE}, \code{promoter_DHS}: the first gene of the
#'     first group; a \code{NONE} first token links to no gene.
#'   \item \code{enhancer_CAGE}, \code{enhancer_DHS}: the first listed
#'     connected gene of the first GeneHancer element; a header without
#'     targets links to no gene. Later targets are not used.
#'   \item \code{ncRNA}: all genes of all groups, restricted to the ncRNA table.
#' }
#' In the FAVOR 2.0 layout the group boundaries of a two-category record are
#' not serialized, and such a record is then linked to every gene it names
#' under each of its masks. The link is then restricted to the mask's gene
#' table, without fallback to another gene. \code{rule = "staar"} reproduces
#' STAARpipeline 0.9.8's literal string operations mask by mask (the whole
#' category string must equal one of the values the function tests, so an
#' \code{upstream;downstream} record joins neither of those two masks, while
#' the UTR function accepts \code{UTR5;UTR3} and the ncRNA function
#' \code{ncRNA_exonic;splicing}; the UTR text before the first "("; the
#' promoter first token after a split on "(",
#' ",", ";" and "-"; the ncRNA first ";" group with one greedy removal of its
#' parenthesized text and its first three "," tokens) and exists for the
#' validation that enumerates the differences.
#'
#' @param gds An open \code{SeqVarGDSClass} of ONE chromosome, or a path to
#'   such a file. A file holding several chromosomes stops with a message.
#' @param chr The chromosome, as the GDS names it (for example \code{22}).
#' @param category One of the eight mask names.
#' @param filter_spec The run's \code{\link{variant_filter}}: its QC label,
#'   pass value and variant type are applied here. Its MAF, MAC and
#'   \code{min_variants} are NOT applied here; they belong to the per-set
#'   extraction, as for every other region type. Its annotation clauses are
#'   ignored.
#' @param variant_index A \code{glow_variant_index} of the GDS from
#'   \code{\link{build_variant_index}}, or NULL (built here).
#' @param context A context from an earlier call on the same file (the
#'   attribute \code{context} of a returned table), or NULL (built here). A
#'   scan builds one per chromosome and passes it to every mask, so the
#'   chromosome-wide vectors are read once.
#' @param rule \code{"glow"} (default) or \code{"staar"}.
#' @param gene_list NULL (the bundled table of the mask's family, on
#'   \code{chr}) or a character vector of gene names to restrict the links to.
#' @param Annotation_dir Character path prefix for annotations in the GDS.
#' @param Annotation_name_catalog Data frame mapping annotation names to GDS
#'   sub-paths. NULL uses the bundled catalog.
#' @param verbose Integer. 1 prints the selection counts.
#'
#' @return A data.frame with one row per gene with at least one linked record,
#'   in gene-table order: \code{region_id} and \code{label} (the gene),
#'   \code{chr}, \code{start} and \code{end} (the first and last position of
#'   the linked records, a descriptive span that is wide for an enhancer far
#'   from its gene), \code{n_records} (the linked records), \code{category}
#'   (the mask), \code{gene_rank} (the gene's row in its table) and the list
#'   column \code{variant_rows} (the linked records' row numbers in the
#'   chromosome's variant index, unique and increasing). Attributes:
#'   \code{context} (reusable for the next mask), \code{n_selected_records}
#'   (the records the mask selected before linkage) and
#'   \code{n_linked_records} (the (gene, record) links kept).
#'
#' @examples
#' \dontrun{
#' idx <- build_variant_index("chr22.gds")
#' spec <- variant_filter(rare_maf_cutoff = 0.01)
#' up <- define_regions_noncoding("chr22.gds", 22, "upstream", spec, idx)
#' ctx <- attr(up, "context")                  # reuse for the next mask
#' enh <- define_regions_noncoding("chr22.gds", 22, "enhancer_DHS", spec, idx,
#'                                 context = ctx)
#' vset <- extract_variant_set("chr22.gds", up[1, ], spec, variant_index = idx)
#' }
#'
#' @seealso \code{\link{noncoding_filter}}, \code{\link{extract_variant_set}},
#'   \code{\link{define_regions_gene}}
#' @export
define_regions_noncoding <- function(gds, chr, category, filter_spec,
                                     variant_index = NULL, context = NULL,
                                     rule = c("glow", "staar"), gene_list = NULL,
                                     Annotation_dir = "annotation/info/FunctionalAnnotation",
                                     Annotation_name_catalog = NULL, verbose = 0) {
  rule <- match.arg(rule)
  specs <- .noncoding_mask_specs()
  if (!is.character(category) || length(category) != 1L || !category %in% names(specs)) {
    stop("`category` must be one of: ", paste(names(specs), collapse = ", "), call. = FALSE)
  }
  spec <- specs[[category]]
  if (!inherits(filter_spec, "glow_variant_filter")) {
    stop("`filter_spec` must be a glow_variant_filter from variant_filter().", call. = FALSE)
  }
  if (!requireNamespace("SeqArray", quietly = TRUE)) {
    stop("SeqArray package required. Install with: BiocManager::install('SeqArray')")
  }
  if (is.character(gds)) {
    stopifnot(file.exists(gds))
    gds <- SeqArray::seqOpen(gds, readonly = TRUE, allow.duplicate = TRUE)
    on.exit(SeqArray::seqClose(gds), add = TRUE)
  }
  if (is.null(Annotation_name_catalog)) {
    Annotation_name_catalog <- get("annotation_name_catalog", envir = asNamespace("GLOWr"))
  }
  chr <- as.character(chr)

  # ---- 1. The chromosome context (read once, reused across the masks) ----
  if (is.null(context)) {
    context <- .noncoding_chr_context(gds, filter_spec, Annotation_dir, Annotation_name_catalog)
  } else if (!inherits(context, "glow_noncoding_context")) {
    stop("`context` must come from an earlier define_regions_noncoding() call ",
         "(its 'context' attribute), or be NULL.", call. = FALSE)
  } else if (!identical(context$qc_label, filter_spec$qc_label) ||
             !identical(context$qc_pass_value, filter_spec$qc_pass_value) ||
             !identical(context$variant_type, filter_spec$variant_type)) {
    # The cached pass vector was computed under the context's filter settings.
    stop("`context` was built under another filter (QC label '", context$qc_label,
         "', pass value '", context$qc_pass_value, "', variant type '",
         context$variant_type, "'); this call's filter has QC label '",
         filter_spec$qc_label, "', pass value '", filter_spec$qc_pass_value,
         "', variant type '", filter_spec$variant_type,
         "'. Build the context under the run's filter, or pass NULL.", call. = FALSE)
  }
  if (!identical(context$chr, chr)) {
    stop("The GDS holds chromosome ", context$chr, ", not the requested ", chr, ".",
         call. = FALSE)
  }
  if (is.null(variant_index)) {
    variant_index <- build_variant_index(gds)
  } else if (!inherits(variant_index, "glow_variant_index")) {
    stop("`variant_index` must be a glow_variant_index from build_variant_index(), or NULL.")
  }
  ix <- variant_index[[chr]]
  if (is.null(ix) || length(ix$pos) != context$n) {
    stop("`variant_index` does not describe this GDS (chromosome ", chr, ", ",
         context$n, " records).", call. = FALSE)
  }

  # ---- 2. The mask's gene table on this chromosome ----
  if (is.null(gene_list)) gene_list <- .noncoding_gene_list(spec$gene_list, chr)
  gene_list <- unique(as.character(gene_list))

  # ---- 3. Selection: QC and type, the category under the rule, the presence
  #         vectors, and for a promoter mask the TSS windows ----
  keep <- context$pass & .category_select(context$category, spec, rule)
  for (p in spec$presence) keep <- keep & context[[p]]
  if (isTRUE(spec$promoter)) {
    sel <- which(keep)                              # subset-assign, never recycle
    keep[sel] <- .in_promoter_windows(context$pos[sel], chr, .promoter_windows_hg38())
  }
  rows <- which(keep)
  if (verbose >= 1) {
    message(sprintf("%s chr%s: %d record(s) selected under rule '%s'",
                    category, chr, length(rows), rule))
  }
  if (length(rows) == 0L) {
    return(.empty_noncoding_regions(context, n_selected = 0L, n_linked = 0L))
  }

  # ---- 4. The linkage field, read for the selected records only ----
  SeqArray::seqSetFilter(gds, variant.sel = ix$sel[rows], verbose = FALSE)
  field <- .read_noncoding_field(gds, spec$linkage, Annotation_dir, Annotation_name_catalog)
  SeqArray::seqResetFilter(gds, verbose = FALSE)

  # ---- 5. The links, as a long (record row, gene) table ----
  genes <- if (identical(spec$linkage, "GeneHancer")) {
    g <- .link_genehancer(field, rule)
    lapply(g, function(x) if (is.na(x)) character(0) else x)
  } else {
    .link_gencode(field, context$category[rows], category, rule)
  }
  long <- data.frame(row = rep.int(rows, lengths(genes)),
                     gene = as.character(unlist(genes)), stringsAsFactors = FALSE)
  long <- long[long$gene %in% gene_list, , drop = FALSE]  # the mask's table, no fallback
  long <- unique(long)                                    # repeated (gene, record) pairs
  if (nrow(long) == 0L) {
    return(.empty_noncoding_regions(context, n_selected = length(rows), n_linked = 0L))
  }

  # ---- 6. One region per gene: unique increasing rows, the descriptive span ----
  by_gene   <- split(long$row, long$gene)
  rows_list <- lapply(by_gene, function(r) sort(unique(as.integer(r))))
  genes_out <- names(by_gene)
  pos <- context$pos
  out <- data.frame(
    region_id = genes_out,
    chr       = rep(chr, length(genes_out)),
    start     = vapply(rows_list, function(r) min(pos[r]), integer(1)),
    end       = vapply(rows_list, function(r) max(pos[r]), integer(1)),
    label     = genes_out,
    n_records = lengths(rows_list),
    category  = rep(category, length(genes_out)),
    gene_rank = match(genes_out, gene_list),
    stringsAsFactors = FALSE)
  out$variant_rows <- I(unname(rows_list))
  out <- out[order(out$gene_rank), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "context") <- context
  attr(out, "n_selected_records") <- length(rows)
  attr(out, "n_linked_records") <- nrow(long)
  out
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' The Chromosome-Wide Vectors the Masks Share
#'
#' Reads once, for the whole file: the chromosome (one value, else stop), the
#' positions, the QC and variant-type pass vector of the filter (a missing QC
#' value does not pass, as everywhere in GLOWr), the GENCODE category (NA as
#' empty) and the presence of the CAGE, DHS and GeneHancer annotations as
#' logical vectors (the strings are dropped after the test and the logical
#' vectors are consumed directly). All vectors are in file order, which for a
#' one-chromosome file is the order of the chromosome's variant index.
#'
#' @param gds An open SeqVarGDSClass.
#' @param filter_spec A glow_variant_filter.
#' @param Annotation_dir,Annotation_name_catalog As in define_regions_noncoding().
#' @return A \code{glow_noncoding_context}.
#' @keywords internal
#' @noRd
.noncoding_chr_context <- function(gds, filter_spec, Annotation_dir, Annotation_name_catalog) {
  SeqArray::seqResetFilter(gds, verbose = FALSE)
  chrs <- as.character(SeqArray::seqGetData(gds, "chromosome"))
  u <- unique(chrs)
  if (length(u) != 1L) {
    stop("define_regions_noncoding() needs a GDS of one chromosome; this file holds ",
         length(u), " (", paste(utils::head(u, 5L), collapse = ", "),
         if (length(u) > 5L) ", ..." else "", ").", call. = FALSE)
  }
  n <- length(chrs)
  pos <- as.integer(SeqArray::seqGetData(gds, "position"))
  qc <- tryCatch(SeqArray::seqGetData(gds, filter_spec$qc_label),
                 error = function(e) rep(filter_spec$qc_pass_value, n))
  pass <- qc == filter_spec$qc_pass_value
  pass[is.na(pass)] <- FALSE
  if (filter_spec$variant_type != "variant") {
    is_snv <- .is_snv_from_alleles(SeqArray::seqGetData(gds, "allele"))
    pass <- pass & (if (filter_spec$variant_type == "SNV") is_snv else !is_snv)
  }
  read_field <- function(name)
    .read_noncoding_field(gds, name, Annotation_dir, Annotation_name_catalog)
  category <- read_field("GENCODE.Category")
  category[is.na(category)] <- ""
  presence <- function(name) {
    v <- read_field(name)
    !is.na(v) & nzchar(v)
  }
  structure(list(chr = u, n = n, pos = pos, pass = pass, category = category,
                 cage = presence("CAGE"), dhs = presence("DHS"),
                 genehancer = presence("GeneHancer"),
                 qc_label = filter_spec$qc_label, qc_pass_value = filter_spec$qc_pass_value,
                 variant_type = filter_spec$variant_type),
            class = "glow_noncoding_context")
}


#' Read One Annotation Field of the Masks, Naming It When It Is Missing
#'
#' A missing node otherwise surfaces as SeqArray's path error, which does not
#' say which of the five fields the masks read is absent or that
#' \code{annotate_favor()} writes it.
#'
#' @param name The catalog name (\code{"GENCODE.Category"}, \code{"GENCODE.Info"},
#'   \code{"GeneHancer"}, \code{"CAGE"}, \code{"DHS"}).
#' @return A character vector over the GDS's current variant filter.
#' @keywords internal
#' @noRd
.read_noncoding_field <- function(gds, name, Annotation_dir, Annotation_name_catalog) {
  path <- .resolve_annotation_path(name, Annotation_dir, Annotation_name_catalog)
  v <- tryCatch(SeqArray::seqGetData(gds, path), error = function(e) {
    stop("The aGDS has no annotation node for the field '", name, "' (resolved path '",
         path, "'), one of the five fields the gene-centric non-coding masks read ",
         "(GENCODE.Category, GENCODE.Info, GeneHancer, CAGE, DHS). annotate_favor() ",
         "writes them by default; Annotation_dir and Annotation_name_catalog name ",
         "where they are. SeqArray said: ", conditionMessage(e), call. = FALSE)
  })
  as.character(v)
}


#' The Mask's Gene Table on One Chromosome
#'
#' @param kind \code{"protein_coding"} (the bundled genes_info, in its order)
#'   or \code{"ncrna"} (the bundled ncRNA_gene).
#' @param chr Character chromosome.
#' @return A character vector of gene names.
#' @keywords internal
#' @noRd
.noncoding_gene_list <- function(kind, chr) {
  if (identical(kind, "ncrna")) {
    nc <- get("ncRNA_gene", envir = asNamespace("GLOWr"))
    return(as.character(nc$ncRNA[as.character(nc$chr) == chr]))
  }
  gi <- get("genes_info", envir = asNamespace("GLOWr"))
  as.character(gi$hgnc_symbol[as.character(gi$chromosome_name) == chr])
}


#' The TSS +/- 3 kb Windows of hg38, Memoised for the Session
#'
#' \code{promoters(genes(TxDb.Hsapiens.UCSC.hg38.knownGene), upstream = 3000,
#' downstream = 3000)}, as STAARpipeline computes them: 30,733 windows, since
#' \code{genes()} drops the genes that cannot be represented as a single
#' range (its message is suppressed here, as is the out-of-bound warning on
#' three alternate contigs). The packages are suggested, not imported, and the
#' call stops with an install hint when one is missing.
#'
#' @return A GRanges of windows.
#' @keywords internal
#' @noRd
.promoter_windows_hg38 <- function() {
  w <- .noncoding_cache$promoter_windows
  if (!is.null(w)) return(w)
  for (pkg in c("GenomicFeatures", "TxDb.Hsapiens.UCSC.hg38.knownGene",
                "GenomicRanges", "IRanges")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("The promoter masks need the ", pkg, " package. Install it with ",
           "BiocManager::install('", pkg, "').", call. = FALSE)
    }
  }
  txdb <- getExportedValue("TxDb.Hsapiens.UCSC.hg38.knownGene",
                           "TxDb.Hsapiens.UCSC.hg38.knownGene")
  g <- suppressMessages(GenomicFeatures::genes(txdb))
  w <- withCallingHandlers(
    GenomicFeatures::promoters(g, upstream = 3000L, downstream = 3000L),
    warning = function(cnd) {
      if (grepl("out-of-bound", conditionMessage(cnd))) invokeRestart("muffleWarning")
    })
  .noncoding_cache$promoter_windows <- w
  w
}


#' Which Positions Lie Inside the Windows of a Chromosome
#'
#' @param pos Integer positions (possibly none).
#' @param chr Character chromosome, with or without the "chr" prefix.
#' @param windows A GRanges of windows (UCSC-style sequence names).
#' @return One non-missing logical per position; \code{logical(0)} for none.
#' @keywords internal
#' @noRd
.in_promoter_windows <- function(pos, chr, windows) {
  if (length(pos) == 0L) return(logical(0))
  seqname <- if (grepl("^chr", chr)) chr else paste0("chr", chr)
  w <- windows[as.character(GenomicRanges::seqnames(windows)) == seqname]
  if (length(w) == 0L) return(rep(FALSE, length(pos)))
  q <- IRanges::IRanges(start = as.integer(pos), width = 1L)
  as.logical(IRanges::countOverlaps(q, IRanges::ranges(w)) > 0L)
}


#' The 0-Row Non-Coding Regions Table
#' @keywords internal
#' @noRd
.empty_noncoding_regions <- function(context, n_selected, n_linked) {
  out <- data.frame(region_id = character(0), chr = character(0),
                    start = integer(0), end = integer(0), label = character(0),
                    n_records = integer(0), category = character(0),
                    gene_rank = integer(0), stringsAsFactors = FALSE)
  out$variant_rows <- I(list())
  attr(out, "context") <- context
  attr(out, "n_selected_records") <- as.integer(n_selected)
  attr(out, "n_linked_records") <- as.integer(n_linked)
  out
}
