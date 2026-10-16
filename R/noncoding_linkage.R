# ==============================================================================
# File: noncoding_linkage.R
# Purpose: the gene-linkage rules of the gene-centric non-coding masks (design
#   4.6). A mask selects records by their GENCODE category or by the presence of
#   a CAGE, DHS or GeneHancer annotation, and links each selected record to the
#   genes its annotation names. The `glow` rule reads each field by its meaning
#   (balanced parenthesized details removed first, ";" between the gene groups
#   of a multi-category record, "," between the genes of a group, "NONE" a
#   sentinel, hyphenated names kept whole), so that FAVOR v1 strings and GLOWr's
#   FAVOR 2.0 serialization give the same links wherever the information is
#   present in both. The `staar` rule reproduces STAARpipeline 0.9.8's literal
#   string operations mask by mask, for the validation gate that enumerates the
#   differences; it is not built from the `glow` cleanup.
#
# INTERNAL HELPERS:
#   - .noncoding_mask_specs()   The eight mask specifications
#   - .category_select()        The records a mask selects, under a rule
#   - .strip_details()          Remove every balanced parenthesized detail
#   - .gene_groups()            A GENCODE.Info string as groups of gene names
#   - .link_gencode()           Genes linked through GENCODE.Info, both rules
#   - .link_gencode_staar()     The literal STAARpipeline GENCODE operations
#   - .first_gene_token()       The promoter masks' link under the glow rule
#   - .link_genehancer()        The first connected gene of a GeneHancer record

#################### INTERNAL HELPER FUNCTIONS ####################

#' The Eight Gene-Centric Non-Coding Mask Specifications
#'
#' One entry per mask: the gene table it draws from (\code{"protein_coding"},
#' the bundled \code{genes_info}, or \code{"ncrna"}, the bundled
#' \code{ncRNA_gene}); the GENCODE categories it selects under the \code{glow}
#' rule (\code{categories}, matched as tokens of the ";"-separated category
#' list) and under the literal rule (\code{literal_categories}, exact
#' equality with the whole category string, as the STAARpipeline function
#' tests it); the annotations that must be present (\code{presence}: the
#' context's logical vectors \code{cage}, \code{dhs}, \code{genehancer}); the
#' linkage field; whether the selection is intersected with the TSS windows
#' (\code{promoter}); and, for the ncRNA mask, that all gene groups of a record
#' are linked (\code{all_groups}).
#'
#' @return A named list of eight specifications.
#' @keywords internal
#' @noRd
.noncoding_mask_specs <- function() {
  gencode <- function(gene_list, categories, literal, all_groups = FALSE)
    list(gene_list = gene_list, categories = categories, literal_categories = literal,
         presence = NULL, linkage = "GENCODE.Info", promoter = FALSE,
         all_groups = all_groups)
  list(
    upstream      = gencode("protein_coding", "upstream", "upstream"),
    downstream    = gencode("protein_coding", "downstream", "downstream"),
    UTR           = gencode("protein_coding", c("UTR3", "UTR5"),
                            c("UTR3", "UTR5", "UTR5;UTR3")),
    promoter_CAGE = list(gene_list = "protein_coding", categories = NULL,
                         literal_categories = NULL, presence = "cage",
                         linkage = "GENCODE.Info", promoter = TRUE, all_groups = FALSE),
    promoter_DHS  = list(gene_list = "protein_coding", categories = NULL,
                         literal_categories = NULL, presence = "dhs",
                         linkage = "GENCODE.Info", promoter = TRUE, all_groups = FALSE),
    enhancer_CAGE = list(gene_list = "protein_coding", categories = NULL,
                         literal_categories = NULL, presence = c("cage", "genehancer"),
                         linkage = "GeneHancer", promoter = FALSE, all_groups = FALSE),
    enhancer_DHS  = list(gene_list = "protein_coding", categories = NULL,
                         literal_categories = NULL, presence = c("dhs", "genehancer"),
                         linkage = "GeneHancer", promoter = FALSE, all_groups = FALSE),
    ncRNA         = gencode("ncrna", c("ncRNA_exonic", "ncRNA_splicing"),
                            c("ncRNA_exonic", "ncRNA_exonic;splicing", "ncRNA_splicing"),
                            all_groups = TRUE)
  )
}


#' The Records a Mask Selects by GENCODE Category
#'
#' Under \code{"glow"} a record is selected when the ";"-separated tokens of
#' its category contain one of the mask's categories (so
#' \code{"upstream;downstream"} joins both the upstream and the downstream
#' mask, and \code{"ncRNA_exonic;splicing"} the ncRNA mask through its first
#' token). Under \code{"staar"} the whole category string must equal one of the
#' literal values the STAARpipeline function tests. A mask without categories
#' (the promoter and enhancer masks) selects every record here; their presence
#' vectors are applied by the caller.
#'
#' @param category Character vector of GENCODE.Category values (NA as empty).
#' @param spec One entry of \code{\link{.noncoding_mask_specs}}.
#' @param rule \code{"glow"} or \code{"staar"}.
#' @return A logical vector, one value per record, never NA.
#' @keywords internal
#' @noRd
.category_select <- function(category, spec, rule) {
  if (is.null(spec$categories)) return(rep(TRUE, length(category)))
  category <- as.character(category)
  category[is.na(category)] <- ""
  if (identical(rule, "staar")) return(category %in% spec$literal_categories)
  pattern <- paste0("(^|;)(", paste(spec$categories, collapse = "|"), ")(;|$)")
  grepl(pattern, category)
}


#' Remove Every Balanced Parenthesized Detail From a String
#'
#' The details of a GENCODE.Info string (transcript ids, exon numbers, HGVS
#' text, distances) sit in parentheses and contain commas, so they are removed
#' before any split on ",". Innermost pairs are removed first, repeatedly, so a
#' nested pair is handled too; an unbalanced parenthesis is left as it is.
#'
#' @param x Character vector (NA stays NA).
#' @return Character vector of the same length.
#' @keywords internal
#' @noRd
.strip_details <- function(x) {
  x <- as.character(x)
  repeat {
    y <- gsub("\\([^()]*\\)", "", x)
    if (identical(y, x)) break
    x <- y
  }
  x
}


#' Split GENCODE.Info Strings Into Groups of Gene Names
#'
#' Details removed, ";" separates the groups (one per token of the category
#' list), "," the genes inside a group. The \code{NONE} sentinel of intergenic
#' records and empty tokens are dropped; hyphenated names stay whole.
#'
#' @param info Character vector of GENCODE.Info values.
#' @return A list (one element per record) of lists of character vectors (one
#'   per group), possibly empty.
#' @keywords internal
#' @noRd
.gene_groups <- function(info) {
  info <- .strip_details(info)
  info[is.na(info)] <- ""
  lapply(strsplit(info, ";", fixed = TRUE), function(groups) {
    lapply(groups, function(g) {
      tok <- trimws(strsplit(g, ",", fixed = TRUE)[[1]])
      tok[nzchar(tok) & tok != "NONE"]
    })
  })
}


#' Genes Linked to Each Record Through GENCODE.Info
#'
#' @param info Character, the GENCODE.Info of the selected records.
#' @param category Character, their GENCODE.Category (same length).
#' @param mask The mask name (its specification gives the categories and the
#'   literal rule).
#' @param rule \code{"glow"}: the groups at the positions of the mask's
#'   categories in the ";"-split category list (all groups for the ncRNA mask);
#'   when the info has fewer groups than the category list has tokens (the
#'   FAVOR 2.0 layout, where the group boundaries are not serialized), every
#'   gene of the list, which is that layout's documented ambiguity policy; for
#'   the promoter masks the first gene of the first group, and a \code{NONE}
#'   first token links to no gene. \code{"staar"}: the literal operations of
#'   the mask's STAARpipeline function (\code{\link{.link_gencode_staar}}).
#' @return A list of character vectors, one per record, each possibly empty.
#' @keywords internal
#' @noRd
.link_gencode <- function(info, category, mask, rule = c("glow", "staar")) {
  rule <- match.arg(rule)
  spec <- .noncoding_mask_specs()[[mask]]
  if (is.null(spec)) stop("Unknown mask: '", mask, "'.", call. = FALSE)
  stopifnot(length(info) == length(category))
  if (rule == "staar") return(.link_gencode_staar(info, mask))
  if (isTRUE(spec$promoter)) return(.first_gene_token(info))
  groups <- .gene_groups(info)
  if (isTRUE(spec$all_groups)) {
    return(lapply(groups, function(g) unique(as.character(unlist(g)))))
  }
  category <- as.character(category)
  category[is.na(category)] <- ""
  tokens <- strsplit(category, ";", fixed = TRUE)
  mapply(function(g, ct) {
    if (length(g) == 0L) return(character(0))
    if (length(g) < length(ct)) return(unique(as.character(unlist(g))))   # FAVOR 2.0 layout
    at <- which(ct %in% spec$categories)
    at <- at[at <= length(g)]
    unique(as.character(unlist(g[at])))
  }, groups, tokens, SIMPLIFY = FALSE, USE.NAMES = FALSE)
}


#' The Literal STAARpipeline GENCODE.Info Operations, Mask by Mask
#'
#' Reproduces \code{legacy-materials/STAARpipeline/R/} 0.9.8 for the
#' validation gate: upstream and downstream split the info on "," and link
#' every token (\code{upstream.R:62-74}); UTR takes the text before the first
#' "(" (\code{UTR.R:61-66}); the promoter masks split on "(", ",", ";" and "-"
#' and take the first token, so a hyphenated name is cut
#' (\code{promoter_CAGE.R:35-36}); ncRNA takes the first ";" group, removes
#' everything from its first "(" to its last ")" in one greedy substitution,
#' and links the first three "," tokens (\code{ncRNA.R:113-130}). Where the
#' native code would fail on an empty string (the promoter \code{[[1]]} of an
#' empty split) this returns no gene.
#'
#' @param info Character, the GENCODE.Info of the selected records.
#' @param mask The mask name.
#' @return A list of character vectors, one per record.
#' @keywords internal
#' @noRd
.link_gencode_staar <- function(info, mask) {
  info <- as.character(info)
  info[is.na(info)] <- ""
  drop_na <- function(x) x[!is.na(x)]
  switch(mask,
    upstream = , downstream =
      lapply(strsplit(info, "[,]"), function(z) drop_na(z)),
    UTR =
      lapply(strsplit(info, "[(]"), function(z) drop_na(z[1L])),
    promoter_CAGE = , promoter_DHS =
      lapply(strsplit(info, "\\(|\\,|;|-"), function(z) drop_na(z[1L])),
    ncRNA =
      lapply(strsplit(info, "[;]"), function(z) {
        g <- gsub("\\(.*\\)", "", z[1L])
        drop_na(strsplit(g, ",")[[1L]][1:3])
      }),
    stop("No literal GENCODE rule for mask '", mask, "'.", call. = FALSE))
}


#' The Promoter Masks' Link Under the glow Rule
#'
#' The gene the variant is annotated to: the first gene of the first group,
#' details removed, the hyphenated name kept whole. A \code{NONE} or empty
#' first token links to no gene, with no fallback to the next token.
#'
#' @param info Character, the GENCODE.Info of the selected records.
#' @return A list of character vectors (length 0 or 1), one per record.
#' @keywords internal
#' @noRd
.first_gene_token <- function(info) {
  info <- .strip_details(info)
  info[is.na(info)] <- ""
  first_group <- vapply(strsplit(info, ";", fixed = TRUE),
                        function(z) if (length(z)) z[1L] else "", character(1))
  first_token <- vapply(strsplit(first_group, ",", fixed = TRUE),
                        function(z) if (length(z)) trimws(z[1L]) else "", character(1))
  lapply(first_token, function(t) if (nzchar(t) && t != "NONE") t else character(0))
}


#' The First Listed Connected Gene of a GeneHancer Record
#'
#' A GeneHancer string reads
#' \code{Name=<score>;genehancer_id=<id>;connected_gene=<gene>;score=<score>;...}
#' with the targets in the element's order; a second element, rare, starts
#' after \code{,<score>;genehancer_id=}. Under \code{"glow"} the string is
#' cut at that separator and the first \code{connected_gene} of the first
#' element is returned (design D5); a header without targets gives NA. Under
#' \code{"staar"} the fourth "=" token is cut at the next ";"
#' (\code{enhancer_CAGE.R:38-39}), which is the same gene on every record with
#' a target; a string with fewer than four "=" tokens, where the native code
#' fails, gives NA.
#'
#' @param gh Character vector of GeneHancer values.
#' @param rule \code{"glow"} or \code{"staar"}.
#' @return A character vector of the same length, NA where no gene is linked.
#' @keywords internal
#' @noRd
.link_genehancer <- function(gh, rule = c("glow", "staar")) {
  rule <- match.arg(rule)
  gh <- as.character(gh)
  if (rule == "staar") {
    return(vapply(gh, function(s) {
      if (is.na(s)) return(NA_character_)
      tok <- strsplit(s, "=", fixed = TRUE)[[1L]]
      if (length(tok) < 4L) return(NA_character_)
      strsplit(tok[4L], ";", fixed = TRUE)[[1L]][1L]
    }, character(1), USE.NAMES = FALSE))
  }
  first_element <- sub(",[0-9.]+;genehancer_id=.*$", "", gh)
  out <- rep(NA_character_, length(gh))
  hit <- !is.na(first_element) & grepl("connected_gene=[^;]+", first_element)
  if (any(hit)) {
    m <- regmatches(first_element[hit], regexpr("connected_gene=[^;]+", first_element[hit]))
    out[hit] <- sub("^connected_gene=", "", m)
  }
  out
}
