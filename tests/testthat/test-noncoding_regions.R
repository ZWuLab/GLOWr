# Tests for noncoding_regions.R: define_regions_noncoding() on the fixture
# aGDS with its planted linkage strings (helper-create-test-agds.R).

.agds <- tempfile(fileext = ".gds")
create_test_agds(.agds)
withr::defer(unlink(.agds), teardown_env())

.spec <- variant_filter(rare_maf_cutoff = 0.5, min_variants = 1L)
.idx  <- build_variant_index(.agds)
.pos  <- .idx[["22"]]$pos
.genes <- c("GENE_A", "GENE_B", "GENE_C", "GENE_D", "GENE_E")
.cache <- GLOWr:::.noncoding_cache

.regs <- function(mask, rule = "glow", gene_list = .genes, ...) {
  define_regions_noncoding(.agds, 22, mask, .spec, .idx, rule = rule,
                           gene_list = gene_list, ...)
}
.rows_of <- function(r, gene) r$variant_rows[[match(gene, r$label)]]
# The fixture's non-PASS records: 3, 8, 18, 43, 53, 63, 73, 85, 115, 123.

test_that("the regions table has the documented shape with unique increasing candidates", {
  r <- .regs("upstream")
  expect_equal(names(r), c("region_id", "chr", "start", "end", "label", "n_records",
                           "category", "gene_rank", "variant_rows"))
  expect_equal(r$label, c("GENE_A", "GENE_B"))
  expect_equal(r$region_id, r$label)
  expect_equal(r$category, c("upstream", "upstream"))
  expect_equal(r$gene_rank, 1:2)
  expect_true(is.list(r$variant_rows))
  for (i in seq_len(nrow(r))) {
    v <- r$variant_rows[[i]]
    expect_true(is.integer(v))
    expect_true(all(diff(v) > 0))
    expect_equal(r$n_records[i], length(v))
    expect_equal(r$start[i], min(.pos[v]))
    expect_equal(r$end[i], max(.pos[v]))
  }
  expect_s3_class(attr(r, "context"), "glow_noncoding_context")
  expect_equal(attr(r, "n_selected_records"), 10L)      # rows 21:30 (all PASS)
  expect_equal(attr(r, "n_linked_records"), 18L)        # 10 to GENE_A, 8 to GENE_B
})

test_that("the many-to-many links, the two-category records and the rules on the fixture", {
  # upstream: rows 21:30 name GENE_A; 21:28 ("GENE_A,GENE_B") also GENE_B; 29:30
  # ("GENE_A;GENE_B", upstream;downstream) are upstream of GENE_A only.
  r <- .regs("upstream")
  expect_equal(.rows_of(r, "GENE_A"), 21:30)
  expect_equal(.rows_of(r, "GENE_B"), 21:28)
  s <- .regs("upstream", rule = "staar")                 # 29:30 join neither mask
  expect_equal(.rows_of(s, "GENE_A"), 21:28)
  expect_equal(.rows_of(s, "GENE_B"), 21:28)
  # downstream: 61:70 minus the non-PASS 63 name GENE_B and GENE_C; 29:30 add GENE_B.
  r <- .regs("downstream")
  expect_equal(.rows_of(r, "GENE_B"), c(29:30, 61:62, 64:70))
  expect_equal(.rows_of(r, "GENE_C"), c(61:62, 64:70))
  s <- .regs("downstream", rule = "staar")
  expect_equal(.rows_of(s, "GENE_B"), c(61:62, 64:70))
  # UTR: 71:80 minus 73; the multi-gene 76:78 and the UTR5;UTR3 79:80 also name GENE_C.
  r <- .regs("UTR")
  expect_equal(.rows_of(r, "GENE_B"), c(71:72, 74:80))
  expect_equal(.rows_of(r, "GENE_C"), 76:80)
  s <- .regs("UTR", rule = "staar")                      # first gene only
  expect_equal(s$label, "GENE_B")
  expect_equal(.rows_of(s, "GENE_B"), c(71:72, 74:80))
  # ncRNA: the bundled chromosome 22 table holds POTEH-AS1 and LINC01297; row 40
  # (ncRNA_exonic;splicing, "POTEH-AS1;LINC01297(...)") links both under glow,
  # the first group only under staar.
  r <- .regs("ncRNA", gene_list = NULL)
  expect_equal(r$label, c("POTEH-AS1", "LINC01297"))
  expect_equal(.rows_of(r, "POTEH-AS1"), 36:40)
  expect_equal(.rows_of(r, "LINC01297"), 40L)
  expect_true(all(r$gene_rank == match(r$label, GLOWr:::.noncoding_gene_list("ncrna", "22"))))
  s <- .regs("ncRNA", rule = "staar", gene_list = NULL)
  expect_equal(s$label, "POTEH-AS1")
})

test_that("the enhancer masks need both annotations and link the first GeneHancer target", {
  # CAGE on 22:26 and 69, GeneHancer on 21:23, 61:63 and 69 (a header without targets).
  r <- .regs("enhancer_CAGE")
  expect_equal(attr(r, "n_selected_records"), 3L)        # 22, 23, 69
  expect_equal(r$label, "GENE_A")
  expect_equal(.rows_of(r, "GENE_A"), 22:23)
  # DHS on 21, 30, 31:35, 61, 70; with GeneHancer: 21 (GENE_A first) and 61 (GENE_C first).
  r <- .regs("enhancer_DHS")
  expect_equal(r$label, c("GENE_A", "GENE_C"))
  expect_equal(.rows_of(r, "GENE_A"), 21L)
  expect_equal(.rows_of(r, "GENE_C"), 61L)
  expect_equal(.regs("enhancer_DHS", rule = "staar")$label, c("GENE_A", "GENE_C"))
})

test_that("the promoter masks intersect the selection with the windows by subset-assignment", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")
  # Of the five CAGE records 22:26, windows that hold positions 1 and 5 of the
  # selection (rows 22 and 26): a recycled logical would keep 22 and drop 26.
  w <- GenomicRanges::GRanges(c("chr22", "chr22", "chr1"),
                              IRanges::IRanges(c(.pos[22], .pos[26], 5L), c(.pos[22], .pos[26], 9L)))
  assign("promoter_windows", w, envir = .cache)
  on.exit(rm("promoter_windows", envir = .cache), add = TRUE)
  r <- .regs("promoter_CAGE")
  expect_equal(r$label, "GENE_A")                         # first gene of "GENE_A,GENE_B"
  expect_equal(.rows_of(r, "GENE_A"), c(22L, 26L))
  expect_equal(attr(r, "n_selected_records"), 2L)
  # DHS records 21, 30, 31:35; windows on rows 21 and 35.
  w2 <- GenomicRanges::GRanges(c("chr22", "chr22"),
                               IRanges::IRanges(c(.pos[21], .pos[35]), c(.pos[21], .pos[35])))
  assign("promoter_windows", w2, envir = .cache)
  r <- .regs("promoter_DHS")
  expect_equal(.rows_of(r, "GENE_A"), c(21L, 35L))
  # No window on the chromosome: nothing selected, an empty table of the right shape.
  assign("promoter_windows", GenomicRanges::GRanges("chr1", IRanges::IRanges(1L, 2L)), envir = .cache)
  r <- .regs("promoter_DHS")
  expect_equal(nrow(r), 0L)
  expect_equal(names(r), c("region_id", "chr", "start", "end", "label", "n_records",
                           "category", "gene_rank", "variant_rows"))
  expect_equal(attr(r, "n_selected_records"), 0L)
  # The window helper: one logical per position, none missing, none for no position.
  expect_equal(GLOWr:::.in_promoter_windows(integer(0), "22", w), logical(0))
  expect_equal(GLOWr:::.in_promoter_windows(c(.pos[22], .pos[23], .pos[26]), "22", w),
               c(TRUE, FALSE, TRUE))
  expect_equal(GLOWr:::.in_promoter_windows(c(.pos[22]), "chr22", w), TRUE)
  expect_equal(GLOWr:::.in_promoter_windows(c(.pos[22]), "21", w), FALSE)
})

test_that("the context: QC and type pass, presence as logicals, reuse gives the same table", {
  g <- SeqArray::seqOpen(.agds, readonly = TRUE)
  on.exit(SeqArray::seqClose(g))
  ctx <- GLOWr:::.noncoding_chr_context(g, .spec, "annotation/info/FunctionalAnnotation",
                                        get("annotation_name_catalog", envir = asNamespace("GLOWr")))
  expect_s3_class(ctx, "glow_noncoding_context")
  expect_equal(ctx$chr, "22")
  expect_equal(ctx$n, 125L)
  expect_false(ctx$pass[3]); expect_true(ctx$pass[4])      # the non-PASS record 3
  for (p in c("cage", "dhs", "genehancer")) {
    expect_type(ctx[[p]], "logical"); expect_false(anyNA(ctx[[p]]))
  }
  expect_equal(which(ctx$cage), c(22:26, 64:66, 69))
  expect_equal(which(ctx$genehancer), c(21:23, 61:63, 69))
  expect_equal(ctx$category[29], "upstream;downstream")
  r1 <- .regs("UTR")
  r2 <- define_regions_noncoding(.agds, 22, "UTR", .spec, .idx, context = ctx, gene_list = .genes)
  expect_equal(r1[, names(r1) != "variant_rows"], r2[, names(r2) != "variant_rows"])
  expect_equal(r1$variant_rows, r2$variant_rows)
  # A presence value of NA counts as absent (it never recycles or passes).
  ctx_na <- ctx; ctx_na$cage[22] <- NA
  r <- define_regions_noncoding(.agds, 22, "enhancer_CAGE", .spec, .idx, context = ctx_na, gene_list = .genes)
  expect_equal(.rows_of(r, "GENE_A"), 23L)
  # The variant-type setting is applied here as in the extraction.
  ctx_ind <- GLOWr:::.noncoding_chr_context(g, variant_filter(variant_type = "Indel", rare_maf_cutoff = 0.5),
                                            "annotation/info/FunctionalAnnotation",
                                            get("annotation_name_catalog", envir = asNamespace("GLOWr")))
  expect_false(any(ctx_ind$pass))                          # the fixture holds SNVs only
})

test_that("gene_list restricts without fallback and the gene order is the table's", {
  r <- .regs("upstream", gene_list = "GENE_B")
  expect_equal(r$label, "GENE_B")
  expect_equal(.rows_of(r, "GENE_B"), 21:28)
  r <- .regs("upstream", gene_list = c("GENE_B", "GENE_A"))
  expect_equal(r$label, c("GENE_B", "GENE_A"))
  expect_equal(r$gene_rank, 1:2)
  r <- .regs("upstream", gene_list = "GENE_Z")
  expect_equal(nrow(r), 0L)
  expect_equal(attr(r, "n_selected_records"), 10L)
  expect_equal(attr(r, "n_linked_records"), 0L)
})

test_that("input validation: the mask name, the filter, the chromosome, the index, a two-chromosome GDS", {
  expect_error(.regs("plof"), "must be one of")
  expect_error(define_regions_noncoding(.agds, 22, "UTR", list(a = 1)), "glow_variant_filter")
  expect_error(define_regions_noncoding(.agds, 21, "UTR", .spec, .idx), "holds chromosome 22")
  bad_idx <- .idx; bad_idx[["22"]]$pos <- bad_idx[["22"]]$pos[1:10]
  expect_error(define_regions_noncoding(.agds, 22, "UTR", .spec, bad_idx), "does not describe")
  expect_error(define_regions_noncoding(.agds, 22, "UTR", .spec, .idx, context = list()), "context")
  # Two chromosomes in one file.
  two <- tempfile(fileext = ".gds")
  file.copy(.agds, two)
  on.exit(unlink(two), add = TRUE)
  g <- gdsfmt::openfn.gds(two, readonly = FALSE)
  gdsfmt::delete.gdsn(gdsfmt::index.gdsn(g, "chromosome"))   # compressed: replace, not overwrite
  gdsfmt::add.gdsn(g, "chromosome", val = c(rep("22", 100L), rep("21", 25L)),
                   compress = "LZMA_ra", closezip = TRUE)
  gdsfmt::closefn.gds(g)
  expect_error(define_regions_noncoding(two, 22, "UTR", .spec), "one chromosome")
})

test_that("a candidate-list region extracts the same set as its interval, columns aligned", {
  up <- .regs("upstream")
  vs <- extract_variant_set(.agds, up[up$label == "GENE_A", ], .spec, variant_index = .idx, verbose = 0)
  vi <- extract_variant_set(.agds, list(chr = "22", start = min(.pos[21:30]), end = max(.pos[21:30])),
                            .spec, variant_index = .idx, verbose = 0)
  expect_equal(vs$n_total_in_region, 10L)
  expect_equal(vs$variant_info, vi$variant_info)
  expect_equal(vs$G, vi$G)
  # A subset of the candidates (record 25 fails the cohort-MAC filter, so two
  # of the three survive): the genotype columns follow the variant ids.
  sub <- extract_variant_set(.agds, list(chr = "22", start = 1L, end = 2L, label = "sub",
                                         variant_rows = c(21L, 25L, 30L)),
                             .spec, variant_index = .idx, verbose = 0)
  expect_equal(sub$n_total_in_region, 3L)
  expect_equal(sub$variant_info$variant_id, intersect(c(21L, 25L, 30L), vi$variant_info$variant_id))
  expect_equal(sub$G, vi$G[, match(sub$variant_info$variant_id, vi$variant_info$variant_id)])
})

# ---- 2026-10-16: the implementation review round ----

test_that("a context built under another filter is refused, and a missing field is named", {
  g <- SeqArray::seqOpen(.agds, readonly = TRUE)
  on.exit(SeqArray::seqClose(g))
  cat_ <- get("annotation_name_catalog", envir = asNamespace("GLOWr"))
  ctx_ind <- GLOWr:::.noncoding_chr_context(g, variant_filter(variant_type = "Indel", rare_maf_cutoff = 0.5),
                                            "annotation/info/FunctionalAnnotation", cat_)
  expect_equal(ctx_ind$qc_pass_value, .spec$qc_pass_value)
  expect_error(define_regions_noncoding(.agds, 22, "UTR", .spec, .idx, context = ctx_ind, gene_list = .genes),
               "built under another filter")
  expect_error(define_regions_noncoding(.agds, 22, "upstream", .spec, .idx, gene_list = .genes,
                                        Annotation_dir = "annotation/info/NoSuchFolder"),
               "no annotation node for the field 'GENCODE.Category'")
})

test_that("the FAVOR 2.0 layout: the same links where the information survives, every gene of a two-category record otherwise", {
  f2 <- tempfile(fileext = ".gds")
  create_test_agds(f2, layout = "favor2")
  on.exit(unlink(f2), add = TRUE)
  g2 <- SeqArray::seqOpen(f2, readonly = TRUE)
  info2 <- as.character(SeqArray::seqGetData(g2, "annotation/info/FunctionalAnnotation/genecode_comprehensive_info"))
  SeqArray::seqClose(g2)
  expect_equal(info2[29], "GENE_A,GENE_B")                    # the ";" group boundary is gone
  expect_equal(info2[79], "GENE_B,GENE_C")                    # the transcript details are gone
  expect_equal(info2[40], "POTEH-AS1,LINC01297")
  idx2 <- build_variant_index(f2)
  regs2 <- function(mask, ...) define_regions_noncoding(f2, 22, mask, .spec, idx2, gene_list = .genes, ...)
  same <- function(a, b) {
    expect_equal(b$label, a$label)
    expect_equal(b$variant_rows, a$variant_rows)
  }
  # The masks whose strings carry the same genes in both layouts.
  same(.regs("UTR"), regs2("UTR"))
  same(.regs("enhancer_CAGE"), regs2("enhancer_CAGE"))
  same(.regs("enhancer_DHS"), regs2("enhancer_DHS"))
  same(.regs("ncRNA", gene_list = NULL),
       define_regions_noncoding(f2, 22, "ncRNA", .spec, idx2, gene_list = NULL))
  # The two-category records 29:30 ("GENE_A,GENE_B" under upstream;downstream)
  # link both genes under both masks: GENE_B gains them upstream, GENE_A downstream.
  up <- regs2("upstream")
  expect_equal(.rows_of(up, "GENE_A"), 21:30)
  expect_equal(.rows_of(up, "GENE_B"), 21:30)
  dn <- regs2("downstream")
  expect_equal(.rows_of(dn, "GENE_A"), 29:30)
  expect_equal(.rows_of(dn, "GENE_B"), c(29:30, 61:62, 64:70))
  expect_equal(.rows_of(dn, "GENE_C"), c(61:62, 64:70))
  # Extraction from a candidate list of the FAVOR 2.0 file: the same set as the
  # interval over those records, the genotype columns aligned to the variant ids.
  pos2 <- idx2[["22"]]$pos
  vs <- extract_variant_set(f2, up[up$label == "GENE_B", ], .spec, variant_index = idx2, verbose = 0)
  vi <- extract_variant_set(f2, list(chr = "22", start = min(pos2[21:30]), end = max(pos2[21:30])),
                            .spec, variant_index = idx2, verbose = 0)
  expect_equal(vs$n_total_in_region, 10L)
  expect_equal(vs$variant_info, vi$variant_info)
  expect_equal(vs$G, vi$G)
})
