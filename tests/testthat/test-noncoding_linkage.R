# Tests for noncoding_linkage.R: the gene-linkage rules of the gene-centric
# non-coding masks, on fixtures of both string layouts (FAVOR v1 strings taken
# from the chromosome 22 census of the planning review, and FAVOR 2.0 strings
# produced by .favor2_serialize() on constructed inputs).

link  <- GLOWr:::.link_gencode
gh    <- GLOWr:::.link_genehancer
sel   <- GLOWr:::.category_select
specs <- GLOWr:::.noncoding_mask_specs()

# ---- the cleanup primitives ----

test_that(".strip_details removes balanced details only, including nested ones and commas inside", {
  strip <- GLOWr:::.strip_details
  expect_equal(strip("POTEH(ENST00000000001.1:c.*1704C>T)"), "POTEH")
  expect_equal(strip("GENE_B(ENST1:c.*100A>G,ENST2:c.*95A>G)"), "GENE_B")
  expect_equal(strip("A(x(y)z),B"), "A,B")
  expect_equal(strip("NONE(dist=NONE),RF00004(dist=216905)"), "NONE,RF00004")
  expect_equal(strip("A(unbalanced,B"), "A(unbalanced,B")
  expect_equal(strip(c("A", NA, "")), c("A", NA, ""))
})

test_that(".gene_groups splits groups and genes, drops NONE and keeps hyphens", {
  groups <- GLOWr:::.gene_groups(c(
    "A;B,C", "POTEH-AS1;POTEH(ENST1)", "NONE(dist=NONE),RF00004(dist=216905)",
    "GGT1(ENST1);GGT1(ENST2)", "", NA))
  expect_equal(groups[[1]], list("A", c("B", "C")))
  expect_equal(groups[[2]], list("POTEH-AS1", "POTEH"))
  expect_equal(groups[[3]], list("RF00004"))
  expect_equal(groups[[4]], list("GGT1", "GGT1"))
  expect_equal(groups[[5]], list())
  expect_equal(groups[[6]], list())
})

# ---- selection ----

test_that(".category_select: glow matches tokens, staar matches the literal strings", {
  cat <- c("upstream", "upstream;downstream", "downstream", "UTR3", "UTR5", "UTR5;UTR3",
           "ncRNA_exonic", "ncRNA_exonic;splicing", "ncRNA_splicing", "exonic;splicing",
           "intronic", "", NA)
  expect_equal(sel(cat, specs$upstream, "glow"),   cat %in% c("upstream", "upstream;downstream"))
  expect_equal(sel(cat, specs$upstream, "staar"),  cat %in% "upstream")
  expect_equal(sel(cat, specs$downstream, "glow"), cat %in% c("downstream", "upstream;downstream"))
  expect_equal(sel(cat, specs$UTR, "glow"),        cat %in% c("UTR3", "UTR5", "UTR5;UTR3"))
  expect_equal(sel(cat, specs$UTR, "staar"),       cat %in% c("UTR3", "UTR5", "UTR5;UTR3"))
  expect_equal(sel(cat, specs$ncRNA, "glow"),
               cat %in% c("ncRNA_exonic", "ncRNA_exonic;splicing", "ncRNA_splicing"))
  expect_equal(sel(cat, specs$ncRNA, "staar"),
               cat %in% c("ncRNA_exonic", "ncRNA_exonic;splicing", "ncRNA_splicing"))
  # "upstreamX" is not the token "upstream"; a presence mask selects everything.
  expect_false(sel("upstreamX", specs$upstream, "glow"))
  expect_equal(sel(cat, specs$enhancer_CAGE, "glow"), rep(TRUE, length(cat)))
  expect_false(anyNA(sel(cat, specs$upstream, "glow")))
})

# ---- the glow rule on v1 strings ----

test_that("glow: upstream/downstream link the group at the category's position", {
  expect_equal(link("A,B", "upstream", "upstream")[[1]], c("A", "B"))
  expect_equal(link("A;B", "upstream;downstream", "upstream")[[1]], "A")
  expect_equal(link("A;B", "upstream;downstream", "downstream")[[1]], "B")
  expect_equal(link("A,B;C", "upstream;downstream", "upstream")[[1]], c("A", "B"))
  expect_equal(link("A;B,C", "upstream;downstream", "downstream")[[1]], c("B", "C"))
  expect_equal(link("", "upstream", "upstream")[[1]], character(0))
  expect_equal(link(NA, "upstream", "upstream")[[1]], character(0))
})

test_that("glow: UTR links all genes of the UTR groups, details with commas never split", {
  expect_equal(link("POTEH(ENST00000000001.1:c.*1704C>T)", "UTR3", "UTR")[[1]], "POTEH")
  expect_equal(link("ESS2(ENST00000252137.11:c.*1704C>T),TSSK2(ENST00000399635.4:c.*16G>A)",
                    "UTR3", "UTR")[[1]], c("ESS2", "TSSK2"))
  expect_equal(link("GENE_B(ENST1:c.*100A>G,ENST2:c.*95A>G)", "UTR3", "UTR")[[1]], "GENE_B")
  expect_equal(link("GGT1(ENST1);GGT1(ENST2)", "UTR5;UTR3", "UTR")[[1]], "GGT1")   # repeated token once
  expect_equal(link("A(ENST1);B(ENST2),C(ENST3)", "UTR5;UTR3", "UTR")[[1]], c("A", "B", "C"))
})

test_that("glow: promoter masks take the first gene of the first group; NONE links nothing", {
  for (m in c("promoter_CAGE", "promoter_DHS")) {
    expect_equal(link("POTEH-AS1;POTEH(ENST1)", "ncRNA_exonic;splicing", m)[[1]], "POTEH-AS1")
    expect_equal(link("GENE_A,GENE_B", "upstream", m)[[1]], "GENE_A")
    expect_equal(link("NONE(dist=NONE),RF00004(dist=216905)", "intergenic", m)[[1]], character(0))
    expect_equal(link("", "", m)[[1]], character(0))
    expect_equal(link("A(x)", "exonic", m)[[1]], "A")
  }
})

test_that("glow: ncRNA links all genes of all groups", {
  expect_equal(link("THAP7-AS1;TUBA3FP(ENST00000422086.5:exon3:c.714-2A>G,ENST00000292748.7:exon4:c.862-2A>G)",
                    "ncRNA_exonic;splicing", "ncRNA")[[1]], c("THAP7-AS1", "TUBA3FP"))
  expect_equal(link("A,B,C,D", "ncRNA_exonic", "ncRNA")[[1]], c("A", "B", "C", "D"))
  expect_equal(link("A(tx:1),B(tx:2)", "ncRNA_exonic", "ncRNA")[[1]], c("A", "B"))
})

# ---- the staar rule reproduces the literal operations ----

test_that("staar: the literal operations, mask by mask", {
  # upstream/downstream: split on "," and link every token (no detail removal)
  expect_equal(link("A,B", "upstream", "upstream", rule = "staar")[[1]], c("A", "B"))
  expect_equal(link("A;B", "upstream;downstream", "upstream", rule = "staar")[[1]], "A;B")
  # UTR: the text before the first "(", so a multi-gene record links its first gene only
  expect_equal(link("ESS2(ENST00000252137.11:c.*1704C>T),TSSK2(ENST00000399635.4:c.*16G>A)",
                    "UTR3", "UTR", rule = "staar")[[1]], "ESS2")
  expect_equal(link("A(ENST1);B(ENST2)", "UTR5;UTR3", "UTR", rule = "staar")[[1]], "A")
  # promoter: split on "(", ",", ";", "-" and take the first token, so a hyphenated name is cut
  expect_equal(link("POTEH-AS1;POTEH(ENST1)", "ncRNA_exonic;splicing", "promoter_CAGE",
                    rule = "staar")[[1]], "POTEH")
  expect_equal(link("NONE(dist=NONE),RF00004(dist=216905)", "intergenic", "promoter_DHS",
                    rule = "staar")[[1]], "NONE")
  expect_equal(link("", "", "promoter_DHS", rule = "staar")[[1]], character(0))
  # ncRNA: first ";" group, one greedy removal from its first "(" to its last ")", three tokens
  expect_equal(link("THAP7-AS1;TUBA3FP(ENST00000422086.5:exon3:c.714-2A>G,ENST00000292748.7:exon4:c.862-2A>G)",
                    "ncRNA_exonic;splicing", "ncRNA", rule = "staar")[[1]], "THAP7-AS1")
  expect_equal(link("A(tx:1),B(tx:2)", "ncRNA_exonic", "ncRNA", rule = "staar")[[1]], "A")
  expect_equal(link("A,B,C,D", "ncRNA_exonic", "ncRNA", rule = "staar")[[1]], c("A", "B", "C"))
  expect_equal(link("", "ncRNA_exonic", "ncRNA", rule = "staar")[[1]], character(0))
})

# ---- GeneHancer ----

test_that("GeneHancer: the first listed connected gene of the first element, both rules", {
  two_targets <- "Name=0.2;genehancer_id=GH22I010683;connected_gene=G1;score=0.23;connected_gene=G2;score=0.15"
  two_elements <- paste0("Name=0.2;genehancer_id=GH22I046625;connected_gene=GRAMD4;score=21.11;",
                         "connected_gene=GC22M046648;score=0.34,0.7;genehancer_id=GH22I046627;",
                         "connected_gene=GRAMD4;score=550.23;connected_gene=GC22M046648;score=0.34")
  header_only <- "Name=0.10;genehancer_id=GH22E000003"
  first_empty <- "Name=0.2;genehancer_id=GHA,0.7;genehancer_id=GHB;connected_gene=X;score=1.00"
  x <- c(two_targets, two_elements, header_only, "", NA)
  expect_equal(gh(x, "glow"),  c("G1", "GRAMD4", NA, NA, NA))
  expect_equal(gh(x, "staar"), c("G1", "GRAMD4", NA, NA, NA))
  # A first element without targets: glow links nothing; the literal fourth
  # "=" token is the second element's id, which no gene table contains.
  expect_equal(gh(first_empty, "glow"), NA_character_)
  expect_equal(gh(first_empty, "staar"), "GHB")
})

# ---- the FAVOR 2.0 layout: the same links where the information survives ----

test_that("FAVOR 2.0 serialization: the glow rule gives the same links, with the all-genes policy", {
  fr <- list(
    .gi_genes = list(c("ESS2", "TSSK2"), c("A", "B"), "POTEH-AS1", c("THAP7-AS1", "TUBA3FP"),
                     c("A", "B", "C", "D")),
    .gh_id = c("GH22I010683", "GH22E000003"), .gh_score = c(0.2, 0.1),
    .gh_targets = list(data.frame(gene = c("G1", "G2"), score = c(0.23, 0.15)),
                       data.frame(gene = character(0), score = numeric(0))))
  s_info <- GLOWr:::.favor2_serialize(fr, 1:5, "genecode_comprehensive_info")$genecode_comprehensive_info
  s_gh   <- GLOWr:::.favor2_serialize(fr, 1:2, "genehancer")$genehancer
  expect_equal(s_info, c("ESS2,TSSK2", "A,B", "POTEH-AS1", "THAP7-AS1,TUBA3FP", "A,B,C,D"))
  # UTR multi-gene record: both genes in both layouts.
  expect_equal(link(s_info[1], "UTR3", "UTR")[[1]], c("ESS2", "TSSK2"))
  # A two-category record loses its group boundary: every gene under each mask.
  expect_equal(link(s_info[2], "upstream;downstream", "upstream")[[1]], c("A", "B"))
  expect_equal(link(s_info[2], "upstream;downstream", "downstream")[[1]], c("A", "B"))
  # Promoter: the first gene; ncRNA: all genes, as in v1.
  expect_equal(link(s_info[3], "ncRNA_exonic", "promoter_CAGE")[[1]], "POTEH-AS1")
  expect_equal(link(s_info[4], "ncRNA_exonic;splicing", "ncRNA")[[1]], c("THAP7-AS1", "TUBA3FP"))
  expect_equal(link(s_info[5], "ncRNA_exonic", "ncRNA")[[1]], c("A", "B", "C", "D"))
  # GeneHancer keeps the v1 layout: the first target, NA for a header without targets.
  expect_equal(s_gh[2], "Name=0.10;genehancer_id=GH22E000003")
  expect_equal(gh(s_gh, "glow"), c("G1", NA))
  expect_equal(gh(s_gh, "staar"), c("G1", NA))
  # The literal UTR rule on the serialized multi-gene string matches no gene.
  expect_equal(link(s_info[1], "UTR3", "UTR", rule = "staar")[[1]], "ESS2,TSSK2")
})

test_that("the mask specifications are complete and consistent", {
  expect_equal(names(specs), c("upstream", "downstream", "UTR", "promoter_CAGE", "promoter_DHS",
                               "enhancer_CAGE", "enhancer_DHS", "ncRNA"))
  for (s in specs) {
    expect_true(s$gene_list %in% c("protein_coding", "ncrna"))
    expect_true(s$linkage %in% c("GENCODE.Info", "GeneHancer"))
    expect_true(is.null(s$presence) || all(s$presence %in% c("cage", "dhs", "genehancer")))
  }
  expect_true(specs$ncRNA$all_groups)
  expect_equal(specs$ncRNA$gene_list, "ncrna")
  expect_error(link("A", "upstream", "nonsense"), "Unknown mask")
})
