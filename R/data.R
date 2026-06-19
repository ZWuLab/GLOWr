########## Package Data Documentation ##########
#
# This file documents the lazy-loaded data objects included with the GLOWr package.
# These datasets are used for examples and demonstrations of B estimation.
#
# EXPORTED DATASETS (documented below, accessible via data()):
#   - ALS_snvs_B_training   Pre-processed ALS GWAS summary statistics for B training
#   - BMD_snvs_B_training   Pre-processed BMD GWAS summary statistics for B training
#
# INTERNAL DATASETS (in R/sysdata.rda, not exported):
#   - genes_info              18,445 protein-coding gene boundaries (autosomes)
#   - ncRNA_gene              21,104 ncRNA gene names by chromosome
#   - annotation_name_catalog 19-row map of annotation names to aGDS node paths
#   See data-raw/build_sysdata.R for full provenance and rebuild instructions.

#' ALS SNVs for B Parameter Training
#'
#' @description
#' Pre-processed GWAS summary statistics for Amyotrophic Lateral Sclerosis (ALS)
#' from literature studies. This dataset has been cleaned and standardized for
#' use in B parameter estimation (effect size as a function of MAF).
#'
#' @format A data frame with 75 rows (variants) and 9 columns:
#' \describe{
#'   \item{rsID}{SNP identifier (character)}
#'   \item{CHR}{Chromosome, autosomes only (numeric 1-22)}
#'   \item{POS}{Genomic position in GRCh37 (numeric)}
#'   \item{MAF}{Minor allele frequency, 0 < MAF \eqn{\le} 0.5 (numeric)}
#'   \item{P}{P-value, 0 < P \eqn{\le} 1 (numeric)}
#'   \item{N}{Sample size, N \eqn{\ge} 1000 (numeric)}
#'   \item{BETA}{Effect size as log odds ratio (numeric)}
#'   \item{TRAIT_TYPE}{"binary" for all variants (character)}
#'   \item{STUDY}{Study accession identifier (character)}
#' }
#'
#' @details
#' ## Processing
#' Prepared from raw ALS literature data using \code{\link{prepare_B_training_data}}
#' with the following quality control:
#' \itemize{
#'   \item Autosomal variants only (CHR 1-22)
#'   \item Minimum sample size: N \eqn{\ge} 1000
#'   \item Non-diagnostic traits excluded (age of onset, survival, C9orf72 interaction)
#'   \item Nicolas A et al. (2018) study excluded (see note below)
#'   \item Duplicates removed (kept variant with largest N)
#'   \item Rows with NA in MAF, P, N, or BETA removed
#' }
#'
#' The Nicolas A study is excluded because its variants overlap substantially with
#' other included studies but report from a different analysis model. This validated
#' 75-SNV dataset produces B estimates with R-squared = 0.50 using the beta method,
#' and correlation > 0.9999 with the original reference formula.
#'
#' Provenance: prepared from the curated 75-SNV ALS set; the raw input
#' \code{inst/extdata/ALS-known-SNPs-raw.xlsx} ships with the package.
#'
#' @source GWAS Catalog \url{https://www.ebi.ac.uk/gwas/}
#'
#' @examples
#' # Load the data
#' data(ALS_snvs_B_training)
#'
#' # Examine structure
#' str(ALS_snvs_B_training)
#' head(ALS_snvs_B_training)
#'
#' \dontrun{
#' # Use in B estimation
#' B <- get_B(
#'   training_trait = "binary",
#'   training_MAF = ALS_snvs_B_training$MAF,
#'   training_BETA = ALS_snvs_B_training$BETA,
#'   target_trait = "binary",
#'   target_MAF = seq(0.01, 0.4, by = 0.05),
#'   target_case_prop = 0.5
#' )
#' }
#'
#' @seealso \code{\link{BMD_snvs_B_training}}, \code{\link{prepare_B_training_data}},
#'   \code{\link{get_B}}
"ALS_snvs_B_training"

#' BMD SNVs for B Parameter Training
#'
#' @description
#' Pre-processed GWAS summary statistics for Bone Mineral Density (BMD) from
#' large-scale literature studies. This dataset has been cleaned and standardized
#' for use in B parameter estimation (effect size as a function of MAF).
#'
#' @format A data frame with 1790 rows (variants) and 10 columns:
#' \describe{
#'   \item{rsID}{SNP identifier (character)}
#'   \item{CHR}{Chromosome (numeric)}
#'   \item{POS}{Genomic position (numeric)}
#'   \item{MAF}{Minor allele frequency, 0 < MAF \eqn{\le} 0.5 (numeric)}
#'   \item{P}{P-value, 0 < P \eqn{\le} 1 (numeric)}
#'   \item{P_mlog10}{Negative log10 p-value, i.e., -log10(P) (numeric)}
#'   \item{N}{Sample size, N \eqn{\ge} 1000 (numeric)}
#'   \item{BETA}{Effect size (numeric, may contain NA)}
#'   \item{TRAIT_TYPE}{"continuous" for all variants (character)}
#'   \item{STUDY}{Study identifier (character)}
#' }
#'
#' @details
#' ## Processing
#' Prepared from raw BMD/osteoporosis literature data (CSV from GWAS Catalog)
#' using \code{\link{prepare_B_training_data}} with the following quality control:
#' \itemize{
#'   \item Three high-quality studies included: "An atlas of genetic influences
#'     on osteoporosis" (Morris 2019), "Identification of 153 new loci",
#'     "Identification of 613 new loci"
#'   \item Binary trait "bone fracture" excluded
#'   \item Minimum sample size: N \eqn{\ge} 1000
#'   \item Rows with NA in MAF, P, or N removed
#'   \item Duplicates NOT removed (retains variants from different studies)
#' }
#'
#' This is a continuous trait dataset. For cross-trait prediction (e.g., BMD to
#' osteoporosis risk), use the p-value method in \code{\link{get_B}}.
#'
#' Provenance: prepared from the curated BMD/osteoporosis literature summary
#' statistics (raw inputs are maintained with the source study, not bundled).
#'
#' @source GWAS Catalog \url{https://www.ebi.ac.uk/gwas/}
#'
#' @references
#' Morris JA, Kemp JP, Youlten SE, et al. (2019). An atlas of genetic influences
#' on osteoporosis in humans and mice. \emph{Nature Genetics}, 51(2), 258-266.
#' \doi{10.1038/s41588-018-0302-x}
#'
#' @examples
#' # Load the data
#' data(BMD_snvs_B_training)
#'
#' # Examine structure
#' str(BMD_snvs_B_training)
#' head(BMD_snvs_B_training)
#'
#' \dontrun{
#' # Use in B estimation (continuous trait, p-value method)
#' B <- get_B(
#'   training_trait = "continuous",
#'   training_MAF = BMD_snvs_B_training$MAF,
#'   training_P = BMD_snvs_B_training$P,
#'   training_N = BMD_snvs_B_training$N,
#'   target_trait = "binary",      # Cross-trait: BMD -> osteoporosis risk
#'   target_MAF = seq(0.01, 0.4, by = 0.05),
#'   target_case_prop = 0.20
#' )
#' }
#'
#' @seealso \code{\link{ALS_snvs_B_training}}, \code{\link{prepare_B_training_data}},
#'   \code{\link{get_B}}
"BMD_snvs_B_training"
