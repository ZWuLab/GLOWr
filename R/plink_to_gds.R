########## PLINK to GDS Conversion Utility ##########
#
# This file provides a wrapper function for converting PLINK binary format
# (.bed/.bim/.fam) to Genomic Data Structure (GDS) format, which serves as
# the canonical data format for GLOWpipeline.
#
# EXPORTED FUNCTIONS:
#   - plink_to_gds(): Convert PLINK to GDS with optional chromosome splitting
#
# INTERNAL FUNCTIONS:
#   - .get_chromosomes_from_bim(): Parse .bim file for unique chromosomes
#   - .convert_chromosome_names(): Apply 23->X, 24->Y conversion (matches seqBED2GDS)
#   - .validate_regions(): Validate and normalize region specification
#   - .split_gds_by_chromosome(): Split GDS by chromosome using SeqArray
#   - .filter_gds_by_regions(): Filter GDS to specific genomic regions
#   - .extract_prefix_from_path(): Extract file prefix from path

#################### EXPORTED MAIN FUNCTIONS ####################

#' Convert PLINK Binary Files to GDS Format
#'
#' @description
#' Converts PLINK binary format files (.bed/.bim/.fam) to Genomic Data Structure
#' (GDS) format, which serves as the unified data format for GLOWpipeline.
#' Supports optional chromosome-based splitting or region filtering for
#' parallelized downstream analysis (e.g., FAVOR annotation).
#'
#' @param plink_prefix Character. Path to PLINK files without extension.
#'   The function expects three files: \code{{prefix}.bed}, \code{{prefix}.bim},
#'   \code{{prefix}.fam}. If the user provides a path with extension
#'   (e.g., "data.bed"), the extension is stripped automatically.
#' @param output_gds Character. Path for the output GDS file, or a directory path.
#'   If a directory is provided in default mode (no splitting/filtering), the
#'   function auto-generates the filename from the PLINK prefix (e.g., providing
#'   "output/" with plink_prefix "data/mydata" creates "output/mydata.gds").
#'   In split/filter modes, this is always treated as a directory path.
#' @param split_by_chr Logical. If TRUE, create separate GDS files per chromosome.
#'   In this case, \code{output_gds} is treated as a directory path. Default: FALSE.
#' @param chromosomes Character or numeric vector. Specific chromosomes to include.
#'   Only used when \code{split_by_chr = TRUE}. NULL (default) includes all
#'   chromosomes found in the input data.
#' @param regions List. Genomic regions to filter. Mutually exclusive with
#'   \code{split_by_chr}. See Details for format specification. Default: NULL.
#' @param output_prefix Character. Override the file prefix for output files.
#'   Default: NULL (extracted from \code{plink_prefix}).
#' @param chr.conv Logical. If TRUE (default), convert chromosome codes:
#'   23 -> X, 24 -> Y, 25 -> XY, 26 -> MT. This matches the default behavior
#'   of \code{\link[SeqArray]{seqBED2GDS}}. When \code{split_by_chr = TRUE},
#'   the same conversion is applied to chromosome names from the .bim file
#'   to ensure filtering works correctly.
#' @param parallel Logical, integer, or cluster object. Controls parallel
#'   processing for \code{\link[SeqArray]{seqBED2GDS}}:
#'   \itemize{
#'     \item \code{FALSE} (default): Serial processing
#'     \item \code{TRUE}: Auto-detect number of cores
#'     \item Integer: Specific number of cores
#'     \item Cluster object: For distributed computing
#'   }
#' @param verbose Integer. Verbosity level: 0=silent, 1=messages (default).
#' @param ... Additional arguments passed to \code{\link[SeqArray]{seqBED2GDS}}.
#'   Useful options include:
#'   \itemize{
#'     \item \code{compress.geno}: Compression for genotype data (default "LZMA_RA")
#'     \item \code{compress.annotation}: Compression for annotations (default "LZMA_RA")
#'     \item \code{optimize}: Optimize GDS access efficiency (default TRUE)
#'     \item \code{digest}: Add hash codes for integrity (default TRUE)
#'   }
#'
#' @return
#' \itemize{
#'   \item Default mode: Character string with path to created GDS file.
#'   \item Split/region mode: Character vector of created file paths.
#' }
#'
#' @details
#' \strong{Splitting and Filtering:}
#'
#' This function supports three modes of operation:
#'
#' \enumerate{
#'   \item \strong{Default mode}: Convert PLINK to a single GDS file (original behavior).
#'   \item \strong{Split by chromosome}: Create per-chromosome GDS files for parallel processing.
#'   \item \strong{Region filter}: Create GDS file(s) for specific genomic regions.
#' }
#'
#' \code{split_by_chr} and \code{regions} are mutually exclusive; specifying both
#' raises an error.
#'
#' \strong{Region Specification:}
#'
#' The \code{regions} parameter accepts two formats:
#'
#' Single region (shorthand):
#' \preformatted{regions = list(chr = "21", start = 1000000, end = 5000000)}
#'
#' Multiple regions (list of lists):
#' \preformatted{regions = list(
#'   list(chr = "1", start = 1000000, end = 2000000),
#'   list(chr = "2", start = 5000000, end = 6000000)
#' )}
#'
#' \strong{Output Naming:}
#'
#' \tabular{ll}{
#'   Mode \tab Naming Pattern \cr
#'   Default \tab \code{output_gds} (user-specified path) \cr
#'   Split by chr \tab \code{{output_dir}/{prefix}_chr{CHR}.gds} \cr
#'   Region filter \tab \code{{output_dir}/{prefix}_chr{CHR}_{START}_{END}.gds}
#' }
#'
#' \strong{Chromosome Naming:}
#'
#' Chromosome names are preserved from the input .bim file. If input uses "X",
#' output uses "chrX"; if input uses "23", output uses "chr23".
#'
#' PLINK convention reference:
#' \tabular{lll}{
#'   Numeric \tab Letter \tab Description \cr
#'   1-22 \tab 1-22 \tab Autosomes \cr
#'   23 \tab X \tab X chromosome \cr
#'   24 \tab Y \tab Y chromosome \cr
#'   25 \tab XY \tab Pseudo-autosomal region \cr
#'   26 \tab MT \tab Mitochondrial
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic conversion (unchanged from original)
#' gds_path <- plink_to_gds(
#'   plink_prefix = "data/my_data",
#'   output_gds = "data/my_data.gds"
#' )
#'
#' # Example 2: Split by chromosome for parallel annotation
#' gds_files <- plink_to_gds(
#'   plink_prefix = "data/my_data",
#'   output_gds = "output/",
#'   split_by_chr = TRUE
#' )
#' # Creates: output/my_data_chr1.gds, output/my_data_chr2.gds, ...
#'
#' # Example 3: Split specific chromosomes only
#' gds_files <- plink_to_gds(
#'   plink_prefix = "data/my_data",
#'   output_gds = "output/",
#'   split_by_chr = TRUE,
#'   chromosomes = c(21, 22)
#' )
#' # Creates: output/my_data_chr21.gds, output/my_data_chr22.gds
#'
#' # Example 4: Filter to a specific genomic region
#' gds_path <- plink_to_gds(
#'   plink_prefix = "data/my_data",
#'   output_gds = "output/",
#'   regions = list(chr = "21", start = 1000000, end = 5000000)
#' )
#' # Creates: output/my_data_chr21_1000000_5000000.gds
#'
#' # Example 5: Multiple regions
#' gds_files <- plink_to_gds(
#'   plink_prefix = "data/my_data",
#'   output_gds = "output/",
#'   regions = list(
#'     list(chr = "1", start = 1000000, end = 2000000),
#'     list(chr = "2", start = 5000000, end = 6000000)
#'   )
#' )
#'
#' # Example 6: Custom output prefix
#' gds_files <- plink_to_gds(
#'   plink_prefix = "data/my_data",
#'   output_gds = "output/",
#'   split_by_chr = TRUE,
#'   output_prefix = "als_study"
#' )
#' # Creates: output/als_study_chr1.gds, output/als_study_chr2.gds, ...
#' }
#'
#' @references
#' Zheng X, et al. (2012). A high-performance computing toolset for relatedness
#' and principal component analysis of SNP data. Bioinformatics 28(24):3326-3328.
#' DOI: 10.1093/bioinformatics/bts606
#'
#' @seealso
#' \code{\link[SeqArray]{seqBED2GDS}} for the underlying conversion function
#' \code{\link{prepare_PI_control_data}} for control variant preparation using GDS
#'
#' For HPC batch annotation, see \code{inst/scripts/annotate_favor_hpc.sh}
#'
#' @export
plink_to_gds <- function(plink_prefix,
                         output_gds,
                         split_by_chr = FALSE,
                         chromosomes = NULL,
                         regions = NULL,
                         output_prefix = NULL,
                         chr.conv = TRUE,
                         parallel = FALSE,
                         verbose = 1,
                         ...) {

  # ========== Step 1: Validate Inputs ==========

  if (!is.character(plink_prefix) || length(plink_prefix) != 1) {
    stop("plink_prefix must be a single character string")
  }

  if (!is.character(output_gds) || length(output_gds) != 1) {
    stop("output_gds must be a single character string")
  }

  if (!is.logical(split_by_chr) || length(split_by_chr) != 1) {
    stop("split_by_chr must be a single logical value")
  }

  # ========== Step 2: Check Mutual Exclusivity ==========

  # split_by_chr and regions cannot both be specified
  if (split_by_chr && !is.null(regions)) {
    stop("split_by_chr and regions are mutually exclusive. ",
         "Use split_by_chr=TRUE for chromosome splitting, or ",
         "regions for genomic region filtering, but not both.")
  }

  # ========== Step 3: Resolve PLINK Paths ==========

  # Remove extension if user provided one (.bed, .bim, or .fam)
  plink_prefix_clean <- sub("\\.(bed|bim|fam)$", "", plink_prefix, ignore.case = TRUE)

  # Construct expected file paths
  bed_file <- paste0(plink_prefix_clean, ".bed")
  bim_file <- paste0(plink_prefix_clean, ".bim")
  fam_file <- paste0(plink_prefix_clean, ".fam")

  # ========== Step 4: Validate File Existence ==========

  # Check which files exist
  files_exist <- c(
    bed = file.exists(bed_file),
    bim = file.exists(bim_file),
    fam = file.exists(fam_file)
  )

  # If any file is missing, provide informative error
  if (!all(files_exist)) {
    missing_files <- names(files_exist)[!files_exist]
    missing_paths <- c(
      bed = bed_file,
      bim = bim_file,
      fam = fam_file
    )[missing_files]

    stop(sprintf(
      "Missing required PLINK files:\n  %s\n\nExpected prefix: %s\nPlease ensure all three files (.bed, .bim, .fam) exist.",
      paste(paste0("- ", missing_paths), collapse = "\n  "),
      plink_prefix_clean
    ))
  }

  # ========== Step 5: Determine Output Mode and Validate ==========

  # Determine if we're in splitting/filtering mode
  is_split_mode <- split_by_chr || !is.null(regions)

  if (is_split_mode) {
    # In split/filter mode, output_gds should be a directory
    # Ensure directory exists (create if needed is not our job - user should create)
    if (!dir.exists(output_gds)) {
      # Try to interpret as directory path and check parent exists
      if (dir.exists(dirname(output_gds))) {
        # Create the directory
        dir.create(output_gds, recursive = FALSE)
        if (verbose >= 1) {
          message(sprintf("Created output directory: %s", output_gds))
        }
      } else {
        stop(sprintf("Output directory does not exist and parent is also missing: %s",
                     output_gds))
      }
    }
    output_dir <- output_gds
  } else {
    # Standard mode: output_gds can be a file path OR a directory
    # If directory, auto-generate filename from PLINK prefix

    if (dir.exists(output_gds)) {
      # output_gds is a directory - auto-generate filename
      auto_prefix <- .extract_prefix_from_path(plink_prefix_clean)
      output_gds <- file.path(output_gds, paste0(auto_prefix, ".gds"))
      if (verbose >= 1) {
        message(sprintf("Output directory provided; using filename: %s", basename(output_gds)))
      }
    }

    output_dir <- dirname(output_gds)
    if (output_dir != "." && !dir.exists(output_dir)) {
      stop(sprintf("Output directory does not exist: %s", output_dir))
    }

    # Warn if output file already exists (will be overwritten)
    if (file.exists(output_gds) && verbose >= 1) {
      warning(sprintf("Output GDS file already exists and will be overwritten: %s",
                      output_gds), call. = FALSE)
    }
  }

  # ========== Step 6: Determine Output Prefix ==========

  # Extract prefix from plink_prefix if not provided
  if (is.null(output_prefix)) {
    output_prefix <- .extract_prefix_from_path(plink_prefix_clean)
  }

  # ========== Step 7: Handle Split/Filter Modes ==========

  if (split_by_chr) {
    # --- SPLIT BY CHROMOSOME MODE ---

    if (verbose >= 1) {
      message("Converting PLINK to GDS with chromosome splitting...")
      message(sprintf("  Input:  %s (.bed/.bim/.fam)", plink_prefix_clean))
      message(sprintf("  Output: %s/{prefix}_chr{N}.gds", output_dir))
    }

    # Get chromosome list from bim file or user specification
    available_chrs <- .get_chromosomes_from_bim(bim_file)

    if (!is.null(chromosomes)) {
      # User specified chromosomes - validate they exist in data
      chromosomes <- as.character(chromosomes)
      missing_chrs <- setdiff(chromosomes, available_chrs)
      if (length(missing_chrs) > 0) {
        warning(sprintf("Requested chromosomes not found in data: %s",
                        paste(missing_chrs, collapse = ", ")), call. = FALSE)
      }
      target_chrs <- intersect(chromosomes, available_chrs)
      if (length(target_chrs) == 0) {
        stop("None of the requested chromosomes are present in the input data")
      }
    } else {
      target_chrs <- available_chrs
    }

    if (verbose >= 1) {
      message(sprintf("  Chromosomes to process: %s", paste(target_chrs, collapse = ", ")))
    }

    # Convert to temp GDS first
    temp_gds <- tempfile(fileext = ".gds")
    on.exit(unlink(temp_gds), add = TRUE)

    if (verbose >= 1) {
      message("  Step 1: Converting full PLINK to temporary GDS...")
    }

    suppressMessages({
      SeqArray::seqBED2GDS(
        bed.fn = bed_file,
        fam.fn = fam_file,
        bim.fn = bim_file,
        out.gdsfn = temp_gds,
        chr.conv = chr.conv,
        parallel = parallel,
        ...
      )
    })

    # Apply chromosome name conversion to match GDS chromosome names
    target_chrs <- .convert_chromosome_names(target_chrs, chr.conv)

    # Split by chromosome
    if (verbose >= 1) {
      message("  Step 2: Splitting by chromosome...")
    }

    result_files <- .split_gds_by_chromosome(
      gds_path = temp_gds,
      output_dir = output_dir,
      prefix = output_prefix,
      chromosomes = target_chrs,
      verbose = verbose
    )

    if (verbose >= 1) {
      message(sprintf("  Created %d chromosome-specific GDS files", length(result_files)))
    }

    return(result_files)

  } else if (!is.null(regions)) {
    # --- REGION FILTER MODE ---

    # Validate and normalize regions
    regions <- .validate_regions(regions)

    if (verbose >= 1) {
      message("Converting PLINK to GDS with region filtering...")
      message(sprintf("  Input:  %s (.bed/.bim/.fam)", plink_prefix_clean))
      message(sprintf("  Regions: %d region(s)", length(regions)))
    }

    # Convert to temp GDS first
    temp_gds <- tempfile(fileext = ".gds")
    on.exit(unlink(temp_gds), add = TRUE)

    if (verbose >= 1) {
      message("  Step 1: Converting full PLINK to temporary GDS...")
    }

    suppressMessages({
      SeqArray::seqBED2GDS(
        bed.fn = bed_file,
        fam.fn = fam_file,
        bim.fn = bim_file,
        out.gdsfn = temp_gds,
        chr.conv = chr.conv,
        parallel = parallel,
        ...
      )
    })

    # Apply chromosome name conversion to region chromosomes for filtering
    for (i in seq_along(regions)) {
      regions[[i]]$chr <- .convert_chromosome_names(regions[[i]]$chr, chr.conv)
    }

    # Filter by regions
    if (verbose >= 1) {
      message("  Step 2: Filtering by regions...")
    }

    result_files <- .filter_gds_by_regions(
      gds_path = temp_gds,
      regions = regions,
      output_dir = output_dir,
      prefix = output_prefix,
      verbose = verbose
    )

    if (verbose >= 1) {
      message(sprintf("  Created %d region-specific GDS files", length(result_files)))
    }

    return(result_files)

  } else {
    # --- DEFAULT MODE (original behavior) ---

    if (verbose >= 1) {
      message("Converting PLINK to GDS...")
      message(sprintf("  Input:  %s (.bed/.bim/.fam)", plink_prefix_clean))
      message(sprintf("  Output: %s", output_gds))
    }

    # Call seqBED2GDS
    if (verbose == 0) {
      suppressMessages({
        SeqArray::seqBED2GDS(
          bed.fn = bed_file,
          fam.fn = fam_file,
          bim.fn = bim_file,
          out.gdsfn = output_gds,
          chr.conv = chr.conv,
          parallel = parallel,
          ...
        )
      })
    } else {
      SeqArray::seqBED2GDS(
        bed.fn = bed_file,
        fam.fn = fam_file,
        bim.fn = bim_file,
        out.gdsfn = output_gds,
        chr.conv = chr.conv,
        parallel = parallel,
        ...
      )
    }

    # Verify output file was created
    if (!file.exists(output_gds)) {
      stop("GDS conversion failed: output file not created")
    }

    # Report summary
    if (verbose >= 1) {
      gds <- SeqArray::seqOpen(output_gds, readonly = TRUE)
      tryCatch({
        n_variants <- length(SeqArray::seqGetData(gds, "variant.id"))
        n_samples <- length(SeqArray::seqGetData(gds, "sample.id"))
        message(sprintf("  Created GDS with %d variants and %d samples",
                        n_variants, n_samples))
      }, finally = {
        SeqArray::seqClose(gds)
      })
    }

    return(output_gds)
  }
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' Extract Unique Chromosomes from PLINK .bim File
#'
#' @param bim_path Path to .bim file
#' @return Character vector of unique chromosome names
#' @keywords internal
#' @noRd
.get_chromosomes_from_bim <- function(bim_path) {
  # Read only first column (chromosome) from .bim file
  # .bim format: CHR  SNP_ID  CM  POS  A1  A2 (tab-delimited, no header)
  bim_data <- utils::read.table(
    bim_path,
    header = FALSE,
    colClasses = c("character", rep("NULL", 5)),  # Only read first column
    sep = "\t"
  )

  # Get unique chromosomes, maintaining order of appearance
  unique_chrs <- unique(bim_data[[1]])

  return(unique_chrs)
}


#' Convert Chromosome Names to Match seqBED2GDS Conversion
#'
#' When \code{chr.conv = TRUE}, seqBED2GDS converts numeric sex/mito
#' chromosomes to letter codes: 23->X, 24->Y, 25->XY, 26->MT.
#' This function applies the same conversion to chromosome names
#' from .bim files to ensure filtering works correctly.
#'
#' @param chrs Character vector of chromosome names
#' @param chr.conv Logical. If TRUE, apply conversion; if FALSE, return unchanged
#' @return Character vector with converted names
#' @keywords internal
#' @noRd
.convert_chromosome_names <- function(chrs, chr.conv) {
  if (!chr.conv) return(chrs)

  # Apply same conversion as seqBED2GDS chr.conv=TRUE
  mapping <- c("23" = "X", "24" = "Y", "25" = "XY", "26" = "MT")

  # Use vectorized ifelse for conversion
  converted <- ifelse(chrs %in% names(mapping), mapping[chrs], chrs)

  return(converted)
}


#' Validate and Normalize Region Specification
#'
#' @param regions List with region specification (single or list-of-lists)
#' @return Normalized list of lists, each with chr, start, end
#' @keywords internal
#' @noRd
.validate_regions <- function(regions) {
  if (!is.list(regions)) {
    stop("regions must be a list")
  }

  # Check if it's a single region (shorthand format)
  # Shorthand: list(chr = "21", start = 1000000, end = 5000000)
  if (all(c("chr", "start", "end") %in% names(regions))) {
    # Normalize to list-of-lists
    regions <- list(regions)
  }

  # Validate each region
  for (i in seq_along(regions)) {
    region <- regions[[i]]

    if (!is.list(region)) {
      stop(sprintf("Region %d must be a list with chr, start, end", i))
    }

    if (!all(c("chr", "start", "end") %in% names(region))) {
      stop(sprintf("Region %d must have 'chr', 'start', and 'end' elements", i))
    }

    # Validate chr
    if (!is.character(region$chr) && !is.numeric(region$chr)) {
      stop(sprintf("Region %d: chr must be character or numeric", i))
    }
    # Convert to character for consistency
    regions[[i]]$chr <- as.character(region$chr)

    # Validate start and end
    if (!is.numeric(region$start) || length(region$start) != 1) {
      stop(sprintf("Region %d: start must be a single numeric value", i))
    }
    if (!is.numeric(region$end) || length(region$end) != 1) {
      stop(sprintf("Region %d: end must be a single numeric value", i))
    }

    # Ensure start <= end
    if (region$start > region$end) {
      stop(sprintf("Region %d: start (%d) cannot be greater than end (%d)",
                   i, region$start, region$end))
    }

    # Convert to integer for file naming consistency
    regions[[i]]$start <- as.integer(region$start)
    regions[[i]]$end <- as.integer(region$end)
  }

  return(regions)
}


#' Split GDS File by Chromosome
#'
#' @param gds_path Path to source GDS file
#' @param output_dir Output directory
#' @param prefix File prefix for output files
#' @param chromosomes Character vector of chromosomes to extract
#' @param verbose Verbosity level
#' @return Character vector of created file paths
#' @keywords internal
#' @noRd
.split_gds_by_chromosome <- function(gds_path, output_dir, prefix, chromosomes, verbose) {

  # Open source GDS
  gds <- SeqArray::seqOpen(gds_path, readonly = TRUE)
  on.exit(SeqArray::seqClose(gds), add = TRUE)

  result_files <- character(0)

  for (chr in chromosomes) {
    # Construct output path
    output_file <- file.path(output_dir, sprintf("%s_chr%s.gds", prefix, chr))

    if (verbose >= 1) {
      message(sprintf("    Processing chromosome %s...", chr))
    }

    # Set filter to this chromosome
    SeqArray::seqSetFilterChrom(gds, include = chr, verbose = FALSE)

    # Check if any variants match
    n_variants <- length(SeqArray::seqGetData(gds, "variant.id"))
    if (n_variants == 0) {
      if (verbose >= 1) {
        message(sprintf("      Skipping chr%s: no variants", chr))
      }
      # Reset filter for next iteration
      SeqArray::seqResetFilter(gds, verbose = FALSE)
      next
    }

    # Export filtered subset
    SeqArray::seqExport(gds, output_file, verbose = FALSE)

    result_files <- c(result_files, output_file)

    if (verbose >= 1) {
      message(sprintf("      Created %s (%d variants)", basename(output_file), n_variants))
    }

    # Reset filter for next iteration
    SeqArray::seqResetFilter(gds, verbose = FALSE)
  }

  return(result_files)
}


#' Filter GDS File by Genomic Regions
#'
#' @param gds_path Path to source GDS file
#' @param regions Normalized list of region specifications
#' @param output_dir Output directory
#' @param prefix File prefix for output files
#' @param verbose Verbosity level
#' @return Character vector of created file paths
#' @keywords internal
#' @noRd
.filter_gds_by_regions <- function(gds_path, regions, output_dir, prefix, verbose) {

  # Open source GDS
  gds <- SeqArray::seqOpen(gds_path, readonly = TRUE)
  on.exit(SeqArray::seqClose(gds), add = TRUE)

  result_files <- character(0)

  for (i in seq_along(regions)) {
    region <- regions[[i]]

    # Construct output path
    output_file <- file.path(
      output_dir,
      sprintf("%s_chr%s_%d_%d.gds", prefix, region$chr, region$start, region$end)
    )

    if (verbose >= 1) {
      message(sprintf("    Processing region %d: chr%s:%d-%d...",
                      i, region$chr, region$start, region$end))
    }

    # First filter by chromosome
    SeqArray::seqSetFilterChrom(gds, include = region$chr, verbose = FALSE)

    # Then filter by position range
    # Get current positions and filter to range
    positions <- SeqArray::seqGetData(gds, "position")
    in_range <- positions >= region$start & positions <= region$end

    if (sum(in_range) == 0) {
      if (verbose >= 1) {
        message(sprintf("      Skipping: no variants in region"))
      }
      SeqArray::seqResetFilter(gds, verbose = FALSE)
      next
    }

    # Apply position filter by getting current variant IDs and subsetting
    current_variant_ids <- SeqArray::seqGetData(gds, "variant.id")
    target_variant_ids <- current_variant_ids[in_range]

    # Reset and apply combined filter
    SeqArray::seqResetFilter(gds, verbose = FALSE)
    SeqArray::seqSetFilter(gds, variant.id = target_variant_ids, verbose = FALSE)

    n_variants <- length(SeqArray::seqGetData(gds, "variant.id"))

    # Export filtered subset
    SeqArray::seqExport(gds, output_file, verbose = FALSE)

    result_files <- c(result_files, output_file)

    if (verbose >= 1) {
      message(sprintf("      Created %s (%d variants)", basename(output_file), n_variants))
    }

    # Reset filter for next iteration
    SeqArray::seqResetFilter(gds, verbose = FALSE)
  }

  return(result_files)
}


#' Extract File Prefix from Path
#'
#' @param path File path (possibly with directory and extension)
#' @return Base name without extension
#' @keywords internal
#' @noRd
.extract_prefix_from_path <- function(path) {
  # Get basename (removes directory)
  base <- basename(path)
  # Remove any remaining extension
  prefix <- sub("\\.[^.]*$", "", base)
  return(prefix)
}
