########## Helper Functions for prepare_B_training_data() ##########
#
# This file contains internal helper functions for smart data conversions,
# filtering, and processing in the prepare_B_training_data() pipeline.

#################### SMART CONVERSION HELPERS ####################

#' Extract Sample Size from Text String
#'
#' @description
#' Extracts numeric sample sizes from text descriptions with intelligent group
#' detection. When multiple groups are described (e.g., men and women, cases
#' and controls), the function sums them. Otherwise, it takes the maximum value.
#'
#' @param text_values Character vector with sample size descriptions
#' @param verbose Logical. If TRUE (default), warns about extraction decisions
#'   and potential ambiguities. Set to FALSE to suppress warnings.
#'
#' @return Numeric vector of sample sizes
#'
#' @details
#' \strong{Pattern Detection and Handling:}
#'
#' The function uses keyword detection to intelligently decide whether to sum
#' or take the maximum of multiple numbers found in text:
#'
#' \strong{Summing Indicators} (numbers are summed):
#' \itemize{
#'   \item Gender indicators: "women", "men", "female", "male"
#'   \item Study design: "cases", "controls"
#'   \item Conjunctions: "and" appearing between numbers (but not "and" alone)
#'   \item Arithmetic: "+" symbol
#' }
#'
#' \strong{Examples of Handled Patterns:}
#' \itemize{
#'   \item "76,067 European ancestry women, 66,420 European ancestry men" -> 142,487 (summed)
#'   \item "5,000 cases + 5,000 controls" -> 10,000 (summed)
#'   \item "10,000 cases and 15,000 controls" -> 25,000 (summed)
#'   \item "394,929 European ancestry individuals" -> 394,929 (single value)
#'   \item "Study with 10,000 in 2019" -> 10,000 (year filtered out automatically)
#'   \item "Data from 2020 with 10,000 women, 15,000 men" -> 25,000 (year filtered, genders summed)
#' }
#'
#' \strong{Automatic Filtering:}
#' The function automatically filters out numbers that are unlikely to be sample sizes:
#' \itemize{
#'   \item Years (1900-2099) are excluded from consideration
#'   \item Very small numbers (<50) are excluded
#'   \item This helps avoid incorrectly summing sample sizes with publication years
#' }
#'
#' \strong{Warnings and Ambiguity Detection:}
#'
#' The function warns when:
#' \itemize{
#'   \item Multiple large numbers (>1000) found without clear summing indicators
#'   \item Numbers are being summed based on pattern detection
#'   \item Extraction might be ambiguous
#' }
#'
#' \strong{IMPORTANT LIMITATIONS:}
#' \itemize{
#'   \item \strong{This is a heuristic approach and may not be 100 percent accurate}
#'   \item Complex text descriptions may be misinterpreted
#'   \item Always verify extracted values, especially for critical analyses
#'   \item When in doubt, manually specify sample sizes using the \code{N}
#'         column in your data
#'   \item The function cannot understand semantic context beyond keywords
#' }
#'
#' \strong{Recommendation:} For publication-quality analyses, always manually
#' verify that extracted sample sizes match the original study descriptions.
#'
#' \strong{Computational Complexity:}
#'
#' O(n * m) where n = number of text values, m = average text length.
#' Pattern matching is done with vectorized regex operations where possible.
#'
#' @examples
#' \dontrun{
#' # Gender-stratified sample
#' extract_sample_size_from_text("76,067 women, 66,420 men")
#' # Returns: 142487 (summed with warning)
#'
#' # Case-control study
#' extract_sample_size_from_text("5,000 cases + 5,000 controls")
#' # Returns: 10000 (summed)
#'
#' # Single cohort
#' extract_sample_size_from_text("394,929 European ancestry individuals")
#' # Returns: 394929
#'
#' # Suppress warnings
#' extract_sample_size_from_text("10,000 men, 15,000 women", verbose = FALSE)
#' # Returns: 25000 (summed, no warning)
#' }
#'
#' @keywords internal
#' @noRd
extract_sample_size_from_text <- function(text_values, verbose = TRUE) {

  # Initialize output
  sample_sizes <- rep(NA_real_, length(text_values))

  # If already numeric, return as-is
  if (is.numeric(text_values)) {
    return(text_values)
  }

  # Convert to character
  text_values <- as.character(text_values)

  # Define keywords that indicate multiple groups should be summed
  gender_keywords <- c("women", "men", "female", "male", "males", "females")
  study_design_keywords <- c("cases", "controls", "case", "control")

  for (i in seq_along(text_values)) {
    txt <- text_values[i]

    # Skip NA or empty
    if (is.na(txt) || nchar(trimws(txt)) == 0) {
      next
    }

    # Try direct numeric conversion first
    direct <- suppressWarnings(as.numeric(txt))
    if (!is.na(direct)) {
      sample_sizes[i] <- direct
      next
    }

    # Extract all numbers from text (remove commas first)
    txt_clean <- gsub(",", "", txt)
    txt_lower <- tolower(txt)

    # Find all numbers
    numbers <- as.numeric(unlist(regmatches(txt_clean, gregexpr("[0-9]+", txt_clean))))

    if (length(numbers) == 0) {
      # No numbers found
      next
    }

    if (length(numbers) == 1) {
      # Single number found - straightforward
      sample_sizes[i] <- numbers[1]
      next
    }

    # Multiple numbers found - need to decide whether to sum or take max
    should_sum <- FALSE
    sum_reason <- ""
    numbers_to_use <- numbers  # By default, use all numbers

    # Filter out likely years (1900-2099) and small numbers (<50) from consideration
    # These are rarely sample sizes
    likely_sample_sizes <- numbers[numbers >= 50 & (numbers < 1900 | numbers > 2099)]

    # If we filtered out some numbers, use only the likely sample sizes
    if (length(likely_sample_sizes) > 0 && length(likely_sample_sizes) < length(numbers)) {
      numbers_to_use <- likely_sample_sizes
    }

    # If only one number remains after filtering, use it
    if (length(numbers_to_use) == 1) {
      sample_sizes[i] <- numbers_to_use[1]
      next
    }

    # Check for explicit arithmetic operator
    if (grepl("\\+", txt)) {
      should_sum <- TRUE
      sum_reason <- "explicit '+' operator"
    }

    # Check for gender indicators
    if (!should_sum && any(sapply(gender_keywords, function(kw) grepl(kw, txt_lower, fixed = TRUE)))) {
      should_sum <- TRUE
      sum_reason <- "gender stratification indicators (women/men)"
    }

    # Check for study design indicators
    if (!should_sum && any(sapply(study_design_keywords, function(kw) grepl(kw, txt_lower, fixed = TRUE)))) {
      should_sum <- TRUE
      sum_reason <- "study design indicators (cases/controls)"
    }

    # Check for "and" appearing between numbers (more sophisticated check)
    # Only trigger if "and" appears between two large numbers (>1000)
    if (!should_sum && grepl("\\band\\b", txt_lower)) {
      # Check if we have multiple large numbers
      large_numbers <- numbers_to_use[numbers_to_use > 1000]
      if (length(large_numbers) >= 2) {
        should_sum <- TRUE
        sum_reason <- "'and' conjunction between large numbers"
      }
    }

    # Apply decision
    if (should_sum) {
      sample_sizes[i] <- sum(numbers_to_use, na.rm = TRUE)

      # Warn user about summing
      if (verbose) {
        warning(sprintf(
          paste0("Sample size extraction: Summed %d values (%s) to get N=%g.\n",
                 "  Original text: \"%s\"\n",
                 "  Reason: %s\n",
                 "  Please verify this is correct!"),
          length(numbers_to_use),
          paste(format(numbers_to_use, big.mark = ",", scientific = FALSE), collapse = " + "),
          sample_sizes[i],
          txt,
          sum_reason
        ), call. = FALSE)
      }
    } else {
      # Take maximum (default behavior for multiple numbers without clear indicators)
      sample_sizes[i] <- max(numbers_to_use, na.rm = TRUE)

      # Warn if multiple large numbers without clear indicators
      large_numbers <- numbers_to_use[numbers_to_use > 1000]
      if (length(large_numbers) > 1 && verbose) {
        warning(sprintf(
          paste0("Sample size extraction: Multiple large numbers found, taking maximum N=%g.\n",
                 "  Original text: \"%s\"\n",
                 "  Found values: %s\n",
                 "  If these represent separate groups that should be summed, please verify!"),
          sample_sizes[i],
          txt,
          paste(format(numbers_to_use, big.mark = ",", scientific = FALSE), collapse = ", ")
        ), call. = FALSE)
      }
    }
  }

  return(sample_sizes)
}


#' Convert Odds Ratio (OR) to Log Odds Ratio (BETA)
#'
#' @description
#' Converts odds ratios to log odds ratios (beta coefficients) for binary traits
#' using the mathematical transformation BETA = log(OR). This function assumes
#' the input values are correctly identified as odds ratios.
#'
#' @param or_values Numeric vector of odds ratios
#' @param warn Logical. If TRUE (default), warns about problematic values
#'
#' @return Numeric vector of log odds ratios (BETA values)
#'
#' @details
#' \strong{Mathematical Transformation:}
#'
#' The function applies the standard transformation:
#' \deqn{BETA = log(OR)}
#'
#' \strong{Special Cases Handled:}
#' \itemize{
#'   \item OR = 1 -> BETA = 0 (no effect)
#'   \item OR < 0 -> NA with warning (invalid - OR must be positive)
#'   \item OR = 0 -> NA with warning (invalid - OR cannot be zero)
#'   \item OR = Inf -> NA with warning (invalid - infinite OR)
#' }
#'
#' \strong{CRITICAL ASSUMPTIONS:}
#' \itemize{
#'   \item \strong{Input values are correctly identified as odds ratios}
#'   \item OR values are from the appropriate allele (not inverted)
#'   \item The column name/content truly represents OR (not relative risk, hazard ratio, etc.)
#'   \item Units are consistent across all variants
#' }
#'
#' \strong{IMPORTANT LIMITATIONS:}
#' \itemize{
#'   \item \strong{The function cannot verify that inputs are actually OR values}
#'   \item Cannot detect if OR values are for the wrong allele direction
#'   \item Does not validate biological plausibility of OR magnitudes
#'   \item Assumes OR values are correctly scaled (not percentage increases, etc.)
#' }
#'
#' \strong{CRITICAL WARNING - User Responsibility:}
#' \itemize{
#'   \item \strong{VERIFY that your OR column contains true odds ratios}
#'   \item Check that OR values make biological sense (typically 0.1 to 10)
#'   \item Ensure effect allele is consistent across variants
#'   \item If data contains a mix of OR and BETA, separate them before conversion
#'   \item For publication-quality work, manually verify conversion is appropriate
#'   \item Consider whether effect directions (allele coding) are consistent
#' }
#'
#' \strong{Common Failure Scenarios:}
#' \itemize{
#'   \item OR column contains beta values (will produce nonsensical results)
#'   \item OR values are percentage increases (e.g., "20 percent" stored as 20, not 1.2)
#'   \item Mixed units in OR column (some per-allele, some per-SD)
#'   \item Relative risk or hazard ratios mislabeled as OR
#'   \item OR values for wrong allele (reference vs. effect allele confusion)
#' }
#'
#' \strong{Recommendation:} Before using automatic OR->BETA conversion:
#' \enumerate{
#'   \item Inspect a sample of OR values to verify they look reasonable (typically 0.5-2.0)
#'   \item Check the data source documentation to confirm the column is truly OR
#'   \item For critical analyses, manually prepare clean BETA values from trusted sources
#'   \item When in doubt, provide BETA directly rather than relying on conversion
#' }
#'
#' \strong{Computational Complexity:}
#'
#' O(n) where n = number of OR values. Single-pass vectorized transformation.
#'
#' @examples
#' \dontrun{
#' # Standard conversion
#' .convert_or_to_beta(c(1.5, 2.0, 0.8))
#' # Returns: c(0.405, 0.693, -0.223)
#'
#' # OR = 1 (no effect)
#' .convert_or_to_beta(1.0)
#' # Returns: 0
#'
#' # Invalid OR (with warning)
#' .convert_or_to_beta(-0.5)
#' # Returns: NA (with warning about invalid value)
#'
#' # Suppress warnings
#' .convert_or_to_beta(c(0, -1, Inf), warn = FALSE)
#' # Returns: c(NA, NA, NA) (no warnings)
#' }
#'
#' @keywords internal
#' @noRd
.convert_or_to_beta <- function(or_values, warn = TRUE) {

  # Validate input
  if (!is.numeric(or_values)) {
    or_values <- suppressWarnings(as.numeric(or_values))
  }

  # Initialize output
  beta_values <- rep(NA_real_, length(or_values))

  # Identify valid OR values (positive, finite, not NA)
  valid_idx <- !is.na(or_values) & or_values > 0 & is.finite(or_values)

  # Convert valid values
  beta_values[valid_idx] <- log(or_values[valid_idx])

  # Report issues if requested
  if (warn) {
    n_invalid <- sum(!valid_idx & !is.na(or_values))
    if (n_invalid > 0) {
      warning(sprintf("Converted OR to BETA: %d invalid OR values (<=0 or Inf) set to NA",
                     n_invalid))
    }
  }

  return(beta_values)
}


#' Extract P-value from Text String
#'
#' @description
#' Extracts numeric p-values from text strings with various formats,
#' including scientific notation and special characters. This is a
#' pattern-based extraction with known limitations.
#'
#' \strong{CRITICAL WARNING - User Responsibility:}
#'
#' This function uses heuristic pattern matching to extract p-values from text.
#' While robust for common formats, it is NOT 100 percent accurate and may fail
#' on unusual or non-standard notation. Users are responsible for verifying
#' that extracted p-values are correct for their specific data.
#'
#' @param pvalue_text Character vector of p-value descriptions
#' @param verbose Logical. If TRUE (default), warns about extraction decisions
#'   and potential ambiguities. Set to FALSE to suppress warnings.
#'
#' @return Numeric vector of p-values
#'
#' @details
#' \strong{Pattern Detection and Handling:}
#'
#' The function attempts to parse various text representations of p-values:
#'
#' \strong{Supported Formats:}
#' \itemize{
#'   \item Scientific notation: "5 x 10^-8" -> 5e-8
#'   \item Computer notation: "5E-8" or "5e-8" -> 5e-8
#'   \item Decimal: "0.00000005" -> 5e-8
#'   \item Inequality: "< 0.001" -> 0.001 (takes conservative bound)
#'   \item Non-significant: "NS" or "n.s." -> NA
#'   \item Unicode characters: Handles multiplication signs, superscripts, etc.
#' }
#'
#' \strong{Examples of Handled Patterns:}
#' \itemize{
#'   \item "5 x 10^-8" -> 5e-8 (scientific notation with special chars)
#'   \item "5E-8" -> 5e-8 (standard scientific notation)
#'   \item "< 5e-8" -> 5e-8 (conservative - uses bound value)
#'   \item "NS" -> NA (non-significant, no numeric value)
#' }
#'
#' \strong{IMPORTANT LIMITATIONS:}
#' \itemize{
#'   \item \strong{This is a pattern-based heuristic and may fail on unusual formats}
#'   \item Cannot interpret semantic context beyond pattern matching
#'   \item Inequality symbols ("<", ">") are handled conservatively (uses bound)
#'   \item Complex or non-standard notations may be missed or misinterpreted
#'   \item The function does NOT validate that extracted values make biological sense
#' }
#'
#' \strong{CRITICAL WARNING - User Responsibility:}
#' \itemize{
#'   \item \strong{ALWAYS verify extracted p-values against original data}
#'   \item Check warnings produced during extraction carefully
#'   \item For publication-quality analyses, manually prepare clean numeric P columns
#'   \item When dealing with non-standard formats, pre-process data manually
#'   \item Extraction failures result in NA values - check for unexpected NAs
#' }
#'
#' \strong{Recommendation:} This is a convenience feature for exploratory analysis.
#' For critical analyses or non-standard data formats, prepare a clean numeric P
#' column in your input data before calling \code{prepare_B_training_data()}.
#'
#' \strong{Computational Complexity:}
#'
#' O(n * m) where n = number of text values, m = average text length.
#' Multiple regex operations per value for pattern matching.
#'
#' @examples
#' \dontrun{
#' # Standard scientific notation
#' .extract_pvalue_from_text("5e-8")
#' # Returns: 5e-8
#'
#' # Text with special characters
#' .extract_pvalue_from_text("5 x 10^-8")
#' # Returns: 5e-8
#'
#' # Inequality (conservative bound)
#' .extract_pvalue_from_text("< 0.001")
#' # Returns: 0.001 (warning: using bound value)
#'
#' # Non-significant
#' .extract_pvalue_from_text("NS")
#' # Returns: NA
#'
#' # Suppress warnings
#' .extract_pvalue_from_text("< 0.001", verbose = FALSE)
#' # Returns: 0.001 (no warning)
#' }
#'
#' @keywords internal
#' @noRd
.extract_pvalue_from_text <- function(pvalue_text, verbose = TRUE) {

  # Initialize output
  pvalues <- rep(NA_real_, length(pvalue_text))

  # Skip if already numeric
  if (is.numeric(pvalue_text)) {
    return(pvalue_text)
  }

  # Convert to character if needed
  pvalue_text <- as.character(pvalue_text)

  # Process each value
  for (i in seq_along(pvalue_text)) {
    txt <- pvalue_text[i]

    # Skip NA or empty
    if (is.na(txt) || nchar(trimws(txt)) == 0) {
      next
    }

    # Handle "NS" or "n.s." (non-significant)
    if (grepl("^(NS|n\\.s\\.|not.sig)", txt, ignore.case = TRUE)) {
      pvalues[i] <- NA
      next
    }

    # Handle "< X" format (use X as conservative estimate)
    has_inequality <- grepl("^\\s*[<>]", txt)
    if (has_inequality) {
      if (verbose) {
        warning(sprintf(
          paste0("P-value extraction: Inequality symbol detected in \"%s\".\n",
                 "  Using bound value as conservative estimate. Verify this is appropriate!"),
          txt
        ), call. = FALSE)
      }
      txt <- sub("^\\s*[<>]\\s*", "", txt)
    }

    # Try direct numeric conversion first
    direct <- suppressWarnings(as.numeric(txt))
    if (!is.na(direct)) {
      pvalues[i] <- direct
      next
    }

    # Handle scientific notation with special characters
    # Pattern: "5 x 10^-8" or "5E-8" or special Unicode chars

    # Replace special characters with standard notation
    txt_clean <- txt
    txt_clean <- gsub("\u00D7", "e", txt_clean, fixed = TRUE)  # Multiplication sign
    txt_clean <- gsub("x", "e", txt_clean, ignore.case = TRUE)
    txt_clean <- gsub("\\^", "", txt_clean)
    txt_clean <- gsub("\u207B", "-", txt_clean, fixed = TRUE)  # Superscript minus
    txt_clean <- gsub("\u2078", "8", txt_clean, fixed = TRUE)  # Superscript 8
    txt_clean <- gsub("\u2212", "-", txt_clean, fixed = TRUE)  # Unicode minus sign
    txt_clean <- gsub("\\s+", "", txt_clean)  # Remove spaces

    # Try conversion again
    cleaned <- suppressWarnings(as.numeric(txt_clean))
    if (!is.na(cleaned)) {
      pvalues[i] <- cleaned
      next
    }

    # Try regex extraction for "coefficient e exponent" pattern
    matches <- regexec("([0-9.]+)\\s*[ex]\\s*(-?[0-9]+)", txt_clean, ignore.case = TRUE)
    if (matches[[1]][1] != -1) {
      parts <- regmatches(txt_clean, matches)[[1]]
      if (length(parts) == 3) {
        coefficient <- as.numeric(parts[2])
        exponent <- as.numeric(parts[3])
        if (!is.na(coefficient) && !is.na(exponent)) {
          pvalues[i] <- coefficient * 10^exponent
        }
      }
    }
  }

  return(pvalues)
}


#' Add TRAIT_TYPE Column to Data (Sophisticated Per-Row Auto-Detection)
#'
#' @description
#' Adds a TRAIT_TYPE column to the data with sophisticated per-row detection
#' using multiple evidence sources. This allows mixed-trait datasets where
#' some variants are from binary trait studies and others from continuous
#' trait studies.
#'
#' @param data data.frame with GWAS summary statistics
#' @param trait_type Character string: "binary", "continuous", or NULL (auto-detect)
#'   When specified, serves as a default that can be overridden by row-level indicators
#' @param verbose Verbosity level
#'
#' @return data.frame with TRAIT_TYPE column added (per-row detection)
#'
#' @details
#' **Sophisticated Per-Row Detection Logic (Conservative Approach):**
#'
#' The function uses multiple evidence sources with decreasing confidence:
#'
#' \strong{High Confidence Binary Indicators:}
#' \enumerate{
#'   \item OR column with non-NA, non-1 value (EXCEPT if column name is ambiguous
#'         like "OR.or.BETA", "OR or BETA")
#'   \item Case Size AND Control Size columns both have values
#'   \item Disease name in trait column (e.g., "diabetes", "sclerosis", "cancer")
#'   \item "susceptibility" keyword in any text column
#' }
#'
#' \strong{High Confidence Continuous Indicators:}
#' \enumerate{
#'   \item Quantitative trait name (e.g., "BMI", "height", "bone mineral density")
#'   \item BETA present with no binary indicators
#' }
#'
#' \strong{Conservative Default:}
#' If no high-confidence indicators are found, trait type is set to "unknown"
#' rather than making a potentially incorrect guess.
#'
#' \strong{Important Notes:}
#' \itemize{
#'   \item This detection is robust but NOT error-proof
#'   \item Manual verification is recommended for critical analyses
#'   \item Trait type is less critical when using P-value-based (vs. BETA-based)
#'         methods for B estimation
#'   \item When in doubt, the function defaults to "unknown" to avoid false classifications
#' }
#'
#' When `trait_type` is explicitly specified, it is applied to all rows
#' directly (short-circuit). When NULL, per-row auto-detection is performed.
#'
#' @keywords internal
#' @noRd
.add_trait_type_column <- function(data, trait_type = NULL, verbose = 1) {

  # Short-circuit: if user explicitly specified trait_type, apply to all rows
  if (!is.null(trait_type)) {
    data$TRAIT_TYPE <- rep(trait_type, nrow(data))
    return(data)
  }

  # Initialize TRAIT_TYPE column
  data$TRAIT_TYPE <- rep("unknown", nrow(data))

  # ========== Check for Available Evidence Sources ==========

  # Check for OR column (but verify it's not ambiguous)
  has_or_col <- "OR" %in% names(data)
  or_col_ambiguous <- FALSE

  if (has_or_col) {
    # Check if OR column name suggests ambiguity
    or_related_cols <- grep("OR", names(data), ignore.case = TRUE, value = TRUE)
    for (col_name in or_related_cols) {
      # Ambiguous patterns: "OR.or.BETA", "OR or BETA", "OR_or_BETA"
      if (grepl("OR[._\\s]*(or|OR)[._\\s]*BETA", col_name, ignore.case = TRUE)) {
        or_col_ambiguous <- TRUE
        if (verbose >= 2) {
          message(sprintf("  Warning: Ambiguous OR column name detected ('%s'). ",
                         "Will not use for trait type detection.", col_name))
        }
        break
      }
    }
  }

  has_beta_col <- "BETA" %in% names(data)

  # Check for case/control columns
  case_control_cols <- .detect_case_control_columns(data)
  has_case_control <- !is.null(case_control_cols$case_col) &&
                      !is.null(case_control_cols$control_col)

  # Check for trait name column
  trait_col <- .detect_trait_column(data)

  # ========== Per-Row Auto-Detection ==========

  for (i in seq_len(nrow(data))) {
    trait_detected <- "unknown"
    confidence <- 0  # Track confidence level

    # ===== Binary Trait Indicators =====

    # Indicator 1: OR value present (non-ambiguous column)
    if (has_or_col && !or_col_ambiguous) {
      if (!is.na(data$OR[i]) && data$OR[i] != 1) {
        trait_detected <- "binary"
        confidence <- max(confidence, 3)  # High confidence
      }
    }

    # Indicator 2: Case and Control sizes present
    if (has_case_control && confidence < 3) {
      case_val <- data[[case_control_cols$case_col]][i]
      control_val <- data[[case_control_cols$control_col]][i]

      if (!is.na(case_val) && !is.na(control_val) &&
          case_val > 0 && control_val > 0) {
        trait_detected <- "binary"
        confidence <- max(confidence, 3)  # High confidence
      }
    }

    # Indicator 3: Disease name in trait column
    if (!is.null(trait_col) && confidence < 3) {
      trait_name <- as.character(data[[trait_col]][i])
      if (!is.na(trait_name) && nchar(trait_name) > 0) {
        disease_detected <- .is_disease_name(trait_name)
        if (disease_detected$is_disease) {
          trait_detected <- "binary"
          confidence <- max(confidence, 2)  # Medium-high confidence
        } else if (disease_detected$is_continuous) {
          trait_detected <- "continuous"
          confidence <- max(confidence, 2)
        }
      }
    }

    # Indicator 4: "susceptibility" keyword in text columns
    if (confidence < 3) {
      susceptibility_detected <- .check_for_susceptibility(data, i)
      if (susceptibility_detected) {
        trait_detected <- "binary"
        confidence <- max(confidence, 2)
      }
    }

    # ===== Continuous Trait Indicators (if no binary evidence) =====

    if (confidence == 0) {
      # Indicator 5: BETA present with no binary indicators
      if (has_beta_col && !is.na(data$BETA[i])) {
        # Only use BETA as continuous indicator if no OR column exists
        # or OR is explicitly NA
        if (!has_or_col || or_col_ambiguous || is.na(data$OR[i])) {
          trait_detected <- "continuous"
          confidence <- 1  # Lower confidence
        }
      }
    }

    # Assign detected trait type (or keep "unknown")
    data$TRAIT_TYPE[i] <- trait_detected
  }

  # ========== Apply User-Specified Default to Unknown Rows ==========

  if (!is.null(trait_type)) {
    unknown_idx <- data$TRAIT_TYPE == "unknown"
    if (any(unknown_idx)) {
      data$TRAIT_TYPE[unknown_idx] <- trait_type
      if (verbose >= 2) {
        message(sprintf("  Applied default trait_type='%s' to %d variants with unknown type",
                       trait_type, sum(unknown_idx)))
      }
    }
  }

  # Report results will be done after filtering in main function
  return(data)
}


#################### TRAIT TYPE DETECTION HELPERS ####################

#' Detect Case and Control Size Columns
#'
#' @description
#' Identifies columns containing case and control sample sizes.
#'
#' @param data data.frame to search
#'
#' @return List with case_col and control_col (NULL if not found)
#'
#' @keywords internal
#' @noRd
.detect_case_control_columns <- function(data) {

  # Common patterns for case size columns
  # Use explicit space, underscore, dot variations
  case_patterns <- c(
    "^Case Size$",            # "Case Size" (exact with space)
    "^Case_Size$",            # "Case_Size"
    "^Case\\.Size$",          # "Case.Size"
    "^CaseSize$",             # "CaseSize"
    "^N cases?$",             # "N cases", "N case"
    "^N_cases?$",             # "N_cases"
    "^Ncases?$",              # "Ncases"
    "^CASE N$",               # "CASE N"
    "^CASE_N$",               # "CASE_N"
    "^Cases$"
  )

  # Common patterns for control size columns
  control_patterns <- c(
    "^Control Size$",         # "Control Size" (exact with space)
    "^Control_Size$",         # "Control_Size"
    "^Control\\.Size$",       # "Control.Size"
    "^ControlSize$",          # "ControlSize"
    "^N controls?$",          # "N controls", "N control"
    "^N_controls?$",          # "N_controls"
    "^Ncontrols?$",           # "Ncontrols"
    "^CONTROL N$",            # "CONTROL N"
    "^CONTROL_N$",            # "CONTROL_N"
    "^Controls$"
  )

  col_names <- names(data)
  case_col <- NULL
  control_col <- NULL

  # Search for case column
  for (pattern in case_patterns) {
    matches <- grep(pattern, col_names, ignore.case = TRUE, value = TRUE)
    if (length(matches) > 0) {
      case_col <- matches[1]
      break
    }
  }

  # Search for control column
  for (pattern in control_patterns) {
    matches <- grep(pattern, col_names, ignore.case = TRUE, value = TRUE)
    if (length(matches) > 0) {
      control_col <- matches[1]
      break
    }
  }

  return(list(case_col = case_col, control_col = control_col))
}


#' Detect Trait Name Column
#'
#' @description
#' Identifies the column containing trait or disease names.
#'
#' @param data data.frame to search
#'
#' @return Column name (NULL if not found)
#'
#' @keywords internal
#' @noRd
.detect_trait_column <- function(data) {

  # Common trait column name patterns (in order of preference)
  trait_patterns <- c(
    "^TRAIT$",
    "^Reported[\\s_.]?trait$",
    "^DISEASE[\\s_.]?TRAIT$",
    "^MAPPED[\\s_.]?TRAIT$",
    "^Trait\\(s\\)$",
    "^PHENOTYPE$"
  )

  col_names <- names(data)

  for (pattern in trait_patterns) {
    matches <- grep(pattern, col_names, ignore.case = TRUE, value = TRUE)
    if (length(matches) > 0) {
      return(matches[1])
    }
  }

  return(NULL)
}


#' Check if Trait Name is a Disease or Continuous Trait
#'
#' @description
#' Uses pattern matching and a knowledge base to classify trait names.
#'
#' @param trait_name Character string with trait name
#'
#' @return List with is_disease (logical) and is_continuous (logical)
#'
#' @keywords internal
#' @noRd
.is_disease_name <- function(trait_name) {

  trait_lower <- tolower(trait_name)

  # ===== Binary Trait / Disease Patterns =====

  # Disease keywords
  disease_keywords <- c(
    "sclerosis",      # ALS, multiple sclerosis, etc.
    "diabetes",       # Type 1, Type 2 diabetes
    "cancer",         # Various cancers
    "carcinoma",
    "tumor",
    "lymphoma",
    "leukemia",
    "disease",        # General disease
    "disorder",
    "syndrome",
    "osteoporosis",
    "arthritis",
    "alzheimer",
    "parkinson",
    "asthma",
    "autism",
    "schizophrenia",
    "depression",
    "hypertension",
    "coronary",       # CAD
    "myocardial",     # MI
    "infarction",
    "stroke",
    "inflammatory",   # IBD, etc.
    "infection",
    "obesity",        # Often treated as binary in GWAS
    "failure",        # Heart failure, kidney failure
    "deficiency",
    "allergic",
    "allergy",
    "fracture"        # Bone fracture
  )

  for (keyword in disease_keywords) {
    if (grepl(keyword, trait_lower)) {
      return(list(is_disease = TRUE, is_continuous = FALSE))
    }
  }

  # ===== Continuous Trait Patterns =====

  # Quantitative trait keywords (exact or partial matches)
  continuous_patterns <- c(
    "bone mineral density",
    "^bmd$",          # Exact BMD
    "\\bbmd\\b",      # BMD as word
    "^bmi$",          # Exact BMI
    "body mass index",
    "height",
    "weight",
    "^ldl$",          # Cholesterol
    "^hdl$",
    "cholesterol",
    "triglyceride",
    "glucose",
    "insulin",
    "^hba1c$",
    "blood pressure",  # Continuous measurement
    "systolic",
    "diastolic",
    "waist",
    "hip ratio",
    "^fev1$",         # Lung function
    "lung function",
    "^egfr$",         # Kidney function
    "creatinine",
    "uric acid",
    "vitamin d",
    "testosterone",
    "estrogen",
    "cortisol",
    "^iq$",           # Cognitive traits
    "cognitive",
    "intelligence",
    "years of education"
  )

  for (pattern in continuous_patterns) {
    if (grepl(pattern, trait_lower)) {
      return(list(is_disease = FALSE, is_continuous = TRUE))
    }
  }

  # ===== Ambiguous - Cannot Classify =====
  return(list(is_disease = FALSE, is_continuous = FALSE))
}


#' Check for "Susceptibility" Keyword in Row
#'
#' @description
#' Searches text columns in a row for the "susceptibility" keyword,
#' which indicates a binary disease trait.
#'
#' @param data data.frame
#' @param row_idx Integer row index
#'
#' @return Logical indicating if "susceptibility" was found
#'
#' @keywords internal
#' @noRd
.check_for_susceptibility <- function(data, row_idx) {

  # Columns to search (typically annotation or description columns)
  search_patterns <- c(
    "annotation",
    "P.value.annotation",
    "P-value annotation",
    "description",
    "note",
    "comment",
    "context"
  )

  # Get character/text columns that match patterns
  col_names <- names(data)
  search_cols <- c()

  for (pattern in search_patterns) {
    matches <- grep(pattern, col_names, ignore.case = TRUE, value = TRUE)
    search_cols <- c(search_cols, matches)
  }

  # Remove duplicates
  search_cols <- unique(search_cols)

  # Search each column for "susceptibility"
  for (col in search_cols) {
    if (col %in% names(data)) {
      value <- as.character(data[[col]][row_idx])
      if (!is.na(value) && grepl("susceptibility", value, ignore.case = TRUE)) {
        return(TRUE)
      }
    }
  }

  return(FALSE)
}


#################### FILTERING HELPERS ####################


#' Apply General Inclusion and Exclusion Filters
#'
#' @description
#' Applies flexible filtering on any columns in the data. Supports both
#' inclusion and exclusion filters with pattern matching or exact matching.
#'
#' @param data data.frame to filter
#' @param include Named list of inclusion filters (column_name = values_to_include)
#' @param exclude Named list of exclusion filters (column_name = values_to_exclude)
#' @param pattern_matching Logical. If TRUE, uses grepl for pattern matching on character columns
#' @param verbose Verbosity level
#'
#' @return Filtered data.frame
#'
#' @details
#' Inclusion filters take precedence. For each column in include, only rows
#' matching at least one value are kept. For columns in exclude (not in include),
#' rows matching any value are removed.
#'
#' Pattern matching (when enabled):
#' - For character columns: uses grepl (case-insensitive partial matching)
#' - For numeric columns: uses exact matching (==, %in%)
#'
#' @keywords internal
#' @noRd
.apply_general_filters <- function(data, include = NULL, exclude = NULL,
                                   pattern_matching = TRUE, verbose = 1) {

  n_before <- nrow(data)

  # ========== Apply Inclusion Filters (takes precedence) ==========

  if (!is.null(include) && length(include) > 0) {
    for (col_name in names(include)) {

      # Check if column exists
      if (!col_name %in% names(data)) {
        if (verbose >= 1) {
          warning(sprintf("Inclusion filter column '%s' not found in data. Skipping.", col_name))
        }
        next
      }

      filter_values <- include[[col_name]]
      col_data <- data[[col_name]]

      # Create logical vector for rows to keep
      keep_idx <- rep(FALSE, nrow(data))

      # Apply filter based on column type and pattern_matching setting
      if (is.character(col_data) && pattern_matching) {
        # Pattern matching for character columns (substring matching)
        # Use fixed=TRUE to treat special regex characters literally
        # Convert to lowercase for case-insensitive matching since ignore.case doesn't work with fixed
        for (val in filter_values) {
          keep_idx <- keep_idx | grepl(tolower(val), tolower(col_data), fixed = TRUE)
        }
      } else {
        # Exact matching (works for both character and numeric)
        keep_idx <- col_data %in% filter_values
      }

      # Handle NAs (exclude them by default)
      keep_idx[is.na(keep_idx)] <- FALSE

      # Apply filter
      n_before_col <- nrow(data)
      data <- data[keep_idx, , drop = FALSE]
      n_removed_col <- n_before_col - nrow(data)

      if (verbose >= 1 && n_removed_col > 0) {
        message(sprintf("  [INCLUDE] Kept %d variants matching '%s' in %s (removed %d)",
                       nrow(data), paste(filter_values, collapse = "' OR '"),
                       col_name, n_removed_col))
      }
    }
  }

  # ========== Apply Exclusion Filters ==========

  if (!is.null(exclude) && length(exclude) > 0) {
    # Get columns to exclude (skip those already handled by include)
    include_cols <- if (!is.null(include)) names(include) else character(0)
    exclude_cols <- setdiff(names(exclude), include_cols)

    for (col_name in exclude_cols) {

      # Check if column exists
      if (!col_name %in% names(data)) {
        if (verbose >= 1) {
          warning(sprintf("Exclusion filter column '%s' not found in data. Skipping.", col_name))
        }
        next
      }

      filter_values <- exclude[[col_name]]
      col_data <- data[[col_name]]

      # Create logical vector for rows to remove
      remove_idx <- rep(FALSE, nrow(data))

      # Apply filter based on column type and pattern_matching setting
      if (is.character(col_data) && pattern_matching) {
        # Pattern matching for character columns (substring matching)
        # Use fixed=TRUE to treat special regex characters literally
        # Convert to lowercase for case-insensitive matching since ignore.case doesn't work with fixed
        for (val in filter_values) {
          remove_idx <- remove_idx | grepl(tolower(val), tolower(col_data), fixed = TRUE)
        }
      } else {
        # Exact matching (works for both character and numeric)
        remove_idx <- col_data %in% filter_values
      }

      # Handle NAs (keep them by default - only remove explicit matches)
      remove_idx[is.na(remove_idx)] <- FALSE

      # Apply filter
      n_before_col <- nrow(data)
      data <- data[!remove_idx, , drop = FALSE]
      n_removed_col <- n_before_col - nrow(data)

      if (verbose >= 1 && n_removed_col > 0) {
        message(sprintf("  [EXCLUDE] Removed %d variants matching '%s' in %s",
                       n_removed_col, paste(filter_values, collapse = "' OR '"),
                       col_name))
      }
    }
  }

  # ========== Summary ==========

  n_after <- nrow(data)
  n_total_removed <- n_before - n_after

  if (verbose >= 1 && n_total_removed > 0) {
    message(sprintf("  Total filtering: %d -> %d variants (%d removed)",
                   n_before, n_after, n_total_removed))
  }

  return(data)
}




#' Remove Duplicates Based on Genomic Position
#'
#' @description
#' Removes duplicate variants based on chromosome and position,
#' keeping the variant with the highest priority value (typically largest N).
#'
#' @param data data.frame with CHR and POS columns
#' @param priority Character string specifying priority column (default: "N")
#' @param verbose Verbosity level
#'
#' @return data.frame with duplicates removed
#'
#' @details
#' Duplicates are identified by identical CHR and POS values.
#' When duplicates exist, the variant with the largest value in the priority
#' column is kept. Also removes duplicates by rsID if present.
#'
#' @keywords internal
#' @noRd
.remove_duplicates_by_position <- function(data, priority = "N", verbose = 1) {

  n_before <- nrow(data)

  # Step 1: Remove duplicates by rsID (if column exists)
  if ("rsID" %in% names(data)) {
    # Sort by rsID and priority (descending)
    if (priority %in% names(data)) {
      data <- data[order(data$rsID, -data[[priority]], na.last = TRUE), ]
    } else {
      data <- data[order(data$rsID), ]
    }

    # Keep first occurrence (highest priority) of each rsID
    data <- data[!duplicated(data$rsID), ]

    n_after_rsid <- nrow(data)
    n_removed_rsid <- n_before - n_after_rsid
  } else {
    n_removed_rsid <- 0
    n_after_rsid <- n_before
  }

  # Step 2: Remove duplicates by CHR + POS (if both columns exist)
  if ("CHR" %in% names(data) && "POS" %in% names(data)) {

    # Create position identifier (chr:pos)
    data$..temp_chr_pos <- paste(data$CHR, data$POS, sep = ":")

    # Sort by position and priority (descending)
    if (priority %in% names(data)) {
      data <- data[order(data$..temp_chr_pos, -data[[priority]], na.last = TRUE), ]
    } else {
      data <- data[order(data$..temp_chr_pos), ]
    }

    # Keep first occurrence (highest priority) of each position
    dup_pos <- duplicated(data$..temp_chr_pos)
    data <- data[!dup_pos, ]

    # Remove temporary column
    data$..temp_chr_pos <- NULL

    n_after_pos <- nrow(data)
    n_removed_pos <- n_after_rsid - n_after_pos
  } else {
    n_removed_pos <- 0
    n_after_pos <- n_after_rsid
  }

  # Report results
  if (verbose >= 1) {
    if (n_removed_rsid > 0 || n_removed_pos > 0) {
      message(sprintf("  Removed duplicates: %d by rsID, %d by CHR:POS (kept highest %s)",
                     n_removed_rsid, n_removed_pos, priority))
    }
  }

  return(data)
}


#################### COLUMN MANAGEMENT HELPERS ####################

#' Reorder Columns to Standard Format
#'
#' @description
#' Reorders columns to a standard format for consistency and readability.
#'
#' @param data data.frame to reorder
#'
#' @return data.frame with reordered columns
#'
#' @details
#' Standard order:
#' 1. Primary: rsID, CHR, POS, MAF, P, P_mlog10, N, BETA, TRAIT_TYPE
#' 2. Study metadata: TRAIT, STUDY, PMID, FIRST_AUTHOR
#' 3. Everything else (alphabetically sorted)
#'
#' @keywords internal
#' @noRd
.reorder_columns_standard <- function(data) {

  # Define standard order
  primary_cols <- c("rsID", "CHR", "POS", "MAF", "P", "P_mlog10", "N", "BETA", "TRAIT_TYPE")
  study_cols <- c("TRAIT", "STUDY", "PMID", "FIRST_AUTHOR")

  # Get current column names
  current_cols <- names(data)

  # Identify which standard columns are present
  present_primary <- intersect(primary_cols, current_cols)
  present_study <- intersect(study_cols, current_cols)

  # Identify remaining columns (not in standard lists)
  remaining_cols <- setdiff(current_cols, c(primary_cols, study_cols))
  remaining_cols <- sort(remaining_cols)  # Alphabetically

  # Create final column order
  final_order <- c(present_primary, present_study, remaining_cols)

  # Reorder
  data <- data[, final_order, drop = FALSE]

  return(data)
}


#' Ensure CHR Column is Numeric and Filter to Autosomes
#'
#' @description
#' Converts CHR column to numeric and removes non-autosomal chromosomes.
#' Retains only variants on chromosomes 1-22 (autosomes).
#'
#' @param data data.frame with CHR column
#' @param verbose Verbosity level
#'
#' @return data.frame with numeric CHR column (1-22 only), non-autosomal variants removed
#'
#' @details
#' \strong{Processing Steps:}
#' \enumerate{
#'   \item Convert CHR to numeric
#'   \item Remove variants on sex chromosomes (X, Y)
#'   \item Remove variants on mitochondrial chromosome (MT)
#'   \item Remove variants with invalid CHR or outside 1-22 range
#'   \item Retain only autosomal variants (CHR 1-22)
#' }
#'
#' Removed: X, Y, MT, and any invalid/unknown chromosomes
#' Retained: 1-22 (autosomes only)
#'
#' \strong{IMPORTANT LIMITATIONS:}
#' \itemize{
#'   \item \strong{Non-standard chromosome naming may not be recognized}
#'   \item Assumes human chromosomes (22 autosomes)
#'   \item Cannot handle species-specific chromosome conventions
#'   \item Sex chromosomes are completely removed (including pseudoautosomal regions)
#'   \item Unusual naming conventions may cause unexpected filtering
#' }
#'
#' \strong{CRITICAL WARNING - Data Loss:}
#' \itemize{
#'   \item \strong{All non-autosomal variants are permanently removed}
#'   \item X chromosome variants removed (important for some phenotypes)
#'   \item Y chromosome variants removed (may be relevant for some studies)
#'   \item Mitochondrial variants removed (if present)
#'   \item Cannot be undone - ensure this filtering is appropriate for your analysis
#' }
#'
#' \strong{Common Issues:}
#' \itemize{
#'   \item Non-standard naming not recognized: "chr1" vs "1" vs "Chr1"
#'   \item Pseudoautosomal regions on X/Y incorrectly removed
#'   \item Species with different numbers of chromosomes not supported
#'   \item Mixed naming conventions cause inconsistent filtering
#'   \item Invalid CHR values silently removed (check verbose output)
#' }
#'
#' \strong{Chromosome Recognition Patterns:}
#' \itemize{
#'   \item X chromosome: "X", "CHRX", "CHR23", "23"
#'   \item Y chromosome: "Y", "CHRY", "CHR24", "24"
#'   \item Mitochondrial: "MT", "M", "CHRMT", "CHRM", "CHR25", "25"
#'   \item Autosomes: Numeric values 1-22
#' }
#'
#' \strong{Recommendation:}
#' \itemize{
#'   \item This is the default and recommended option for most GLOW analyses
#'   \item Verify chromosome counts in verbose output to ensure expected filtering
#'   \item For X-linked analyses, use filter_autosomes=FALSE and handle separately
#'   \item Document filtering decisions in analysis records
#'   \item For non-human data, manually prepare CHR column before processing
#' }
#'
#' @keywords internal
#' @noRd
.ensure_chr_numeric <- function(data, verbose = 1) {

  if (!"CHR" %in% names(data)) {
    return(data)
  }

  n_before <- nrow(data)

  # Convert to character first
  chr_orig <- as.character(data$CHR)
  chr_numeric <- rep(NA_real_, length(chr_orig))

  # Identify sex and mitochondrial chromosomes (to be removed)
  is_x <- toupper(chr_orig) %in% c("X", "CHR23", "23", "CHRX")
  is_y <- toupper(chr_orig) %in% c("Y", "CHR24", "24", "CHRY")
  is_mt <- toupper(chr_orig) %in% c("MT", "M", "CHR25", "25", "CHRM", "CHRMT")

  # Convert numeric strings
  chr_numeric <- suppressWarnings(as.numeric(chr_orig))

  # Count non-autosomal chromosomes
  n_sex_chr <- sum(is_x | is_y, na.rm = TRUE)
  n_mt <- sum(is_mt, na.rm = TRUE)
  n_invalid <- sum(is.na(chr_numeric) & !is_x & !is_y & !is_mt)
  n_out_of_range <- sum(!is.na(chr_numeric) & (chr_numeric < 1 | chr_numeric > 22), na.rm = TRUE)

  # Update CHR column with numeric values
  data$CHR <- chr_numeric

  # Filter to keep only autosomal chromosomes (1-22)
  # Remove: sex chromosomes (X, Y), mitochondrial (MT), invalid, and out-of-range
  keep_idx <- !is.na(data$CHR) & data$CHR >= 1 & data$CHR <= 22
  data <- data[keep_idx, , drop = FALSE]

  n_after <- nrow(data)
  n_removed <- n_before - n_after

  # Report what was removed
  if (verbose >= 1 && n_removed > 0) {
    message(sprintf("  Removed %d non-autosomal variants:", n_removed))
    if (n_sex_chr > 0) message(sprintf("    - %d on sex chromosomes (X/Y)", n_sex_chr))
    if (n_mt > 0) message(sprintf("    - %d on mitochondrial chromosome (MT)", n_mt))
    if (n_invalid > 0) message(sprintf("    - %d with invalid CHR", n_invalid))
    if (n_out_of_range > 0) message(sprintf("    - %d with CHR outside 1-22", n_out_of_range))
    message(sprintf("  Retained %d autosomal variants (CHR 1-22)", n_after))
  }

  return(data)
}
