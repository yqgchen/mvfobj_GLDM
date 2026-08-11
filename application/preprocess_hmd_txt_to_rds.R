# =============================================================================
# File: preprocess_hmd_txt_to_rds.R
# Purpose: Parse raw HMD period life table txt files and save one combined
#          (female + male) rds file per country for downstream analysis.
# Dependencies: base R only (read.delim, saveRDS)
# =============================================================================

# -----------------------------------------------------------------------------
# HMD txt-to-rds conversion function
# -----------------------------------------------------------------------------

#' Convert HMD life table txt files to country-level rds files
#'
#' @description
#' Reads HMD period life table txt files (1x1 format, i.e., 1-year age and
#' 1-year period resolution) for females and males separately, extracts the
#' \code{Year}, \code{Age}, and \code{dx} (death count) columns, standardises
#' the open age group label \code{"110+"} to \code{"110"}, appends a
#' \code{sex} factor column, row-binds the two sexes, and saves the result as
#' an rds file named \code{<lowercase_code>_dx.rds} in \code{out_dir}.
#'
#' @param base_dir  Character; root directory that must contain subdirectories
#'                  \code{lt_female/fltper_1x1/} and \code{lt_male/mltper_1x1/}.
#' @param out_dir   Character; directory where output rds files are written.
#'                  Created recursively if it does not exist.
#' @param countries Character vector of HMD country codes to process (e.g.,
#'                  \code{"USA"}, \code{"FRATNP"}). If \code{NULL}, all
#'                  countries found in \code{lt_female/fltper_1x1/} are used.
#'
#' @return Invisibly returns a character vector of the country codes that were
#'         successfully processed.
preprocess_hmd_txt_to_rds <- function(
    base_dir,
    out_dir,
    countries = NULL
) {
  # Construct paths to the HMD per-sex subdirectories
  female_dir <- file.path(base_dir, "lt_female", "fltper_1x1")
  male_dir   <- file.path(base_dir, "lt_male", "mltper_1x1")

  # Validate input directories exist before proceeding
  if (!dir.exists(female_dir)) stop("Female life table directory not found: ", female_dir)
  if (!dir.exists(male_dir))   stop("Male life table directory not found: ", male_dir)
  # Create output directory if needed
  if (!dir.exists(out_dir))    dir.create(out_dir, recursive = TRUE)

  # Auto-detect countries from filenames if the user did not supply a list
  if (is.null(countries)) {
    # Strip the standard HMD suffix to recover the country code
    countries <- sub("\\.fltper_1x1\\.txt$", "",
                     list.files(female_dir, pattern = "\\.fltper_1x1\\.txt$"))
  }

  # Inner function that processes a single country code
  process_one <- function(code) {
    # Build the expected file paths using the HMD naming convention
    f_file <- file.path(female_dir, paste0(code, ".fltper_1x1.txt"))
    m_file <- file.path(male_dir,   paste0(code, ".mltper_1x1.txt"))

    # Warn and skip rather than error if a file is missing
    if (!file.exists(f_file)) { warning("Female file not found: ", f_file); return(NULL) }
    if (!file.exists(m_file)) { warning("Male file not found: ", m_file);   return(NULL) }

    # HMD txt files have one title line before the column header; skip = 1 drops it
    f.data <- read.delim(f_file, skip = 1, header = TRUE, sep = "")
    m.data <- read.delim(m_file, skip = 1, header = TRUE, sep = "")

    # --- Female: extract Year, Age, dx and standardise ---
    female.dxs <- f.data[, c("Year", "Age", "dx")]
    # HMD uses "110+" for the open-ended top age group; convert to "110"
    # so it can be coerced to integer without producing NA
    female.dxs[female.dxs$Age == "110+", "Age"] <- "110"
    female.dxs$age <- as.integer(female.dxs$Age)
    colnames(female.dxs)[1] <- "year"
    # Encode sex as a factor with a fixed level order (female first)
    female.dxs$sex <- factor("female", levels = c("female", "male"))
    female.dxs <- female.dxs[, c("year", "sex", "age", "dx")]

    # --- Male: same processing steps as female ---
    male.dxs <- m.data[, c("Year", "Age", "dx")]
    male.dxs[male.dxs$Age == "110+", "Age"] <- "110"
    male.dxs$age <- as.integer(male.dxs$Age)
    colnames(male.dxs)[1] <- "year"
    male.dxs$sex <- factor("male", levels = c("female", "male"))
    male.dxs <- male.dxs[, c("year", "sex", "age", "dx")]

    # Combine female and male rows and save as rds
    dx <- rbind(female.dxs, male.dxs)
    out_file <- file.path(out_dir, paste0(tolower(code), "_dx.rds"))
    saveRDS(dx, out_file)
    cat(sprintf("  [OK] %s -> %s (%d rows)\n", code, basename(out_file), nrow(dx)))
    return(code)
  }

  # Process each country in sequence; NULL entries indicate failures
  cat(sprintf("Processing %d countries...\n", length(countries)))
  results <- lapply(countries, process_one)
  processed <- unlist(results[!sapply(results, is.null)])
  cat(sprintf("Done. %d/%d countries processed.\n", length(processed), length(countries)))

  invisible(processed)
}

# -----------------------------------------------------------------------------
# Standalone execution block
# -----------------------------------------------------------------------------

# ====================================================================
# Run: convert 35 countries used in the paper
# ====================================================================
# This block executes only when the script is run directly (not sourced),
# detected by checking that the call stack depth is 0
if (sys.nframe() == 0) {
  paper_countries <- c(
    "AUS", "AUT", "BEL", "BGR", "BLR", "CAN", "CHE", "CZE",
    "DEUTE", "DEUTW", "DNK", "ESP", "EST", "FIN", "FRATNP",
    "GBR_NP", "GRC", "HUN", "IRL", "ISL", "ISR", "ITA", "JPN",
    "LTU", "LUX", "LVA", "NLD", "NOR", "NZL_NP", "POL", "PRT",
    "SVK", "SVN", "SWE", "USA"
  )

  preprocess_hmd_txt_to_rds(
    base_dir  = ".",
    out_dir   = "data",
    countries = paper_countries
  )
}
