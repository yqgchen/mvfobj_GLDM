# =============================================================================
# File: analyze_data.R
# Purpose: Run the full GLDM mortality data analysis pipeline—from raw HMD
#          text files through preprocessing, model fitting, and saving results.
# Dependencies: preprocess_hmd_txt_to_rds.R, fctns_preprocess.R,
#               fctns_estimation.R; libraries purrr, frechet, fdadensity,
#               fdapace, pracma, plyr, dplyr, tidyr, reshape2, OPW, minqa
# =============================================================================

# -----------------------------------------------------------------------------
# Package loading
# -----------------------------------------------------------------------------

library(purrr)
library(frechet)
library(fdadensity)
library(fdapace)
library(pracma)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
if (!requireNamespace("OPW", quietly = TRUE)) {
  # Install the OPW package from GitHub if not already available
  devtools::install_github("yqgchen/OPW", ref = "main", upgrade = "never")
}
library(OPW)

library(minqa)

# -----------------------------------------------------------------------------
# Path setup
# -----------------------------------------------------------------------------

# Set working directory to the location of this script (RStudio only)
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  rstudioapi::getSourceEditorContext()$path %>% dirname %>% setwd
}

DATA_DIR <- "data"      # directory containing raw HMD txt life table files
RDS_DIR  <- "data_rds"  # directory for intermediate per-country rds files
OUT_DIR <- "results"    # directory for final output rds files

# Create the results directory if it does not yet exist
if ( !dir.exists(OUT_DIR) ) {
  dir.create( OUT_DIR )
}

# -----------------------------------------------------------------------------
# Global parameters
# -----------------------------------------------------------------------------

START_YR   <- 1989       # first calendar year included in the analysis
END_YR     <- 2016       # last calendar year included in the analysis
MAX_AGE    <- 100        # maximum age (exclusive upper bound for life-table bins)
N_KNOTS    <- 4          # number of interior knots used in pairwise synchronization
REF_SOURCE <- "pooled"   # reference point source: pooled across both sexes

# -----------------------------------------------------------------------------
# Source helper scripts
# -----------------------------------------------------------------------------

source("preprocess_hmd_txt_to_rds.R")  # defines preprocess_hmd_txt_to_rds()
source("fctns_preprocess.R")            # defines data_gen()
source("fctns_estimation.R")            # defines ldm_distn(), get_pxcalign(), get_sxcalign(), etc.

# -----------------------------------------------------------------------------
# Step 1: Convert HMD txt files to country-level rds (run once)
# -----------------------------------------------------------------------------

# Use the USA file as a sentinel to decide whether conversion has already run
rds_check <- file.path(RDS_DIR, "usa_dx.rds")
if ( !file.exists(rds_check) ) {
  cat("Converting HMD txt files to rds...\n")
  # Process the 35 countries used in the paper; DEUTE and DEUTW (East/West
  # Germany) are read separately and merged inside data_gen()
  preprocess_hmd_txt_to_rds(
    base_dir  = DATA_DIR,
    out_dir   = RDS_DIR,
    countries = c(
      "AUS", "AUT", "BEL", "BGR", "BLR", "CAN", "CHE", "CZE",
      "DEUTE", "DEUTW", "DNK", "ESP", "EST", "FIN", "FRATNP",
      "GBR_NP", "GRC", "HUN", "IRL", "ISL", "ISR", "ITA", "JPN",
      "LTU", "LUX", "LVA", "NLD", "NOR", "NZL_NP", "POL", "PRT",
      "SVK", "SVN", "SWE", "USA"
    )
  )
}

# -----------------------------------------------------------------------------
# Step 2: Preprocess data — build trajectories of quantile functions
# -----------------------------------------------------------------------------

data_file <- file.path(OUT_DIR, "mort_preprocessed.rds")
if ( file.exists(data_file) ) {
  # Load cached preprocessed data to avoid rerunning the expensive LFR step
  dat <- readRDS(data_file)
} else {
  # data_gen() converts life-table death counts into histogram, density, and
  # quantile-function representations, then smooths via local Fréchet regression
  dat <- data_gen(start_yr = START_YR, end_yr = END_YR, data_dir = RDS_DIR, max_age = MAX_AGE)
  # Attach human-readable country names (keyed by the three-letter abbreviations
  # in dat$subj_labs) for use in plots
  dat$subj_names <- c(
    "Australia","Austria","Belgium","Bulgaria","Belarus","Canada",
    "Switzerland","Czech Republic","Germany","Denmark","Spain","Estonia","Finland",
    "France","United Kingdom","Greece","Hungary","Ireland","Iceland",
    "Israel","Italy","Japan","Lithuania","Luxembourg","Latvia","Netherlands",
    "Norway","New Zealand","Poland","Portugal","Slovakia","Slovenia",
    "Sweden","United States"
  ) %>% set_names(dat$subj_labs)
  # Use atomic save to prevent partial writes from corrupting the cache file
  safe_save_rds( obj = dat, path = data_file )
}

# -----------------------------------------------------------------------------
# Step 3: GLDM estimation
# -----------------------------------------------------------------------------

# Use the LFR-smoothed quantile functions (LqfLFR) as model input
Lqf    <- dat$LqfLFR
qSup   <- dat$qSup      # quantile support grid on [0, 1] (length 201)
obsGrid <- dat$obsGrid  # calendar years 1989–2016

# Density support grid for converting estimated quantile functions back to
# densities; 201 equally spaced points from 0 to MAX_AGE
dSup <- seq(0, MAX_AGE, length.out = 201)

# Options passed to frechet::DenFMean() and frechet::LocDenReg() throughout the
# estimation: Epanechnikov kernel, cross-validated bandwidth, bounded support
optns <- list(
  qSup = qSup,
  lower = 0, upper = MAX_AGE,
  dSup = dSup,
  kernel = "epan", bwReg = "CV"
)

## fit the LDM ----
# ldm_distn() estimates omega_0, amplitude factors A, subject-level deformation
# maps H_i, component tempos eta_j, the latent template tau, and component-level
# deformation maps Psi_j in a single call.  ref_source = "pooled" computes the
# reference point from the pooled distribution across both sexes and all years.
ldmfit <- ldm_distn(
  obsGrid = obsGrid,
  Lqf = Lqf,
  qSup = qSup,
  ref_source = REF_SOURCE,
  nknots = N_KNOTS,
  optns = optns
)

## compute population-level cross-component alignment map from female to male ----
# pxcalign_f2m(t) = Psi_female^{-1}(Psi_male(t)), capturing the average
# timing lead/lag of males relative to females at the population level
ldmfit$pxcalign_f2m <- with(
  ldmfit,
  get_pxcalign(
    j1 = 'female', j2 = 'male',
    obsGrid = obsGrid,
    Psi = res_Psi$H, PsiInv = res_Psi$Hinv
  )
)

## compute subject-level cross-component alignment map from female to male ----
# sxcalign_f2m[i,] = G_female_i^{-1} o G_male_i(t), giving the country-specific
# female-to-male timing alignment for subject i.  Uses the final global
# deformation functions res_G (not the preliminary res_G_pre).
ldmfit$sxcalign_f2m <- with(
  ldmfit,
  get_sxcalign(
    j1 = 'female', j2 = 'male',
    obsGrid = obsGrid,
    G = res_G$H, Ginv = res_G$Hinv
  )
)

# Save the full fitted model object atomically
safe_save_rds(ldmfit, file.path(OUT_DIR, "mort_ldmfit.rds") )

# -----------------------------------------------------------------------------
# Structure of ldmfit (saved to mort_ldmfit.rds)
# -----------------------------------------------------------------------------
# ldmfit is a list of the following fields:
# res_omega0: A list holding the reference point with the following fields:
#   qout and qSup: Vectors holding the quantile function values and support points;
#   dout and dSup: Vectors holding the density function values and support points;
#   res_PGs: A list holding the results regarding principal geodesics.
# A: A matrix holding the amplitude factors with rows corresponding to subjects and columns to components.
# Lqf_normd: A list holding the quantile functions of normalized processes, where Lqf_normd[[j]][[i]] is a matrix with rows = time points in obsGrid and columns = quantile values on qSup.
# res_H: A list holding the subject-level deformation functions with the following fields:
#   H: A matrix with rows = subjects, columns = time points in obsGrid.
#   Hinv: Inverse of H, same structure.
#   obsGrid: Support grid of time points.
#   Loptns, LtimingWarp: Control options and timing for pairwise synchronization.
# res_eta_pre: A list holding preliminary estimates of component tempo trajectories (same structure as res_eta).
# res_G_pre: A list holding preliminary estimates of global deformation functions with fields:
#   H: List of matrices (one per component) with rows = subjects, cols = time points.
#   Hinv: Inverse of H, same structure.
#   Lqf_normd_aligned: Quantile functions of globally aligned normalized processes.
#   obsGrid, optns, timingWarp: Grid, options, and timing.
# res_tau: A list holding the latent template with the following fields:
#   qout: Matrix with rows = time points, columns = quantile values on qSup.
#   qSup: Quantile support vector.
#   dout: Matrix with rows = time points, columns = density values.
#   dSup: Matrix of density support points (one row per time point).
# res_Psi: A list holding the component-level deformation functions with the following fields:
#   H: List of vectors (one per component) giving Psi_j(t) on obsGrid.
#   Hinv: Inverse of Psi, same structure.
#   obsGrid: Support grid of time points.
# res_eta: A list holding the resolved component tempo trajectories (same structure as res_eta_pre).
# res_G: A list holding the final global deformation functions (same structure as res_G_pre).
# obsGrid: A vector holding the support grid of time points.
# qSup: A vector holding the support points of quantile functions.
# timing: Scalar holding the total computation time.
# pxcalign_f2m: A vector holding the values of the population-level cross-component alignment map from females' to males' on obsGrid.
# sxcalign_f2m: A matrix holding the subject-level cross-component alignment maps from females' to males' with rows = subjects, columns = time points in obsGrid.
