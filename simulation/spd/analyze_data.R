# =============================================================================
# File: analyze_data.R
# Purpose: Monte Carlo simulation comparing GLDM, CM, SRVF_KM, and SRVF_pair
#          alignment methods for SPD matrix-valued functional data.
# Dependencies: simulate_data.R, fctns_estimation.R, fctns_evaluation.R,
#               fctns_srvf.R, ../make_plots.R
#              Libraries: frechet, minqa, Matrix, expm, SMFilter, pracma,
#                         fdasrvf, dplyr
# =============================================================================
#
# Usage (local): source this file in R or run from RStudio.
# Adjust B, sigma grids, and obsGrid as needed.

# -----------------------------------------------------------------------------
# Package loading
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(frechet)    # Frechet regression and CovFMean for SPD data
  library(minqa)      # bobyqa optimizer for pairwise warping
  library(Matrix)     # Sparse matrix utilities
  library(expm)       # Matrix exponential/logarithm
  library(SMFilter)   # FDist2 (squared Frobenius distance)
  library(pracma)     # trapz (trapezoidal integration)
  library(fdasrvf)    # SRVF-based curve alignment baselines
  library(dplyr)      # pipe operator (%>%) for path manipulation
})

rm(list = ls())  # clear workspace before running simulation
# Set working directory to the location of this script (RStudio only)
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  rstudioapi::getSourceEditorContext()$path %>% dirname %>% setwd
}

# -----------------------------------------------------------------------------
# Source helper files
# -----------------------------------------------------------------------------
source("simulate_data.R")    # simulate_data(): generates synthetic SPD trajectories
source("fctns_estimation.R") # gldm_spd_once(): fits GLDM and CM models
source("fctns_evaluation.R") # evaluate_gldm_spd_once(), evaluate_srvf(): computes metrics
source("fctns_srvf.R")       # srvf_step1_km(), srvf_step1_pair(), srvf_pipeline()

# ==================================================================
# Main simulation function
# ==================================================================

#' Run Monte Carlo comparison of GLDM, CM, SRVF_KM, and SRVF_pair methods
#'
#' @description
#' For each simulation seed b = 1..B and each combination of noise parameters
#' (sigma_warp, sigma_dist, sigma_pert), generates synthetic SPD trajectories,
#' fits all four alignment methods, and records ISE/MISE metrics. Returns a
#' list of per-seed result tables for downstream plotting and summarization.
#'
#' @param B          Integer; number of Monte Carlo repetitions (seeds).
#' @param n          Integer; number of subjects per simulation.
#' @param sigma_warp Numeric vector; values of subject-level warp noise to evaluate.
#' @param sigma_dist Numeric vector; values of nuisance time-distortion noise.
#' @param sigma_pert Numeric vector; values of multiplicative SPD perturbation noise.
#' @param axis_variation Numeric; controls eigenvalue spread in the SPD template.
#' @param obsGrid    Numeric vector; common observation time grid on [0,1].
#' @param nknots     Integer; number of interior knots for piecewise-linear warping.
#'
#' @return A named list of length B. Each element is a list with:
#'   - seed: the seed index b
#'   - table: a data.frame (n_combos rows x all metric columns) for that seed
run_srvf_comparison <- function(
    B = 250,
    n = 30,
    sigma_warp = c(0, 0.5, 1),
    sigma_dist = c(0, 0.5, 1),
    sigma_pert = c(0, 0.5, 1),
    axis_variation = 10,
    obsGrid = seq(0, 1, length.out = 51),
    nknots = 4
) {
  total_start <- Sys.time()

  # Build the full factorial grid of noise parameter combinations
  param_grid <- expand.grid(
    sigma_warp = sigma_warp,
    sigma_dist = sigma_dist,
    sigma_pert = sigma_pert,
    stringsAsFactors = FALSE
  )
  n_combos <- nrow(param_grid)

  # -----------------------------------------------------------------------------
  # Define output metric column names for each method
  # -----------------------------------------------------------------------------

  # GLDM and CM error metrics (computed jointly from gldm_spd_once output)
  col_gldm <- c("sWMISE",           # subject-level warp MISE (GLDM)
                 "gWMISEp", "gWMISEr", "gWMISEr_cm",  # global warp MISE: primitive / refined / CM
                 "cWMISE", "cWMISE_cm",                # component-level warp MISE: GLDM / CM
                 "TISE", "TISE_cm",                    # template ISE: GLDM / CM
                 "PMISE", "PMISE_cm")                  # aligned trajectory MISE: GLDM / CM
  # SRVF Karcher-mean baseline metrics
  col_km   <- c("sWMISE_srvf_km",
                 "gWMISEp_srvf_km", "gWMISEr_srvf_km",
                 "cWMISE_srvf_km", "TISE_srvf_km",
                 "PMISE_srvf_km")
  # SRVF pairwise alignment baseline metrics
  col_pair <- c("sWMISE_srvf",
                 "gWMISEp_srvf", "gWMISEr_srvf",
                 "cWMISE_srvf", "TISE_srvf",
                 "PMISE_srvf")
  all_cols <- c(col_gldm, col_km, col_pair)

  results <- vector("list", B)

  # -----------------------------------------------------------------------------
  # Outer loop: Monte Carlo seeds
  # -----------------------------------------------------------------------------
  for (b in seq_len(B)) {
    cat(sprintf("\n========== Seed %d/%d ==========\n", b, B))
    seed_start <- Sys.time()

    # Pre-allocate error matrix for all parameter combinations in this seed
    errs <- matrix(NA_real_, nrow = n_combos, ncol = length(all_cols),
                   dimnames = list(NULL, all_cols))

    # ---- Inner loop: parameter combinations ----
    for (g in seq_len(n_combos)) {
      sw <- param_grid$sigma_warp[g]  # subject-level warp noise level
      sd <- param_grid$sigma_dist[g]  # nuisance distortion noise level
      sp <- param_grid$sigma_pert[g]  # SPD perturbation noise level

      cat(sprintf("  [%d/%d] W=%.1f D=%.1f P=%.1f ... ", g, n_combos, sw, sd, sp))

      # Fix RNG for reproducible data generation across methods
      RNGversion("3.5.0"); set.seed(2025L)
      simdata <- simulate_data(n = n, sigma_warp = sw, sigma_dist = sd,
                               sigma_pert = sp, axis_variation = axis_variation,
                               obsGrid = obsGrid, seed = b)

      # ---- GLDM + CM: fit both methods jointly via gldm_spd_once ----
      gldmfit <- gldm_spd_once(obsGrid = obsGrid, Y = simdata$Y,
                               nknots = nknots, method = "nearest",
                               doLFR = TRUE, seed = b)
      # Compute all ISE/MISE metrics for GLDM and CM
      gldm_eval <- evaluate_gldm_spd_once(gldmfit, simdata, obsGrid)

      # Store GLDM and CM metrics into the error matrix
      errs[g, "sWMISE"]      <- gldm_eval$Hmise           # subject warp MISE (both methods share H)
      errs[g, "gWMISEp"]     <- gldm_eval$gWMISE_gldm_p   # global warp MISE (primitive G)
      errs[g, "gWMISEr"]     <- gldm_eval$gWMISE_gldm_r   # global warp MISE (refined G, GLDM)
      errs[g, "gWMISEr_cm"]  <- gldm_eval$gWMISE_cm       # global warp MISE (CM)
      errs[g, "cWMISE"]      <- gldm_eval$Pmise_new        # component warp MISE (GLDM)
      errs[g, "cWMISE_cm"]   <- gldm_eval$Pmise_old        # component warp MISE (CM)
      errs[g, "TISE"]        <- gldm_eval$Tise_new         # template ISE (GLDM)
      errs[g, "TISE_cm"]     <- gldm_eval$Tise_old         # template ISE (CM)
      errs[g, "PMISE"]        <- gldm_eval$Xmise_aligned_p  # aligned trajectory MISE (GLDM, pre-refinement G)
      errs[g, "PMISE_cm"]     <- gldm_eval$Xmise_aligned_old # aligned trajectory MISE (CM)

      # Extract amplitude-normalized and smoothed trajectories for SRVF baselines
      Xstar <- gldmfit$norm_data  # amplitude-normalized trajectories X*
      Xhat  <- gldmfit$Xhat       # smoothed (LocCovReg) trajectories

      # ---- SRVF_KM: Karcher-mean-based SRVF alignment ----
      tryCatch({
        # Step 1: estimate H_i and G_pre via multivariate Karcher mean
        step1_km <- srvf_step1_km(Xstar, obsGrid)
        # Steps 2-4: estimate tau, Psi, G_ref using the shared SRVF pipeline
        srvf_km <- srvf_pipeline(Xstar, Xhat, obsGrid,
                                 step1_km$res_H, step1_km$res_G_pre)
        km_eval <- evaluate_srvf(srvf_km, simdata, obsGrid)
        errs[g, "sWMISE_srvf_km"]  <- km_eval$sWMISE
        errs[g, "gWMISEp_srvf_km"] <- km_eval$gWMISE_p
        errs[g, "gWMISEr_srvf_km"] <- km_eval$gWMISE_r
        errs[g, "cWMISE_srvf_km"]  <- km_eval$cWMISE
        errs[g, "TISE_srvf_km"]    <- km_eval$TISE
        errs[g, "PMISE_srvf_km"]    <- km_eval$PMISE_pre
      }, error = function(e) {
        cat(sprintf("SRVF_KM error: %s ", e$message))
      })

      # ---- SRVF_pair: pairwise-alignment-based SRVF ----
      tryCatch({
        # Step 1: estimate H_i and G_pre via pairwise curve alignment
        step1_pair <- srvf_step1_pair(Xstar, obsGrid)
        # Steps 2-4: shared SRVF pipeline (same as KM, different Step 1 outputs)
        srvf_pair <- srvf_pipeline(Xstar, Xhat, obsGrid,
                                   step1_pair$res_H, step1_pair$res_G_pre)
        pair_eval <- evaluate_srvf(srvf_pair, simdata, obsGrid)
        errs[g, "sWMISE_srvf"]  <- pair_eval$sWMISE
        errs[g, "gWMISEp_srvf"] <- pair_eval$gWMISE_p
        errs[g, "gWMISEr_srvf"] <- pair_eval$gWMISE_r
        errs[g, "cWMISE_srvf"]  <- pair_eval$cWMISE
        errs[g, "TISE_srvf"]    <- pair_eval$TISE
        errs[g, "PMISE_srvf"]    <- pair_eval$PMISE_pre
      }, error = function(e) {
        cat(sprintf("SRVF_pair error: %s ", e$message))
      })

      cat("done.\n")
    }

    # Build table for this seed: bind parameter columns with error matrix
    seed_table <- cbind(
      data.frame(seed = b,
                 sigma_warp = param_grid$sigma_warp,
                 sigma_dist = param_grid$sigma_dist,
                 sigma_pert = param_grid$sigma_pert),
      as.data.frame(errs)
    )

    results[[b]] <- list(seed = b, table = seed_table)
    names(results)[b] <- sprintf("seed_%03d", b)

    seed_dur <- difftime(Sys.time(), seed_start, units = "mins")
    cat(sprintf("  Seed %d done in %.1f min.\n", b, as.numeric(seed_dur)))
  }

  total_dur <- difftime(Sys.time(), total_start, units = "hours")
  cat(sprintf("\nAll done in %.1f hours.\n", as.numeric(total_dur)))

  return(results)
}

# ==================================================================
# Run (adjust parameters as needed)
# ==================================================================
# This example runs the "bothvar" case (amplitude + phase variation).
# For phase-only (no amplitude variation), replace gldm_spd_once with
# gldm_spd_phaseonly in the main loop above, and rename the output
# file to "mc_tables_spd_phaseonly.rds".
res <- run_srvf_comparison(
  B = 2,                           # set to 250 for full run
  n = 30,
  sigma_warp = c(0, 0.5),         # subset for testing; use c(0, 0.5, 1) for full
  sigma_dist = c(0.5),
  sigma_pert = c(0.5),
  axis_variation = 10,
  obsGrid = seq(0, 1, length.out = 51),
  nknots = 4
)

# -----------------------------------------------------------------------------
# Save results
# -----------------------------------------------------------------------------
dir.create("../results", showWarnings = FALSE)  # create results directory if needed
saveRDS(res, "../results/mc_tables_spd_bothvar.rds")
cat("Saved: ../results/mc_tables_spd_bothvar.rds\n")

# To generate plots, run simulation/make_plots.R with data_type = "spd" and scenario = "bothvar".
