# =============================================================================
# File: analyze_lap.R
# Purpose: Run Monte Carlo comparison of GLDM, CM, SRVF_KM, and SRVF_pair
#          methods on Laplacian-valued functional data simulations.
# Dependencies: simulate_lap.R (or simulate_lap_phaseonly.R), lnr.R,
#               fctns_estimation_lap.R, fctns_evaluation.R, fctns_srvf_lap.R
# =============================================================================

# Suppress startup messages from attached packages for cleaner output
suppressPackageStartupMessages({
  library(frechet)    # Frechet mean and local covariance regression
  library(minqa)      # BOBYQA optimizer for warping estimation
  library(Matrix)     # Sparse matrix utilities
  library(expm)       # Matrix exponential (used in SPD geometry helpers)
  library(SMFilter)   # Provides FDist2() for squared Frobenius distance
  library(pracma)     # trapz() for numerical integration
  library(fdasrvf)    # SRVF-based multivariate curve alignment
  library(dplyr)      # Pipe operator %>% used for setwd
})

# Clear workspace and set working directory to the script's location
rm(list = ls())
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  rstudioapi::getSourceEditorContext()$path %>% dirname %>% setwd
}

# Load simulation, estimation, and evaluation helpers
source("simulate_lap.R")          # make_tau(), simulate_from_tau()
source("lnr.R")                   # lnr() local network regression (Zhou & Müller 2022)
source("fctns_estimation_lap.R")  # gldm_lap_once() and all estimation helpers
source("fctns_evaluation.R")      # compute_H_mise(), evaluate_gldm_lap_once(), etc.
source("fctns_srvf_lap.R")        # srvf_step1_km_lap(), srvf_pipeline_lap(), etc.

# =============================================================================
# Main simulation loop
# =============================================================================

#' Run Monte Carlo comparison of GLDM, CM, SRVF_KM, and SRVF_pair
#'
#' @description
#' For each Monte Carlo seed b = 1..B, generates Laplacian-valued functional
#' data using simulate_from_tau(), fits GLDM/CM and both SRVF baselines, and
#' records all ISE/MISE metrics per parameter combination.
#'
#' @param B        Integer; number of Monte Carlo replications.
#' @param n        Integer; number of subjects per replication.
#' @param sigma_warp Numeric vector; levels of subject-level warp noise.
#' @param sigma_dist Numeric vector; levels of nuisance time-distortion noise.
#' @param sigma_pert Numeric vector; levels of edge-weight perturbation noise.
#' @param obsGrid  Numeric vector; common observation time grid on [0, 1].
#' @param nnode    Integer; number of nodes in the graph.
#' @param p_edge   Numeric in (0,1); Bernoulli edge probability for the mask.
#'
#' @return A named list of length B; each element is a list with:
#'   - seed: the MC seed used
#'   - table: data.frame with one row per parameter combination and metric columns
run_srvf_comparison <- function(
    B = 250,
    n = 30,
    sigma_warp = c(0, 0.5, 1),
    sigma_dist = c(0, 0.5, 1),
    sigma_pert = c(0, 0.5, 1),
    obsGrid = seq(0, 1, length.out = 51),
    nnode = 10,
    p_edge = 0.5
) {
  total_start <- Sys.time()

  # Create the Laplacian template tau and edge mask ONCE with a fixed seed so
  # the graph structure is identical across all MC replications and noise levels.
  tau_obj <- make_tau(obsGrid = obsGrid, nnode = nnode, p_edge = p_edge, seed = 123)
  tau <- tau_obj$tau   # d x d x T template trajectory
  mask <- tau_obj$mask # d x d binary edge-presence mask

  # Build the full factorial grid over noise parameters
  param_grid <- expand.grid(
    sigma_warp = sigma_warp,
    sigma_dist = sigma_dist,
    sigma_pert = sigma_pert,
    stringsAsFactors = FALSE
  )
  n_combos <- nrow(param_grid)

  # Column names for all metric columns across methods
  # GLDM/CM columns: subject warp, global warp (primitive/refined), component
  #   warp, template, and pre-alignment PMISE
  col_gldm <- c("sWMISE",
                 "gWMISEp", "gWMISEr", "gWMISEr_cm",
                 "cWMISE", "cWMISE_cm",
                 "TISE", "TISE_cm",
                 "PMISE", "PMISE_cm")
  col_km   <- c("sWMISE_srvf_km",
                 "gWMISEp_srvf_km", "gWMISEr_srvf_km",
                 "cWMISE_srvf_km", "TISE_srvf_km",
                 "PMISE_srvf_km")
  col_pair <- c("sWMISE_srvf",
                 "gWMISEp_srvf", "gWMISEr_srvf",
                 "cWMISE_srvf", "TISE_srvf",
                 "PMISE_srvf")
  all_cols <- c(col_gldm, col_km, col_pair)

  results <- vector("list", B)

  for (b in seq_len(B)) {
    cat(sprintf("\n========== Seed %d/%d ==========\n", b, B))
    seed_start <- Sys.time()

    # Pre-allocate error matrix for all combos and all metric columns
    errs <- matrix(NA_real_, nrow = n_combos, ncol = length(all_cols),
                   dimnames = list(NULL, all_cols))

    for (g in seq_len(n_combos)) {
      # Extract noise levels for this combination
      sw <- param_grid$sigma_warp[g]
      sd <- param_grid$sigma_dist[g]
      sp <- param_grid$sigma_pert[g]

      cat(sprintf("  [%d/%d] W=%.1f D=%.1f P=%.1f ... ", g, n_combos, sw, sd, sp))

      # Simulate Laplacian-valued functional data for this noise combination
      simdata <- simulate_from_tau(
        tau = tau, obsGrid = obsGrid, mask = mask,
        n = n, sigma_warp = sw, sigma_dist = sd,
        sigma_pert = sp, seed = b
      )

      # ---- GLDM + CM ----
      # Fit the GLDM (and simultaneous CM baseline) on observed trajectories Y
      gldmfit <- gldm_lap_once(obsGrid = obsGrid, Y = simdata$Y,
                               method = "nearest",
                               doLFR = TRUE, seed = b)
      # Evaluate all GLDM/CM metrics against the known ground truth
      gldm_eval <- evaluate_gldm_lap_once(gldmfit, simdata, obsGrid)

      errs[g, "sWMISE"]      <- gldm_eval$Hmise
      errs[g, "gWMISEp"]     <- gldm_eval$gWMISE_gldm_p
      errs[g, "gWMISEr"]     <- gldm_eval$gWMISE_gldm_r
      errs[g, "gWMISEr_cm"]  <- gldm_eval$gWMISE_cm
      errs[g, "cWMISE"]      <- gldm_eval$Pmise_new
      errs[g, "cWMISE_cm"]   <- gldm_eval$Pmise_old
      errs[g, "TISE"]        <- gldm_eval$Tise_new
      errs[g, "TISE_cm"]     <- gldm_eval$Tise_old
      errs[g, "PMISE"]        <- gldm_eval$Xmise_aligned_p
      errs[g, "PMISE_cm"]     <- gldm_eval$Xmise_aligned_old

      # Extract amplitude-normalized and smoothed trajectories for SRVF baselines
      Xstar <- gldmfit$norm_data  # amplitude-normalized X*_{ij}(t)
      Xhat  <- gldmfit$Xhat       # lnr-smoothed X_{ij}(t)

      # ---- SRVF_KM ----
      # Step 1 uses multivariate Karcher mean from fdasrvf for H_i and G_ij
      tryCatch({
        step1_km <- srvf_step1_km_lap(Xstar, obsGrid, mask)
        srvf_km <- srvf_pipeline_lap(Xstar, Xhat, obsGrid,
                                     step1_km$res_H, step1_km$res_G_pre, mask)
        km_eval <- evaluate_srvf_lap(srvf_km, simdata, obsGrid)
        errs[g, "sWMISE_srvf_km"]  <- km_eval$sWMISE
        errs[g, "gWMISEp_srvf_km"] <- km_eval$gWMISE_p
        errs[g, "gWMISEr_srvf_km"] <- km_eval$gWMISE_r
        errs[g, "cWMISE_srvf_km"]  <- km_eval$cWMISE
        errs[g, "TISE_srvf_km"]    <- km_eval$TISE
        errs[g, "PMISE_srvf_km"]    <- km_eval$PMISE_pre
      }, error = function(e) {
        cat(sprintf("SRVF_KM error: %s ", e$message))
      })

      # ---- SRVF_pair ----
      # Step 1 uses pairwise curve_pair_align averaged over all reference pairs
      tryCatch({
        step1_pair <- srvf_step1_pair_lap(Xstar, obsGrid, mask)
        srvf_pair <- srvf_pipeline_lap(Xstar, Xhat, obsGrid,
                                       step1_pair$res_H, step1_pair$res_G_pre, mask)
        pair_eval <- evaluate_srvf_lap(srvf_pair, simdata, obsGrid)
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

    # Combine the parameter columns and error columns into a single data.frame
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

# =============================================================================
# Execute: run the Monte Carlo simulation
# =============================================================================
# This script runs the "bothvar" case (amplitude + phase variation in tau).
# For the phase-only variant (A_{ij} = 1 for all i,j, no amplitude noise),
# replace source("simulate_lap.R") above with source("simulate_lap_phaseonly.R")
# and change the output filename to "mc_tables_nw_phaseonly.rds".
res <- run_srvf_comparison(
  B = 2,                           # set to 250 for full run
  n = 30,
  sigma_warp = c(0, 0.5),         # subset for testing; use c(0, 0.5, 1) for full
  sigma_dist = c(0.5),
  sigma_pert = c(0.5),
  obsGrid = seq(0, 1, length.out = 51),
  nnode = 10,
  p_edge = 0.5
)

# Save results to a versioned RDS file in the shared results directory
dir.create("../results", showWarnings = FALSE)
saveRDS(res, "../results/mc_tables_nw_bothvar.rds")
cat("Saved: ../results/mc_tables_nw_bothvar.rds\n")

# To generate plots, run simulation/make_plots.R with data_type = "nw" and scenario = "bothvar".
