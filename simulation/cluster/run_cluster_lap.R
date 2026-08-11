# run_cluster_lap.R
# 4 methods: GLDM, CM, SRVF_KM, SRVF_pair (Laplacian/network case)
# Usage: Rscript run_cluster_lap.R --seed=1 --sigma_warp=0.5 --sigma_dist=0.5 --sigma_pert=0.5

suppressPackageStartupMessages({
  library(frechet); library(minqa); library(Matrix)
  library(expm); library(SMFilter); library(pracma)
  library(fdasrvf)
})
cat("[R] packages loaded.\n"); flush.console()

source("simulate_lap.R")
source("lnr.R")
source("fctns_estimation_lap.R")
source("fctns_evaluation.R")
source("fctns_srvf_lap.R")
cat("[R] Functions sourced.\n"); flush.console()

# ==================================================================
# Main
# ==================================================================
get_arg <- function(key, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", key, "="), "", hit[1])
}

seed <- as.integer(get_arg("seed", "1"))
w <- as.numeric(get_arg("sigma_warp", "0.5"))
d <- as.numeric(get_arg("sigma_dist", "0.5"))
p_noise <- as.numeric(get_arg("sigma_pert", "0.5"))
result_root <- get_arg("result_root", "results_nw_bothvar")

obsGrid <- seq(0, 1, length.out = 51)

cat(sprintf("[R] seed=%d, W=%.1f, D=%.1f, P=%.1f\n", seed, w, d, p_noise))

# ---- Create mask (fixed seed=123, same for all runs) ----
tau_obj <- make_tau(obsGrid = obsGrid, nnode = 10, p_edge = 0.5, seed = 123)
tau <- tau_obj$tau
mask <- tau_obj$mask

# ---- Timing: START ----
t_start <- proc.time()

# ---- Simulate data ----
simdata <- simulate_from_tau(
  tau = tau, obsGrid = obsGrid, mask = mask,
  n = 30, sigma_warp = w, sigma_dist = d,
  sigma_pert = p_noise, seed = seed
)

# ---- GLDM + CM ----
cat("[GLDM/CM] Running...\n"); flush.console()
gldmfit <- gldm_lap_once(obsGrid = obsGrid, Y = simdata$Y,
                         nknots = 4, method = "nearest", doLFR = TRUE, seed = seed)
gldm_metrics <- evaluate_gldm_lap_once(gldmfit, simdata, obsGrid)

# ---- Get Xstar and Xhat from gldmfit ----
Xstar <- gldmfit$norm_data
Xhat  <- gldmfit$Xhat

# ---- SRVF_KM ----
cat("[SRVF_KM] Running...\n"); flush.console()
step1_km <- tryCatch(srvf_step1_km_lap(Xstar, obsGrid, mask), error = function(e) {
  cat("[SRVF_KM] Step 1 error:", e$message, "\n"); NULL
})
srvf_km_metrics <- NULL
if (!is.null(step1_km)) {
  srvf_km <- tryCatch(
    srvf_pipeline_lap(Xstar, Xhat, obsGrid, step1_km$res_H, step1_km$res_G_pre, mask),
    error = function(e) { cat("[SRVF_KM] Pipeline error:", e$message, "\n"); NULL }
  )
  if (!is.null(srvf_km)) {
    srvf_km_metrics <- tryCatch(
      evaluate_srvf_lap(srvf_km, simdata, obsGrid),
      error = function(e) { cat("[SRVF_KM] Eval error:", e$message, "\n"); NULL }
    )
  }
}

# ---- SRVF_pair ----
cat("[SRVF_pair] Running...\n"); flush.console()
step1_pair <- tryCatch(srvf_step1_pair_lap(Xstar, obsGrid, mask), error = function(e) {
  cat("[SRVF_pair] Step 1 error:", e$message, "\n"); NULL
})
srvf_pair_metrics <- NULL
if (!is.null(step1_pair)) {
  srvf_pair <- tryCatch(
    srvf_pipeline_lap(Xstar, Xhat, obsGrid, step1_pair$res_H, step1_pair$res_G_pre, mask),
    error = function(e) { cat("[SRVF_pair] Pipeline error:", e$message, "\n"); NULL }
  )
  if (!is.null(srvf_pair)) {
    srvf_pair_metrics <- tryCatch(
      evaluate_srvf_lap(srvf_pair, simdata, obsGrid),
      error = function(e) { cat("[SRVF_pair] Eval error:", e$message, "\n"); NULL }
    )
  }
}

time_total <- (proc.time() - t_start)[3]

# ---- Collect all metrics ----
safe_val <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  as.numeric(x[1])
}
sv <- function(lst, nm) if (!is.null(lst)) lst[[nm]] else NA_real_

all_metrics <- list(
  # sWMISE
  sWMISE             = gldm_metrics$Hmise,
  sWMISE_srvf        = sv(srvf_pair_metrics, "sWMISE"),
  sWMISE_srvf_km     = sv(srvf_km_metrics, "sWMISE"),
  # gWMISEp (primitive G)
  gWMISEp            = gldm_metrics$gWMISE_gldm_p,
  gWMISEp_srvf       = sv(srvf_pair_metrics, "gWMISE_p"),
  gWMISEp_srvf_km    = sv(srvf_km_metrics, "gWMISE_p"),
  # gWMISEr (refined G)
  gWMISEr            = gldm_metrics$gWMISE_gldm_r,
  gWMISEr_cm         = gldm_metrics$gWMISE_cm,
  gWMISEr_srvf       = sv(srvf_pair_metrics, "gWMISE_r"),
  gWMISEr_srvf_km    = sv(srvf_km_metrics, "gWMISE_r"),
  # cWMISE
  cWMISE             = gldm_metrics$Pmise_new,
  cWMISE_cm          = gldm_metrics$Pmise_old,
  cWMISE_srvf        = sv(srvf_pair_metrics, "cWMISE"),
  cWMISE_srvf_km     = sv(srvf_km_metrics, "cWMISE"),
  # TISE
  TISE               = gldm_metrics$Tise_new,
  TISE_cm            = gldm_metrics$Tise_old,
  TISE_srvf          = sv(srvf_pair_metrics, "TISE"),
  TISE_srvf_km       = sv(srvf_km_metrics, "TISE"),
  # PMISE (primitive G)
  PMISE               = gldm_metrics$Xmise_aligned_p,
  PMISE_cm            = gldm_metrics$Xmise_aligned_old,
  PMISE_srvf          = sv(srvf_pair_metrics, "PMISE_pre"),
  PMISE_srvf_km       = sv(srvf_km_metrics, "PMISE_pre")
)

# ---- Save CSV ----
out_dir <- file.path(result_root, sprintf("W%.2f_D%.2f_P%.2f", w, d, p_noise))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_file <- file.path(out_dir, sprintf("seed_%03d.csv", seed))
df <- data.frame(seed = seed, sigma_warp = w, sigma_dist = d, sigma_pert = p_noise,
                 as.data.frame(lapply(all_metrics, safe_val)))
write.csv(df, out_file, row.names = FALSE)

cat(sprintf("\n[DONE] seed=%d | saved: %s (%.0fs)\n", seed, out_file, time_total))
