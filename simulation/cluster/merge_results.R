#' Merge CSV results from cluster runs into mc_tables rds format
#'
#' @description
#' Reads per-seed CSV files from cluster result directories and combines them
#' into the mc_tables rds format used by make_plots.R for generating boxplots.
#'
#' @param result_dir Directory containing W*_D*_P*/ subdirectories with seed CSVs.
#' @param out_file   Path to save the output rds file.
#' @param n_seeds    Number of seeds (default 250).
#'
#' @return Invisibly returns the mc_tables list.

merge_cluster_results <- function(result_dir, out_file, n_seeds = 250) {
  combo_dirs <- list.dirs(result_dir, recursive = FALSE, full.names = TRUE)
  cat(sprintf("Found %d combo directories in %s\n", length(combo_dirs), result_dir))

  # Columns to drop (timing/metadata from cluster runs)
  drop_cols <- c("time_gldm_cm_sec", "time_srvf_km_sec", "time_srvf_pair_sec",
                 "time_total_sec", "timestamp")

  results <- vector("list", n_seeds)
  for (s in seq_len(n_seeds)) {
    sf <- sprintf("seed_%03d.csv", s)
    rows <- list()
    for (combo_dir in combo_dirs) {
      f <- file.path(combo_dir, sf)
      if (file.exists(f)) {
        rows[[length(rows) + 1]] <- read.csv(f, stringsAsFactors = FALSE)
      }
    }
    if (length(rows) > 0) {
      tbl <- do.call(rbind, rows)
      tbl <- tbl[, !names(tbl) %in% drop_cols]
      results[[s]] <- list(seed = s, table = tbl)
    }
    names(results)[s] <- sprintf("seed_%03d", s)
  }

  n_complete <- sum(sapply(results, function(x) {
    if (!is.null(x)) nrow(x$table) == 27 else FALSE
  }))
  cat(sprintf("Seeds with 27 combos: %d/%d\n", n_complete, n_seeds))

  saveRDS(results, out_file)
  cat(sprintf("Saved: %s\n", out_file))
  invisible(results)
}

# ====================================================================
# Usage examples:
# merge_cluster_results("results_spd_bothvar", "mc_tables_spd_bothvar.rds")
# merge_cluster_results("results_spd_phaseonly", "mc_tables_spd_phaseonly.rds")
# merge_cluster_results("results_nw_bothvar", "mc_tables_nw_bothvar.rds")
# merge_cluster_results("results_nw_phaseonly", "mc_tables_nw_phaseonly.rds")
# ====================================================================
