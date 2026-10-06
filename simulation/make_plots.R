# =============================================================================
# File: make_plots.R
# Purpose: Generate boxplot figures comparing GLDM, CM, and SRVF methods
#          across Monte Carlo replicates and noise-level combinations.
# Inputs:  results/mc_tables_{data_type}_{scenario}.rds  -- merged MC results
#          (produced by cluster/merge_results.R)
# Outputs: figures/{scenario}_{data_type}_boxplot_{metric}.pdf for each metric
#          in {sWMISE, gWMISEp, gWMISEr, cWMISE, TISE, PMISE}
# Dependencies: purrr, plyr, dplyr, tidyr, reshape2, ggplot2, patchwork,
#               cowplot, colorspace, ggrepel, readr, magrittr, scales,
#               rgl, corrplot
# =============================================================================

# -----------------------------------------------------------------------------
# Package loading
# -----------------------------------------------------------------------------

library(purrr)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(patchwork)
library(cowplot)
library(colorspace)
library(ggrepel)
library(readr)
library(magrittr)
library(scales)   # parse_format() for parsed axis labels
library(rgl)
library(corrplot)

# -----------------------------------------------------------------------------
# Path setup and scenario selection
# -----------------------------------------------------------------------------

# Set working directory to the script's location (RStudio only)
rstudioapi::getSourceEditorContext()$path %>% dirname %>% setwd

OUT_DIR <- "results"   # directory containing merged Monte Carlo rds files
FIG_DIR <- "figures"   # directory where PDF figures will be written
if ( !dir.exists(FIG_DIR) ) {
  dir.create( FIG_DIR )
}

# data_type selects the object space:  "spd" = SPD matrices (Case 1),
#                                       "nw"  = graph Laplacians (Case 2)
# data_type <- "nw"
data_type <- "spd"

# scenario selects the variation setting:  "bothvar"   = amplitude + phase,
#                                           "phaseonly" = phase variation only
# scenario <- "phaseonly"
scenario <- "bothvar"

# -----------------------------------------------------------------------------
# Load and reshape simulation results
# -----------------------------------------------------------------------------

# Each rds file is a list (one entry per noise combination) where each entry
# contains a "table" data frame with one row per MC replicate and one column
# per evaluation metric.
res <- readRDS( file.path( OUT_DIR, paste0("mc_tables_",data_type,"_",scenario,".rds") ) )

# Bind all noise-combination tables into a single data frame and add
# human-readable factor labels for the three noise parameters.
df <- lapply( res, `[[`, 'table' ) %>% bind_rows() %>%
  mutate(
    sigma_dist_label = paste("sigma[R] ==", as.character(sigma_dist)),
    sigma_pert_label = paste("sigma[P] ==", as.character(sigma_pert)),
    sigma_warp_label = paste("sigma[H] ==", as.character(sigma_warp))
  )

# -----------------------------------------------------------------------------
# Plotting helper
# -----------------------------------------------------------------------------

#' Generate and save a side-by-side boxplot for one or more metric columns
#'
#' @description
#' Pivots the simulation results to long format, filters extreme values (if
#' \code{xmax} is provided), and produces a boxplot faceted by
#' sigma_pert (columns) x sigma_dist (rows), with sigma_warp on the y-axis.
#' Omitted extreme values are written to a CSV in \code{OUT_DIR} for reference.
#'
#' @param variables  Character vector of column names in \code{df} to compare
#'                   (e.g., \code{c("gWMISEr", "gWMISEr_cm", "gWMISEr_srvf")}).
#' @param fill_cols  Character vector of fill colors, one per element of
#'                   \code{variables} (same order).
#' @param fill_labs  Character vector of legend labels for each method; defaults
#'                   to \code{variables} if \code{NULL}.
#' @param var_label  Character; metric name used for the x-axis label and the
#'                   output filename (e.g., \code{"gWMISEr"}).
#' @param xmax       Numeric or \code{NULL}; if provided, values above this
#'                   threshold are excluded from the plot and saved to CSV.
#' @param scales     Passed to \code{ggplot2::facet_grid()} scales argument;
#'                   use \code{"fixed"} (default) or \code{"free_x"}.
#' @param filename_fun Function mapping \code{var_label} to an output filename;
#'                   defaults to \code{"{scenario}_{data_type}_boxplot_{var_label}.pdf"}.
#'
#' @return Invisibly returns \code{NULL}; the figure is saved to \code{FIG_DIR}.
make_boxplot_for_multiple_variables <- function (
    variables = c("Tise_old", "Tise_new"),
    fill_cols = c("grey40","white"),
    fill_labs = NULL,
    var_label="TISE", xmax = NULL, scales = "fixed",
    filename_fun = function (var_label) { paste0(scenario,"_",data_type,"_boxplot_", var_label, ".pdf" ) }
) {
  # Pivot the selected metric columns to long format for grouped boxplots
  df_long <- df %>%
    pivot_longer(
      cols = variables,
      names_to = "type",
      values_to = "value"
    )
  if( !is.null(xmax) ) {
    # Exclude extreme values that would compress the main distribution visually;
    # the omitted rows are written to CSV so they are not silently dropped.
    dfpl <- df_long %>% filter(value <= xmax)
    df_omitted <- df_long %>% filter(value > xmax)
    write.csv(
      df_omitted,
      file.path( OUT_DIR, paste0(data_type,"_omitted_data_", var_label, ".csv" ) ),
      row.names = FALSE
    )
    cat("Omitted:")
    print( table(df_omitted$type) )
  } else {
    dfpl <- df_long
  }

  if ( is.null( fill_labs ) ) {
    fill_labs <- variables
  }

  pl <-
    dfpl %>%
    mutate( type = factor(type, levels = variables) ) %>%
    ggplot() +
    geom_boxplot( aes( y = sigma_warp_label, x= value, fill = type), size = 0.125, outlier.size = 0.125 ) +
    scale_y_discrete(labels = parse_format(), limits = rev) +  # parsed Greek labels, reversed order
    facet_grid(
      sigma_dist_label ~ sigma_pert_label,
      labeller = label_parsed, scales = scales
    ) +
    scale_fill_manual(
      values = fill_cols,
      labels = fill_labs
    ) +
    labs( x = var_label, y = NULL, fill = NULL ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "top",
      legend.key.spacing.x = unit(15, "pt"),
      legend.title = element_text(size=rel(1)),
      legend.text = element_text(size=rel(.9)),
      strip.text = element_text(size=rel(1.1)),
      strip.background = element_blank(),
      panel.grid = element_line(linewidth = 0.3, linetype = 2),
      plot.margin = margin(1,3,1,1, unit = "pt"),
      axis.title = element_text(size=rel(1.1)),
      axis.text.x = element_text(size=rel(1)),
      axis.text.y = element_text(size=rel(1.3))
    )
  ggsave(
    filename = filename_fun(var_label),
    plot = pl, path = FIG_DIR,
    width = 9, height = 2+length(variables)
  )
}

# -----------------------------------------------------------------------------
# Tau visualization helpers
# -----------------------------------------------------------------------------

#' Visualize the SPD template trajectory as a sequence of 3D ellipsoids
#'
#' @param tau      d x d x T array; SPD template trajectory on the time grid.
#' @param fig_dir  Character; output directory for the PNG file.
#' @param filename Character; output filename (default \code{"spd_tau.png"}).
#' @param spacing  Numeric; gap between consecutive ellipsoids on the x-axis (default 4).
#' @param alpha    Numeric in (0,1); ellipsoid transparency (default 0.6).
#'
#' @return Invisibly returns NULL; the PNG is saved to \code{file.path(fig_dir, filename)}.
plot_tau_spd <- function(tau, fig_dir = FIG_DIR, filename = "spd_tau.png",
                         spacing = 4, alpha = 0.6) {
  n_matrices <- dim(tau)[3]

  # Blue-to-red color gradient across time slices (t=0 → blue, t=1 → red)
  color_base <- rev(c("#0033FF", "#66DDFF", "#FFCC33", "#FF0000"))
  colors <- colorRampPalette(color_base)(n_matrices)

  # Open rgl window; windowRect controls the pixel dimensions of the PNG output
  open3d(windowRect = c(20, 20, 820, 110))
  bg3d("white")

  for (i in 1:n_matrices) {
    # Represent the i-th SPD matrix as the ellipsoid {x : x'M^{-1}x = 1}
    mesh <- ellipse3d(tau[,,i], t = 1)
    # Translate along x-axis so consecutive ellipsoids don't overlap
    mesh_moved <- translate3d(mesh, x = (i-1) * spacing, y = 0, z = 0)
    shade3d(mesh_moved, col = colors[i], alpha = alpha)
    wire3d(mesh_moved, col = "black", alpha = 0.1)
  }

  view3d(zoom = 0.1)
  rgl.snapshot(file.path(fig_dir, filename), fmt = "png")
  close3d()
}

#' Visualize the Laplacian template trajectory as a sequence of adjacency heatmaps
#'
#' @param tau      d x d x T array; Laplacian template trajectory on the time grid.
#' @param fig_dir  Character; output directory for the PDF file.
#' @param filename Character; output filename (default \code{"nw_tau.pdf"}).
#'
#' @return Invisibly returns NULL; the PDF is saved to \code{file.path(fig_dir, filename)}.
plot_tau_nw <- function(tau, fig_dir = FIG_DIR, filename = "nw_tau.pdf") {
  n_matrices <- dim(tau)[3]
  tau <- tau[,,seq(1, n_matrices, 2)]
  n_matrices <- dim(tau)[3]

  ncol <- ceiling(n_matrices / 1)

  # Greyscale palette: white (low weight) → dark grey (high weight)
  col_palette <- grey(seq(1, 0.5, length = 201))

  # Convert Laplacians to adjacency matrices (A = D - L) and find global
  # color scale limits so all panels share the same scale
  adj_list <- list()
  min_val <- Inf; max_val <- -Inf
  for (i in 1:n_matrices) {
    mat <- diag(diag(tau[,,i])) - tau[,,i]
    adj_list[[i]] <- mat
    r <- range(mat)
    if (r[1] < min_val) min_val <- r[1]
    if (r[2] > max_val) max_val <- r[2]
  }
  global_lim <- c(min_val, max_val)

  # Layout: row 1 = shared color bar legend, row 2 = heatmap panels
  layout_matrix <- matrix(c(rep(1, ncol), seq_len(ncol) + 1), nrow = 2, byrow = TRUE)

  pdf(file.path(fig_dir, filename), width = ncol * 2, height = 2.9)
  layout(layout_matrix, heights = c(0.3, 1))

  # Draw shared color bar
  par(mar = c(1.8, 40, .1, 40))
  plot(c(0,1), c(0,1), type = "n", axes = FALSE, xlab = "", ylab = "", main = "")
  rasterImage(as.raster(matrix(col_palette, nrow = 1)), 0.1, 0.3, 0.9, 0.8)
  axis(1, at = seq(0.1, 0.9, length.out = 5),
       labels = round(seq(min_val, max_val, length.out = 5), 2),
       pos = 0.3, cex.axis = 2.8)

  # Draw one heatmap per (subsampled) time slice
  par(mgp = c(1.8, .5, 0), mar = c(.5, .5, .5, .5))
  for (i in 1:n_matrices) {
    corrplot::corrplot(adj_list[[i]],
                       method = "color", is.corr = FALSE,
                       col = col_palette,
                       title = paste0("t = ", (i-1)/(n_matrices-1)),
                       mar = c(0, 0, 2.2, 0), cex.main = 3, font.main = 1,
                       cl.lim = global_lim,
                       cl.pos = "n",   # hide per-panel color legend (shared above)
                       tl.pos = "n")   # hide row/col labels
    rect(0.5, 0.5, ncol(adj_list[[i]]) + 0.5, ncol(adj_list[[i]]) + 0.5, lwd = 1)
  }
  dev.off()
}

# Fill colors for GLDM (white), CM (grey80), and SRVF (grey40)
colors <- c("white","grey80","grey40") %>%
  set_names( c("GLDM","CM", "SRVF") )

# =============================================================================
# Figures: one PDF per evaluation metric
# =============================================================================

# sWMISE: subject-level warping-function MISE (GLDM vs. SRVF only;
#          CM does not estimate subject-level warping functions)
make_boxplot_for_multiple_variables(
  variables = c("sWMISE", "sWMISE_srvf"),
  fill_cols = colors[c("GLDM","SRVF")] %>% set_names(NULL),
  fill_labs = c("GLDM","SRVF"),
  var_label="sWMISE"
)

# gWMISEp: global warping-function MISE (phase component; GLDM vs. SRVF)
make_boxplot_for_multiple_variables(
  variables = c("gWMISEp", "gWMISEp_srvf"),
  fill_cols = colors[c("GLDM","SRVF")] %>% set_names(NULL),
  fill_labs = c("GLDM","SRVF"),
  var_label="gWMISEp"
)

# gWMISEr: global warping-function MISE (amplitude component; all three methods)
make_boxplot_for_multiple_variables(
  variables = c("gWMISEr", "gWMISEr_cm", "gWMISEr_srvf"),
  fill_cols = colors %>% set_names(NULL),
  fill_labs = names(colors),
  var_label="gWMISEr"
)

# cWMISE: cross-sectional WMISE for estimated component tempos (all three methods)
make_boxplot_for_multiple_variables(
  variables = c("cWMISE", "cWMISE_cm", "cWMISE_srvf"),
  fill_cols = colors %>% set_names(NULL),
  fill_labs = names(colors),
  var_label="cWMISE"
)

# TISE: template ISE — integrated squared error of the estimated template (all three methods)
make_boxplot_for_multiple_variables(
  variables = c("TISE", "TISE_cm", "TISE_srvf"),
  fill_cols = colors %>% set_names(NULL),
  fill_labs = names(colors),
  var_label="TISE"
)

# PMISE: principal geodesic MISE — accuracy of the estimated principal geodesics (all three methods)
make_boxplot_for_multiple_variables(
  variables = c("PMISE", "PMISE_cm", "PMISE_srvf"),
  fill_cols = colors %>% set_names(NULL),
  fill_labs = names(colors),
  var_label="PMISE"
)
