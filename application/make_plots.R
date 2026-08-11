# =============================================================================
# File: make_plots.R
# Purpose: Generate all figures for the GLDM mortality application section,
#          reading fitted model results from results/ and writing PDFs to figures/.
# Inputs:  results/mort_preprocessed.rds  -- preprocessed mortality data
#          results/mort_ldmfit.rds         -- fitted GLDM object
#          GDP_MPD_1988.csv                -- 1988 GDP per capita (Maddison Project)
# Outputs: figures/mort_data.pdf, mort_amp.pdf, mort_template_ref_comp_tempo.pdf,
#          mort_comp_wf.pdf, mort_pxc_align.pdf, mort_subj_wf_per_country.pdf,
#          mort_sxc_align_per_country.pdf, mort_subj_wf_sxc_align_together.pdf,
#          mort_subj_nti.pdf, mort_subj_rti.pdf
# Dependencies: purrr, plyr, dplyr, tidyr, reshape2, ggplot2, patchwork,
#               cowplot, colorspace, ggrepel, readr
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

# -----------------------------------------------------------------------------
# Path setup
# -----------------------------------------------------------------------------

# Set working directory to the script's location (RStudio only)
rstudioapi::getSourceEditorContext()$path %>% dirname %>% setwd

OUT_DIR <- "results"   # directory containing fitted model rds files
FIG_DIR <- "figures"   # directory where PDF figures will be written
if ( !dir.exists(FIG_DIR) ) {
  dir.create( FIG_DIR )
}

# -----------------------------------------------------------------------------
# Load results
# -----------------------------------------------------------------------------

# Preprocessed mortality data: histograms, densities, and quantile functions
dat <- readRDS( file.path(OUT_DIR, "mort_preprocessed.rds") )
# Full GLDM fit: omega_0, A, H_i, eta_j, tau, Psi_j, G_i, and alignment maps
ldmfit <- readRDS( file.path(OUT_DIR, "mort_ldmfit.rds") )

# External covariate: 1988 GDP per capita (Maddison Project Database)
# Used to color-code country curves by economic development
gdp_pc <- read.csv('GDP_MPD_1988.csv')

# -----------------------------------------------------------------------------
# Country name helpers
# -----------------------------------------------------------------------------

# Shorten three verbose names to fit plot labels
shorter_names <- c(
  "Czechia", "UK", "USA"
) %>% set_names(
  c( "Czech Republic", "United Kingdom", "United States" )
)
dat$subj_names <-
  dat$subj_names %>% recode( !!!shorter_names ) %>%
  set_names( dat$subj_labs )
gdp_pc <- gdp_pc %>% mutate( country = country %>% recode(!!!shorter_names) )

# Recode three-letter HMD country codes to display names in a data frame
#
# @param df  Data frame with a column named "country" holding HMD abbreviations.
# @return    The same data frame with "country" replaced by display-ready names.
recode_country <- function (df) {
  df %>% mutate(
    country = country %>% recode( !!! dat$subj_names )
  )
}

# -----------------------------------------------------------------------------
# Shared color scales
# -----------------------------------------------------------------------------

# Blue-to-orange gradient used for year coloring (early = blue, recent = orange)
portland_colors <- c("#003366", "#3366CC", "#66CCFF", "#FFCC33", "#FF6600")
# Reversed palette used for GDP coloring (warm = high income, cool = low income)
portland_colors_rev <- rev( portland_colors )

# Standard (wide) color bars for per-country panel plots
col_scale <- scale_color_gradientn(
  colors = portland_colors,
  guide = guide_colorbar(
    title.position = "top",
    barwidth = 20,
    barheight = 0.6
  )
)
col_scale_rev <- scale_color_gradientn(
  colors = portland_colors_rev,
  guide = guide_colorbar(
    title.position = "top",
    barwidth = 20,
    barheight = 0.6
  )
)

# Narrower color bars for single-panel plots
col_scale_small <- scale_color_gradientn(
  colors = portland_colors,
  guide = guide_colorbar(
    title.position = "top",
    barwidth = 13,
    barheight = 0.4
  )
)
col_scale_rev_small <- scale_color_gradientn(
  colors = portland_colors_rev,
  guide = guide_colorbar(
    title.position = "top",
    barwidth = 13,
    barheight = 0.4
  )
)

# =============================================================================
# Figure: mort_data.pdf
# Purpose: Raw age-at-death density trajectories for each country, colored by
#          year.  Female and male panels are shown side by side with a shared
#          color bar and a shared y-axis range.
# =============================================================================

Ldens <- dat$Ldens
# Flatten the nested list (sex -> country -> year -> density) into a tidy frame
df <- map_dfr(names(Ldens), function(sex) {
  map_dfr(names(Ldens[[sex]]), function(country) {
    map_dfr(names(Ldens[[sex]][[country]]), function(year) {
      data.frame(
        sex = sex,
        country = country,
        year = as.numeric(year),
        x = Ldens[[sex]][[country]][[year]]$x,
        y = Ldens[[sex]][[country]][[year]]$y
      )
    })
  })
})

df <- df %>% recode_country()
df$sex <- factor(df$sex, levels = c("female", "male"))

# Shared y-axis range so female and male panels are directly comparable
y_range <- range( df$y, na.rm = TRUE )

# Shared theme for the two country-grid panels
pltheme <- theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size=rel(1.2)),
    panel.grid = element_line(linewidth = 0.3, linetype = 2),
    plot.margin = margin(1,6,1,1, unit = "pt"),
    axis.title = element_text(size=rel(1.2)),
    axis.text = element_text(size=rel(1))
  )

# Female plot
p_female <- ggplot(filter(df, sex == "female"),
                   aes(x = x, y = y, color = year, group = year)) +
  geom_line(linewidth = 0.3, alpha = 0.9) +
  facet_wrap(~ country, ncol = 4) +
  scale_y_continuous(limits = y_range) +
  col_scale +
  labs(x = "Age-at-death (year)", y = "Density", title = "Females") +
  pltheme

# Male plot
p_male <- ggplot(filter(df, sex == "male"),
                 aes(x = x, y = y, color = year, group = year)) +
  geom_line(linewidth = 0.3, alpha = 0.9) +
  facet_wrap(~ country, ncol = 4) +
  scale_y_continuous(limits = y_range) +
  col_scale +
  labs(x = "Age-at-death (year)", y = "Density", title = "Males") +
  pltheme

# Extract the shared legend from a dummy plot, then assemble the final layout
legend_plot <- ggplot(df, aes(x = x, y = y, color = year)) +
  geom_line() +
  col_scale +
  labs(color = "Year") +
  theme_void() +
  theme(
    legend.position = "top",
    legend.title = element_text(size=rel(1.2)),
    legend.text = element_text(size=rel(1.1)),
    legend.margin = margin(t = -15, b = -15),
    legend.box.margin = margin(-10, 0, -10, 0)
  )

legend <- cowplot::get_legend(legend_plot)

pl <- plot_grid(
  legend,
  plot_grid(p_female, p_male, ncol = 2, align = "v"),
  ncol = 1,
  rel_heights = c(0.08, 1)
)

ggsave( plot = pl, filename = "mort_data.pdf", path = FIG_DIR,
        width = 10.5, height = 12 )

# =============================================================================
# Figure: mort_amp.pdf
# Purpose: Scatter plot of subject-level amplitude factors A_i (females vs.
#          males) with a fitted linear trend.  Each point is one country.
#          Points above the identity line indicate larger female amplitude.
# =============================================================================

# Fit a linear model to quantify the female-male amplitude relationship
lmod <- lm(female~male, data=as.data.frame(ldmfit$A))
lmod %>% summary()

# Convert the amplitude matrix to a tidy data frame for ggplot
df <- data.frame( country = rownames(ldmfit$A), ldmfit$A ) %>%
  recode_country()

pl <- ggplot(df) +
  geom_smooth(aes(x = male, y = female), method = lm, formula = y~x,
              linewidth = .9, color = "grey50", fill = "grey80") +
  geom_point(aes(x = male, y = female), size = 1.5) +
  geom_text_repel(
    aes(x = male, y = female, label = country),
    size = 4.25,
    box.padding = 0.3,
    point.padding = 0.1,
    segment.color = "gray65",
    force = 60,
    force_pull = 0.05,
    min.segment.length = unit(0.2, "lines"),
    max.overlaps = Inf
  ) +
  geom_abline(intercept = 0, slope = 1, linewidth = .9, linetype = "dashed", color = "gray75") +
  labs(
    x = "Males' amplitude factor",
    y = "Females' amplitude factor"
  ) +
  theme_minimal(base_size=11) +
  theme(
    legend.position = "none",
    panel.grid = element_line(linewidth = 0.3, linetype = 2),
    plot.margin = margin(1,3,1,1, unit = "pt"),
    axis.title = element_text(size=rel(1.1)),
    axis.text = element_text(size=rel(1))
  )
ggsave( plot = pl, filename = "mort_amp.pdf", path = FIG_DIR,
        width = 7, height = 4.5 )

# =============================================================================
# Figure: mort_template_ref_comp_tempo.pdf
# Purpose: Three-panel figure showing (top) the latent template tau(t),
#          (middle) females' resolved component tempo eta_female(t), and
#          (bottom) males' resolved component tempo eta_male(t).  Each panel
#          shows density curves colored by year; the dashed curve is the
#          reference point omega_0.
# =============================================================================

# Flatten the component tempo list (one matrix per sex) into a tidy frame
Ldens <- ldmfit$res_eta
df_eta <- map_dfr(names(Ldens$qout), function(sex) {
  map_dfr( seq_along(dat$obsGrid), function (i) {
    data.frame(
      sex = sex,
      Year = dat$obsGrid[i],
      x = Ldens$dSup[[sex]][i,],
      y = Ldens$dout[[sex]][i,]
    )
  })
}) %>%
  mutate( sex = sex %>% recode( female = "Female", male = "Male" ) )

# Flatten the latent template (one matrix) into a tidy frame
Ldens <- ldmfit$res_tau
df_tau <- map_dfr( seq_along(dat$obsGrid), function (i) {
  data.frame(
    Year = dat$obsGrid[i],
    x = Ldens$dSup[i,],
    y = Ldens$dout[i,]
  )
})

# Reference point omega_0 (time-invariant; drawn as a dashed curve in all panels)
df_omega0 <- data.frame(
  x = ldmfit$res_omega0$dSup,
  y = ldmfit$res_omega0$dout
)

# Combine template and component tempo rows; set factor order for facet_grid
df <- bind_rows(
  df_tau %>% mutate(sex = "Template", .before = 1),
  df_eta %>% mutate(sex = paste0(sex,"s' tempo") ),
  .id = NULL
)
df <- df %>%
  mutate( sex = factor( sex, levels = c("Template", "Females' tempo", "Males' tempo") ) )

pl <- ggplot(df) +
  facet_grid(sex ~ . ) +
  geom_line( aes(x = x, y = y, color = Year, group = Year), linewidth = 0.3, alpha = 0.9) +
  col_scale_small +
  geom_line( aes( x = x, y = y ), linewidth = 0.9, linetype = "dashed",
             data = df_omega0 %>% mutate( sex = factor("Template",levels = c("Template", "Females' tempo", "Males' tempo")) ) ) +
  labs(x = "Age-at-death (year)", y = "Density", color = "Year") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    legend.title = element_text(size=rel(1)),
    legend.text = element_text(size=rel(.9)),
    strip.text = element_text(size=rel(1.1)),
    panel.grid = element_line(linewidth = 0.3, linetype = 2),
    plot.margin = margin(1,3,1,1, unit = "pt"),
    axis.title = element_text(size=rel(1.1)),
    axis.text = element_text(size=rel(1))
  )
ggsave( plot = pl, filename = "mort_template_ref_comp_tempo.pdf", path = FIG_DIR,
        width = 5.5, height = 6 )

# =============================================================================
# Figure: mort_comp_wf.pdf
# Purpose: Component-level deformation functions Psi_female(t) and Psi_male(t)
#          plotted against the identity (dashed).  Deviations from the identity
#          reflect the systematic timing shift of each sex relative to the
#          latent template.
# =============================================================================

# Flatten the per-component deformation list into a tidy frame
Ly <- ldmfit$res_Psi$H
df <- map_dfr(names(Ly), function(sex) {
  data.frame(
    sex = sex,
    x = dat$obsGrid,
    y = Ly[[sex]]
  )
}) %>%
  mutate( sex = sex %>% recode( female = "Female", male = "Male" ) )

pl <- ggplot(df, aes(x = x, y = y, color = sex, group = sex)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.7, linetype = "dashed", color = "gray75") +
  geom_line(linewidth = 1, alpha = 0.9) +
  labs(x = "Year", y = "Year", color = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    legend.title = element_text(size=rel(1.1)),
    legend.text = element_text(size=rel(1)),
    strip.text = element_text(size=rel(1.1)),
    panel.grid = element_line(linewidth = 0.3, linetype = 2),
    plot.margin = margin(1,3,1,1, unit = "pt"),
    axis.title = element_text(size=rel(1)),
    axis.text = element_text(size=rel(.9))
  )
ggsave( plot = pl, filename = "mort_comp_wf.pdf", path = FIG_DIR,
        width = 2.3, height = 2.6 )

# =============================================================================
# Figure: mort_pxc_align.pdf
# Purpose: Population-level cross-component alignment map from females to males,
#          pxcalign_f2m(t) = Psi_female^{-1}(Psi_male(t)).  Values above the
#          identity mean males are systematically ahead of females in time.
# =============================================================================

df <- data.frame(
  x = ldmfit$obsGrid,
  y = ldmfit$pxcalign_f2m
)

pl <- ggplot(df, aes(x = x, y = y)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.7, linetype = "dashed", color = "gray75") +
  geom_line(linewidth = 1, alpha = 0.9) +
  labs(x = "Year", y = "Year") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    legend.title = element_text(size=rel(1.1)),
    legend.text = element_text(size=rel(1)),
    strip.text = element_text(size=rel(1.1)),
    panel.grid = element_line(linewidth = 0.3, linetype = 2),
    plot.margin = margin(1,3,1,1, unit = "pt"),
    axis.title = element_text(size=rel(1)),
    axis.text = element_text(size=rel(.9))
  )
ggsave( plot = pl, filename = "mort_pxc_align.pdf", path = FIG_DIR,
        width = 2.3, height = 2.05 )

# =============================================================================
# Figure: mort_subj_wf_per_country.pdf
# Purpose: Subject-level deformation functions H_i(t) for all 34 countries,
#          one panel per country, colored by 1988 GDP per capita.  Countries
#          with higher GDP (warmer color) tend to show stronger timing shifts.
# =============================================================================

# Collect H_i(t) for each country and attach the GDP covariate
dfsw <- map_dfr( seq_along(dat$subj_labs), function (i) {
  data.frame(
    country = dat$subj_labs[i],
    x = ldmfit$obsGrid,
    y = ldmfit$res_H$H[i,]
  )
}) %>%
  recode_country() %>%
  left_join( y = gdp_pc %>% dplyr::select(country,gdp_pc), by = 'country' )

pl <- ggplot(dfsw, aes(x = x, y = y, color = gdp_pc, group = country)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.6, linetype = "dashed", color = "gray75") +
  geom_line(linewidth = 0.8) +
  col_scale_rev_small +
  facet_wrap(~ country, ncol = 6) +
  labs(x = "Year", y = "Year", color = "GDP per capita in 1988") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    legend.title = element_text(size=rel(1)),
    legend.text = element_text(size=rel(.9)),
    strip.text = element_text(size=rel(1.1)),
    panel.grid = element_line(linewidth = 0.3, linetype = 2),
    plot.margin = margin(1,3,1,1, unit = "pt"),
    axis.title = element_text(size=rel(1.1)),
    axis.text = element_text(size=rel(.9))
  )
ggsave( plot = pl, filename = "mort_subj_wf_per_country.pdf", path = FIG_DIR,
        width = 7, height = 8.9 )

# =============================================================================
# Figure: mort_sxc_align_per_country.pdf
# Purpose: Subject-level cross-component alignment maps sxcalign_f2m[i,](t)
#          for all 34 countries, one panel per country.  Each curve shows the
#          country-specific female-to-male timing alignment; curves below the
#          identity indicate that females lead males in longevity timing.
# =============================================================================

# Collect sxcalign_f2m[i,] for each country and attach GDP
dfa <- map_dfr( seq_along(dat$subj_labs), function (i) {
  data.frame(
    country = dat$subj_labs[i],
    x = ldmfit$obsGrid,
    y = ldmfit$sxcalign_f2m[i,]
  )
}) %>%
  recode_country() %>%
  left_join( y = gdp_pc %>% dplyr::select(country,gdp_pc), by = 'country' )

pl <- ggplot(dfa, aes(x = x, y = y, color = gdp_pc, group = country)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.6, linetype = "dashed", color = "gray75") +
  geom_line(linewidth = 0.8) +
  col_scale_rev_small +
  facet_wrap(~ country, ncol = 6) +
  labs(x = "Year", y = "Year", color = "GDP per capita in 1988") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    legend.title = element_text(size=rel(1)),
    legend.text = element_text(size=rel(.9)),
    strip.text = element_text(size=rel(1.1)),
    panel.grid = element_line(linewidth = 0.3, linetype = 2),
    plot.margin = margin(1,3,1,1, unit = "pt"),
    axis.title = element_text(size=rel(1.1)),
    axis.text = element_text(size=rel(.9))
  )
ggsave( plot = pl, filename = "mort_sxc_align_per_country.pdf", path = FIG_DIR,
        width = 7, height = 8.9 )

# =============================================================================
# Figure: mort_subj_wf_sxc_align_together.pdf
# Purpose: Two-panel summary: subject-level deformation functions (left) and
#          subject-level cross-component alignment maps (right), all countries
#          overlaid in one panel each, colored by GDP.
# =============================================================================

# Stack the two data frames and use the list ID as the facet variable
df_comb <- list(
  f1 = dfsw,  # subject-level deformation functions
  f2 = dfa    # subject-level cross-component alignment maps
) %>% bind_rows(.id = 'fctn')

pl <- ggplot(df_comb, aes(x = x, y = y, color = gdp_pc, group = country)) +
  facet_wrap(~ fctn) +
  geom_line(linewidth = 0.8) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.6, linetype = "dashed", color = "gray55") +
  col_scale_rev_small +
  labs(x = "Year", y = "Year", color = "GDP per capita in 1988") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    legend.title = element_text(size=rel(1)),
    legend.text = element_text(size=rel(.9)),
    strip.text = element_blank(),   # panel labels suppressed; explained in caption
    panel.grid = element_line(linewidth = 0.3, linetype = 2),
    panel.spacing = unit(25, "pt"),
    plot.margin = margin(1,3,1,1, unit = "pt"),
    axis.title = element_text(size=rel(1)),
    axis.text = element_text(size=rel(.9))
  )
ggsave( plot = pl, filename = "mort_subj_wf_sxc_align_together.pdf", path = FIG_DIR,
        width = 5, height = 3.25 )

# =============================================================================
# Figure: mort_subj_nti.pdf
# Purpose: Subject-level net tempo index (sNTI) vs. 1988 GDP per capita.
#          sNTI_i = integral of (H_i(t) - t) dt / range(obsGrid)^2; a positive
#          value means the country's overall mortality timing is ahead of the
#          template, a negative value means it lags behind.
# =============================================================================

# Compute sNTI for each country via trapezoidal integration
df <- map_dfr( seq_along(dat$subj_labs), function (i) {
  data.frame(
    country = dat$subj_labs[i],
    int = pracma::trapz( x = ldmfit$obsGrid, y = ldmfit$res_H$H[i,] - ldmfit$obsGrid )/diff(range(ldmfit$obsGrid))^2
  )
}) %>%
  recode_country() %>%
  left_join( y = gdp_pc %>% dplyr::select(country,gdp_pc), by = 'country' )

pl <- ggplot( df ) +
  geom_hline( yintercept = 0, linewidth = 0.8, linetype = "dashed", color = "gray75") +
  geom_point(aes( x = gdp_pc, y = int ) ) +
  geom_text_repel(
    aes( x = gdp_pc, y = int, label = country ),
    size = 5,
    box.padding = 0.3,
    point.padding = 0.1,
    segment.color = "gray65",
    force = 30,
    force_pull = 0.05,
    min.segment.length = unit(0.1, "lines"),
    max.overlaps = Inf
  ) +
  labs(
    x = "GDP per capita in 1988",
    y = "sNTI"
  ) +
  theme_minimal(base_size=11) +
  theme(
    legend.position = "none",
    panel.grid = element_line(linewidth = 0.3, linetype = 2),
    plot.margin = margin(3,1,1,1, unit = "pt"),
    axis.title = element_text(size=rel(1.3)),
    axis.text = element_text(size=rel(1))
  )
ggsave( plot = pl, filename = "mort_subj_nti.pdf", path = FIG_DIR,
        width = 6, height = 5.5 )

# =============================================================================
# Figure: mort_subj_rti.pdf
# Purpose: Negative subject-level relative timing index (sRTI) of males to
#          females vs. 1988 GDP per capita.
#          sRTI_i = integral of (sxcalign_f2m[i,](t) - t) dt / range^2;
#          the sign is flipped so that a larger value means males age faster
#          relative to females (more intuitive display direction).
# =============================================================================

# Compute sRTI for each country and flip the sign for display
df <- map_dfr( seq_along(dat$subj_labs), function (i) {
  data.frame(
    country = dat$subj_labs[i],
    int = pracma::trapz( x = ldmfit$obsGrid, y = ldmfit$sxcalign_f2m[i,] - ldmfit$obsGrid )/diff(range(ldmfit$obsGrid))^2
  )
}) %>%
  recode_country() %>%
  left_join( y = gdp_pc %>% dplyr::select(country,gdp_pc), by = 'country' ) %>%
  mutate( int = -int )  # flip sign: positive = males age faster than females

pl <- ggplot( df ) +
  geom_point( aes( x = gdp_pc, y = int ) ) +
  geom_text_repel(
    aes( x = gdp_pc, y = int, label = country ),
    size = 5,
    box.padding = 0.3,
    point.padding = 0.1,
    segment.color = "gray65",
    force = 20,
    force_pull = 0.05,
    min.segment.length = unit(0.2, "lines"),
    max.overlaps = Inf
  ) +
  labs(
    x = "GDP per capita in 1988",
    y = "Negative sRTI of males to females"
  ) +
  theme_minimal(base_size=11) +
  theme(
    legend.position = "none",
    panel.grid = element_line(linewidth = 0.3, linetype = 2),
    plot.margin = margin(3,1,1,1, unit = "pt"),
    axis.title = element_text(size=rel(1.3)),
    axis.text = element_text(size=rel(1))
  )
ggsave( plot = pl, filename = "mort_subj_rti.pdf", path = FIG_DIR,
        width = 6, height = 5.5 )
