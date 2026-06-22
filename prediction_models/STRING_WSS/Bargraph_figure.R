# =============================================================================
# Joint Model RMSE Bar Plot
# x-axis = covariate set (variant), bars coloured by model type
# =============================================================================

library(tidyverse)
library(ggplot2)
library(ggsci)

# ---- USER SETTINGS ----------------------------------------------------------

csv_dir <- "all_model/WSS_CSVs"    # <-- folder containing your CSVs
out_dir  <- "all_model"   # <-- folder to save the plot

# ---- MODEL NAME MAPPING (filename stem -> display name) ---------------------
# Edit the right-hand side if your filenames differ

model_name_map <- c(
  "Glmnet Model Metrics C Index Eigengenes Int Updated" = "Elastic Net",
  "Joinet Model Metrics C Index Eigengenes Int 2"       = "Joinet Multivariate Elastic Net",
  "Rf Model Metrics C Index Eigengenes Int 2"           = "Random Forest",
  "Rfsrc Model Metrics 2"                               = "Multivariate Random Forest"
)

# ---- VARIANT ORDER & DISPLAY LABELS -----------------------------------------

variant_levels <- c(
  "ALL",
  "Age_Sex_HybridNormControl_PCs",
  "Age_Sex_HybridNormControl",
  "Biomarkers_Age_Sex",
  "Biomarkers_Only"
)

variant_labels <- c(
  "ALL"                           = "All",
  "Age_Sex_HybridNormControl_PCs" = "Age + Gender + PCs",
  "Age_Sex_HybridNormControl"     = "HybridNormControl\n+ Gender + Sex",
  "Biomarkers_Age_Sex"            = "Age + Gender",
  "Biomarkers_Only"               = "Biomarkers Only"
)

# ---- LOAD & COMBINE CSVs ----------------------------------------------------

csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) stop("No CSV files found in: ", csv_dir)
message("Found ", length(csv_files), " CSV file(s).")

df_all <- map_dfr(csv_files, function(f) {
  stem <- tools::file_path_sans_ext(basename(f))
  # Match model name by detecting keywords in the filename (case-insensitive)
  model_label <- case_when(
    str_detect(stem, regex("rfsrc", ignore_case = TRUE))  ~ "Multivariate Random Forest",
    str_detect(stem, regex("joinet", ignore_case = TRUE)) ~ "Joinet Multivariate Elastic Net",
    str_detect(stem, regex("glmnet", ignore_case = TRUE)) ~ "Elastic Net",
    str_detect(stem, regex("^rf_|_rf_|_rf$|^rf$", ignore_case = TRUE)) ~ "Random Forest",
    TRUE ~ stem   # fallback: use raw filename stem
  )
  read_csv(f, show_col_types = FALSE) |>
    mutate(model = model_label)
})

# ---- PREPARE DATA -----------------------------------------------------------

model_order <- c("Elastic Net", "Joinet Multivariate Elastic Net",
                 "Random Forest", "Multivariate Random Forest")

df_plot <- df_all |>
  filter(variant %in% variant_levels) |>
  mutate(
    RMSE    = as.numeric(RMSE),
    variant = factor(variant_labels[variant],
                     levels = unname(variant_labels)),
    model   = factor(model, levels = model_order)
  ) |>
  separate(col = C_index, into = c("c.index", NA), sep = " ± ", convert = TRUE, remove = FALSE)

# ---- COLOURS ----------------------------------------------------------------

npg_cols <- pal_npg("nrc")(4)
model_colours <- c(
  "Elastic Net"          = npg_cols[1],  # #E64B35 deep red
  "Joinet Multivariate Elastic Net"                     = npg_cols[2],  # #4DBBD5 teal
  "Random Forest"              = npg_cols[3],  # #00A087 green
  "Multivariate Random Forest" = npg_cols[4]   # #3C5488 navy
)

# ---- PLOT FUNCTION ----------------------------------------------------------

make_metric_plot <- function(data, metric, title, subtitle = NULL, y_label,
                             y_limits = NULL, y_breaks = waiver()) {
  
  dodge <- position_dodge(width = 0.75)
  
  p <- ggplot(data,
              aes(x = variant, y = .data[[metric]],
                  fill = model, colour = model,
                  group = model)) +
    # Semi-transparent bar to mean
    stat_summary(
      fun      = mean,
      geom     = "bar",
      position = dodge,
      width    = 0.65,
      alpha    = 0.3,
      colour   = NA
    ) +
    # Individual eigengene dots
    geom_point(
      position = position_jitterdodge(
        jitter.width = 0.1,
        dodge.width  = 0.75
      ),
      size  = 1.8,
      alpha = 0.8,
      shape = 16
    ) +
    # Median line
    stat_summary(
      fun       = median,
      geom      = "crossbar",
      position  = dodge,
      width     = 0.65,
      linewidth = 2.0,
      fatten    = 0,
      colour    = "grey15",
      fill      = NA
    ) +
    scale_fill_manual(values   = model_colours, name = "Model") +
    scale_colour_manual(values = model_colours, name = "Model") +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.03)),
      breaks = y_breaks
    ) +
    { if (!is.null(y_limits)) coord_cartesian(ylim = y_limits) } +
    labs(
      title    = title,
      subtitle = subtitle,
      x        = NULL,
      y        = y_label
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x      = element_text(size = 9, lineheight = 1.2,
                                      colour = "black", face = "bold"),
      axis.text.y      = element_text(size = 10, colour = "black", face = "bold"),
      axis.title.y     = element_text(size = 11, face = "bold",
                                      margin = margin(r = 8)),
      axis.title.x     = element_text(size = 11, face = "bold"),
      axis.line        = element_line(colour = "black"),
      axis.ticks       = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position  = "right",
      legend.title     = element_text(face = "bold", size = 10),
      legend.text      = element_text(size = 9),
      legend.key.size  = unit(0.5, "cm"),
      plot.title       = element_text(face = "bold", size = 13, colour = "grey10"),
      plot.subtitle    = element_text(size = 9, colour = "grey50"),
      plot.margin      = margin(12, 16, 12, 12)
    )
  
  p
}

# ---- GENERATE ALL THREE PLOTS -----------------------------------------------

# RMSE (lower = better, range 0–0.05)
p_rmse <- make_metric_plot(
  data      = df_plot,
  metric    = "RMSE",
  title     = "RMSE by Covariate Set and Model",
  subtitle  = "Each dot = one module",
  y_label   = "RMSE",
  y_limits  = c(0, 1.5),
  y_breaks  = seq(0, 1.5, by = 0.5)
)
ggsave(file.path(out_dir, "plot_RMSE.pdf"), plot = p_rmse,
       width = 11, height = 6, device = cairo_pdf)
message("Saved: plot_RMSE.pdf")

# MAE (lower = better, auto range)
p_mae <- make_metric_plot(
  data      = df_plot,
  metric    = "MAE",
  title     = "MAE by Covariate Set and Model",
  subtitle  = "Each dot = one module",
  y_label   = "MAE",
  y_limits  = c(0, max(df_plot$MAE, na.rm = TRUE) * 1.05)
)
ggsave(file.path(out_dir, "plot_MAE.pdf"), plot = p_mae,
       width = 11, height = 6, device = cairo_pdf)
message("Saved: plot_MAE.pdf")

# R² (higher = better, auto-scaled to data)
r2_max <- ceiling(max(df_plot$R2, na.rm = TRUE) * 10) / 10  # round up to nearest 0.1
p_r2 <- make_metric_plot(
  data      = df_plot,
  metric    = "R2",
  title     = "R² by Covariate Set and Model",
  subtitle  = "Each dot = one module",
  y_label   = expression(R^2),
  y_limits  = c(0, r2_max),
  y_breaks  = seq(0, r2_max, by = 0.1)
)
ggsave(file.path(out_dir, "plot_R2.pdf"), plot = p_r2,
       width = 11, height = 4, device = cairo_pdf)
message("Saved: plot_R2.pdf")

# C-index (higher = better, auto-scaled to data)
c_max <- ceiling(max(df_plot$c.index, na.rm = TRUE) * 10) / 10  # round up to nearest 0.1
p_c <- make_metric_plot(
  data      = df_plot,
  metric    = "c.index",
  title     = "Concordance index by Covariate Set and Model",
  subtitle  = "Each dot = one module",
  y_label   = "C-index",
  y_limits  = c(0.5, c_max),
  y_breaks  = seq(0, c_max, by = 0.1)
)
ggsave(file.path(out_dir, "plot_c-index.pdf"), plot = p_c,
       width = 11, height = 4, device = cairo_pdf)
message("Saved: plot_c-index.pdf")

# C-index subset (higher = better, auto-scaled to data)
p_c_small <- make_metric_plot(
  data      = df_plot |> filter(variant %in% c("All", "Age + Gender")),
  metric    = "c.index",
  title     = "Concordance index by Covariate Set and Model",
  subtitle  = "Each dot = one module",
  y_label   = "C-index",
  y_limits  = c(0.5, c_max),
  y_breaks  = seq(0, c_max, by = 0.1)
)
ggsave(file.path(out_dir, "plot_c-index_small.pdf"), plot = p_c_small,
       width = 6, height = 4, device = cairo_pdf)
message("Saved: plot_c-index_small.pdf")

message("\nAll done! Plots saved to: ", normalizePath(out_dir))
