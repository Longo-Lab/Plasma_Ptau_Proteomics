library(dplyr)
library(stringr)
library(knitr)
library(ggplot2)

output <- "C:/Users/audweber/Box/plasma_ptau_proteomics/prediction_models/eigengenes/mrf/z_scale/"

df <- read.csv("C:/Users/audweber/Box/plasma_ptau_proteomics/prediction_models/eigengenes/glmnet/mgaussian/glmnet_model_metrics_gmaussian_z_scale_INT.csv")


df %>%
  mutate(
    variant = factor(variant, levels = unique(variant)),
    C_index_num = as.numeric(str_extract(C_index, "^[0-9.]+"))
  ) %>%
  group_by(variant) %>%
  summarise(
    RMSE    = median(RMSE, na.rm = TRUE),
    MAE     = median(MAE, na.rm = TRUE),
    R2      = median(R2, na.rm = TRUE),
    C_index = median(C_index_num, na.rm = TRUE),
    N     = median(N, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  knitr::kable(
    format = "markdown",
    digits = 3
  )

df_corr <- read.csv("C:/Users/audweber/Box/plasma_ptau_proteomics/prediction_models/eigengenes/glmnet/mgaussian/glmnet_model_correlations_eigengenes_mgaussian_z_scale.csv")


p <- ggplot(df_corr, aes(x = variant, y = Pearson_Correlation)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.3) +
  stat_summary(
    fun = median,
    geom = "point",
    size = 3,
    color = "darkgreen"
  ) +
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = round(..y.., 2)),
    vjust = -1,
    color = "darkgreen",
    size = 4.5,
    fontface = "bold"
  ) +
  scale_y_continuous(limits = c(0, 1)) +   # <-- y-axis fixed to 0–1
  theme_classic() +
  labs(
    x = NULL,
    y = "Correlation (actual vs predicted)"
  )

ggsave(file.path(output, "correlation_boxplot_eigengene_mgaussian_z_scale_INT.png"), plot = p, width = 8, height = 6, dpi = 300)



