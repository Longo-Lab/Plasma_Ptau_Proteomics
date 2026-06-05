library(knitr)
library(dplyr)
library(tidyr)

kable(head(metrics), format = "markdown")
df <- read_excel("your_file.xlsx")


#glmnet markdown

setwd("C:/Users/audweber/Box/plasma_ptau_proteomics/prediction_models/biodomains/glmnet/INT_biomarkers")

df <- read.csv("glmnet_outcomes_by_variant_INT.csv", stringsAsFactors = FALSE)

glmnet_wide <- df %>%
  mutate(rmse_fmt = paste0(outcome, ": ", round(rse_oof, 4))) %>%
  group_by(variant) %>%
  summarise(summary = paste(rse_fmt, collapse = " ")) %>%
  pivot_wider(names_from = variant, values_from = summary) %>%
  mutate(model = "individual glmnets") %>%
  select(model, everything())


metrics <- c("rmse_oob", "rse_oob", "mae", "mape")

# Loop through each metric and create + print markdown table
for (metric in metrics) {
  
  # Dynamically build the metric-specific column
  df_wide <- all_data %>%
    mutate(metric_fmt = paste0(outcome, ": ", round(.data[[metric]], 4))) %>%
    group_by(variant) %>%
    summarise(summary = paste(metric_fmt, collapse = " ")) %>%
    pivot_wider(names_from = variant, values_from = summary) %>%
    mutate(model = paste("individual glmnets —", toupper(metric))) %>%
    select(model, everything())
  
  # Print table title
  cat("\n###", toupper(metric), "\n")
  
  # Print the markdown table
  print(kable(df_wide, format = "markdown"))
}


#rf markdown

setwd("C:/Users/audweber/Box/plasma_ptau_proteomics/prediction_models/biodomains/rf/INT_biomarkers")

file_list <- list.files(path = "INT_biomarkers", pattern = "\\.csv$", full.names = TRUE)


metrics <- c("rmse_oob", "rse_oob", "mae", "mape")

# Loop through metrics
for (metric in metrics) {
  
  rf_wide <- all_data %>%
    mutate(metric_fmt = paste0(outcome, ": ", round(.data[[metric]], 4))) %>%
    group_by(variant) %>%
    summarise(summary = paste(metric_fmt, collapse = " ")) %>%
    pivot_wider(names_from = variant, values_from = summary) %>%
    mutate(model = paste("individual rf —", toupper(metric))) %>%
    select(model, everything())
  
  # Print title and markdown table
  cat("\n###", toupper(metric), "\n")
  print(kable(rf_wide, format = "markdown"))
}



#joinet
setwd("C:/Users/audweber/Box/plasma_ptau_proteomics/prediction_models/biodomains/multi_rf/INT_biomarkers")
all_macro <- read_csv("INT_mrf_macro_by_variant.csv")

macro_metrics <- c("macro_rmse", "macro_r2", "macro_rse", "macro_mae", "macro_mape")

# Create a long format version: one row per (metric × variant)
macro_long <- all_macro %>%
  select(variant, all_of(macro_metrics)) %>%
  pivot_longer(cols = all_of(macro_metrics), names_to = "metric", values_to = "value")

# Round for display
macro_long <- macro_long %>%
  mutate(value = round(value, 4))

# Now pivot wider to make the summary table: one row per metric
macro_wide <- macro_long %>%
  pivot_wider(names_from = variant, values_from = value) %>%
  mutate(model = "macro glmnet") %>%
  select(model, everything())

# Print the final Markdown table
kable(macro_wide, format = "markdown")


