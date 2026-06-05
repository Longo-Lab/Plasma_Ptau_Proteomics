#
## ---------------------------
##
## Script name: biomarker_proteome_correlations.R
##
## Version: 0.0.1
##
## Purpose of script: Show correlations between plasma biomarkers 
## and proteins
##
## Author: Audrey Weber
##
## Date Created: 2025-09-29
##
## 
## 
##
## ---------------------------
##
## Notes:
##   
## Interested in seeing if diagnosis had an effect on these correlations 
##
##  
##
## ---------------------------

## load up the packages we will need: 
library(ggplot2)
library(dplyr)
library(ggpubr)
library(Cairo)





#Filer diagnosis 

diagnosis <- c("HC")
AD <- prote %>% filter(`Diagnosis_group` %in% diagnosis)

#Log transform data
cols_log10 <- AD %>%
  select(all_of(cols)) %>%
  mutate(across(everything(), ~ log10(.x + 1)))

#Graph
log10_nfl <- ggplot(cols_log10, aes(x=PlasmaNfL, y=NEFL.10082.251.3)) +
  geom_smooth(method = "lm", se = FALSE, color = "#69b3a2", linewidth = 0.5) +
  geom_point( 
    color="#69b3a2",
    fill="#69b3a2",
    shape=21,
    alpha=0.75,
    size=0.5,
    stroke = 1
  ) +
  labs(
    x = expression(bold(log[10](PlasmaNfL))),
    y = expression(bold(log[10](NfL)))
  ) + 
  stat_cor(method = "pearson", label.x = 0.8, label.y = 2, size=3, face="bold", colour = "black") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),    
    axis.line = element_line(),        
    axis.ticks = element_line(),
    axis.title.x = element_text(size=14, face="bold", colour = "black"),    
    axis.title.y = element_text(size=14, face="bold", colour = "black"), 
    axis.text.x = element_text(size=12, face="bold", colour = "black"),
    axis.text.y = element_text(size=12, face="bold", colour = "black"),
    axis.line.x = element_line(color="black"),
    axis.line.y = element_line(color="black"),
  )
log10_nfl

#GFAP1 correlation 
log10_GFAP1 <- ggplot(cols_log10, aes(x=PlasmaGFAP, y=GFAP.20126.19.3)) +
  geom_smooth(method = "lm", se = FALSE, color = "salmon", linewidth = 0.5) +
  geom_point( 
    color="salmon",
    fill="salmon",
    shape=21,
    alpha=0.75,
    size=0.5,
    stroke = 1
  ) +
  labs(
    x = expression(bold(log[10](PlasmaGFAP))),
    y = expression(bold(log[10](GFAP1)))
  ) + 
  stat_cor(method = "pearson", label.x = 1.2, label.y = 2.6, size=3, face="bold", colour = "black") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),    
    axis.line = element_line(),        
    axis.ticks = element_line(),
    axis.title.x = element_text(size=14, face="bold", colour = "black"),    
    axis.title.y = element_text(size=14, face="bold", colour = "black"), 
    axis.text.x = element_text(size=12, face="bold", colour = "black"),
    axis.text.y = element_text(size=12, face="bold", colour = "black"),
    axis.line.x = element_line(color="black"),
    axis.line.y = element_line(color="black")
  )
log10_GFAP1

#Logtransformed correlation for GFAP2

log10_GFAP2 <- ggplot(cols_log10, aes(x=PlasmaGFAP, y=GFAP.3034.1.2)) +
  geom_smooth(method = "lm", se = FALSE, color = "sandybrown", linewidth = 0.5) +
  geom_point( 
    color="sandybrown",
    fill="sandybrown",
    shape=21,
    alpha=0.75,
    size=0.5,
    stroke = 1
  ) +
  labs(
    x = expression(bold(log[10](PlasmaGFAP))),
    y = expression(bold(log[10](GFAP2)))
  ) + 
  stat_cor(method = "pearson", label.x = 1.2, label.y = 2.5, size=3, face="bold", colour = "black") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),    
    axis.line = element_line(),        
    axis.ticks = element_line(),
    axis.title.x = element_text(size=14, face="bold", colour = "black"),    
    axis.title.y = element_text(size=14, face="bold", colour = "black"), 
    axis.text.x = element_text(size=12, face="bold", colour = "black"),
    axis.text.y = element_text(size=12, face="bold", colour = "black"),
    axis.line.x = element_line(color="black"),
    axis.line.y = element_line(color="black")
  )
log10_GFAP2

cairo_pdf("AW_09.29.2025_biomarker_correlations_by_diagnosis_HC.pdf", width = 10, height = 4)
print(log10_nfl)  
print(log10_GFAP1) 
print(log10_GFAP2) 
dev.off()


