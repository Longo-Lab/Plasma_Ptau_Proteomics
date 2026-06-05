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
## Date Created: 2025-09-26
##
## 
## 
##
## ---------------------------
##
## Notes:
##   
##   
##
##  
##
## ---------------------------

## load up the packages we will need: 
library(ggplot2)
library(dplyr)
library(ggpubr)
library(Cairo)
library(data.table)

## Windows- set up working directory ------------------
proj.dir <- file.path(
   "C:", "Users", "audweber", "Box", "plasma_ptau_proteomics" 
)
setwd(file.path(proj.dir, "gnpc"))

##Load biomarker and proteome csv 
prote <- fread(
  file.path(
    "C:", 
    "Users", 
    "audweber", 
    "Box", 
    "plasma_ptau_proteomics", 
    "dataProt_SS-205063.merged.QC.csv"
  )
)
# prote <- fread(
#   file.path(
#     "../dataProt_SS-205063.merged.QC.csv"
#   )
# )

# Find proteins and biomarkers --------- 
cols <- names(prote)[names(prote) %like% "GFAP"]
##GFAP
GFAP <- grep("GFAP", names(prote), ignore.case = TRUE)
# column names
gfap_cols <- names(prote)[GFAP]

##NFL
nfl <- grep("NEFL", names(prote), ignore.case = TRUE)
# column names
gfap_cols <- names(prote)[NEFL]

#scatterplots without log transformation 

#NfL correlation 
Nfl <- ggplot(prote, aes(x=PlasmaNfL, y=NEFL.10082.251.3)) +
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
    x = "Plasma NfL",
    y = "NfL"
  ) + 
  stat_cor(method = "pearson", label.x = 2, label.y = 2.5, size=3, face="bold", colour = "black") +
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
Nfl

#GFAP1 correlation 
GFAP1 <- ggplot(prote, aes(x=PlasmaGFAP, y=GFAP.20126.19.3)) +
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
    x = "Plasma GFAP",
    y = "GFAP1"
  ) + 
  stat_cor(method = "pearson", label.x = 0, label.y = 1000, size=3, face="bold", colour = "black") +
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
GFAP1

#Logtransformed correlation for GFAP2

GFAP2 <- ggplot(prote, aes(x=PlasmaGFAP, y=GFAP.3034.1.2)) +
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
    x = "Plasma GFAP",
    y = "GFAP2"
  ) + 
  stat_cor(method = "pearson", label.x = 0, label.y = 2000, size=3, face="bold", colour = "black") +
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
GFAP2

##Log10 transform data for correlations
#Determine what columns you need
cols <- c("PlasmaNfL", "NEFL.10082.251.3", "PlasmaGFAP", "GFAP.20126.19.3", "GFAP.3034.1.2")
#Select them and log transform them in a new dataframe
cols_log10 <- prote %>%
  select(all_of(cols)) %>%
  mutate(across(everything(), ~ log10(.x + 1)))

#Graph each log transformed dataset

#NfL correlation 
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


##run correlation long10 transformed data 

cor.test(cols_log10$PlasmaNfL, cols_log10$NEFL.10082.251.3, method = "pearson", use = "complete.obs")
cor.test(cols_log10$PlasmaGFAP, cols_log10$GFAP.20126.19.3, method = "pearson", use = "complete.obs")
cor.test(cols_log10$PlasmaGFAP, cols_log10$GFAP.3034.1.2, method = "pearson", use = "complete.obs")

#Print all plots in a pdf 


cairo_pdf("AW_09.29.2025_biomarker_correlations.pdf", width = 10, height = 4)
print(Nfl | log10_nfl)  # page 1
print(GFAP1 | log10_GFAP1)  # page 2
print(GFAP2 | log10_GFAP2)  # page 3
dev.off()




