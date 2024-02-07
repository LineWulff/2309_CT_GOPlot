## Metabolic data to biplots for Carolyn Thompson
## Author: Line Wulff
## Date created: 2024.02.05

rm(list = ls())

#### Libraries ####
library(ggplot2)
library(tidyverse)
library(viridis)
library(stringr)

#### General varaibles ####
dir <- "/Users/linewulff/Documents/work/projects/2309_CT_GOPlot"
setwd(dir)
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)

#### Read in data ####
# consist of one data sheet
# also a complimenting excel sheets with samples to combine and compare
metb <- read.csv("/Users/linewulff/Documents/work/projects/2309_CT_GOPlot/inputdata/HA107_all combined_FC2.csv", row.names = 1)
metb[1:5,1:5]

#### Groupings and pullings samples from groupings ####
# example
GF <- c("Pancreas_E12_NOD_GF_2w", "Pancreas_E19_NOD_GF_2w")
HA <- c("Pancreas_E12_NOD_HA_2w",	"Pancreas_E19_NOD_HA_2w")

colnames(metb)
GF <- c(colnames(metb)[startsWith(colnames(metb),GF)])
HA <- c(colnames(metb)[startsWith(colnames(metb),HA)])
samp_n <- length(c(HA,GF))

norm_df <- as.data.frame(t(metb[6,c(GF,HA)]))
norm_df$col <- unlist(str_split(rownames(norm_df),pattern = "_"))[seq(4,samp_n*7,7)]
norm_df$ET <- unlist(str_split(rownames(norm_df),pattern = "_"))[seq(2,samp_n*7,7)]

ggplot(norm_df, aes(x = col, y = `2-5-DIHYDROXYBENZOATE`))+
  geom_boxplot()+
  geom_point(aes(colour=ET))+
  theme_classic()

# now as a loop
