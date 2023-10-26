## GO enrichment Plots for CT
## Creation of radarplot based selected terms by CT.
## Author: Line Wulff
## Date created: 2023.10.26

#### Libraries ####
library(ggplot2)
library(tidyverse)
library(viridis)
library(stringr)

#### General varaibles ####
dir <- "/Users/linewulff/Documents/work/projects/2309_CT_GOPlot"
setwd(dir)
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)

#### Read in data from CT ####
GO_df <- read.csv(paste(dir,"inputdata","Metabolismpathways_quicklook_sh2.csv", sep="/"))
GO_df <- GO_df[1:20,c("Pathway","Strength","Ratio","p.adj")]
head(GO_df)

GO_df$Pathway <- factor(GO_df$Pathway, levels = c(GO_df[order(GO_df$Strength),]$Pathway))
