## GO enrichment Plots for Carolyn Thompson
## Venn diagrams and overlapping DEGs between the two organs
## Author: Line Wulff
## Date created: 2023.10.24

#### Libraries ####
library(ggplot2)
library(tidyverse)
library(viridis)
library(stringr)
library(STRINGdb)

#### General variables ####
dir <- "/Users/linewulff/Documents/work/projects/2309_CT_GOPlot"
setwd(dir)
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)

#### Read in data from CT ####
Pancreas <- read.csv(paste(dir,"inputdata","Panc_all_genes.csv", sep="/"))
head(Pancreas)

SI <- read.csv(paste(dir,"inputdata","SI_all_genes.csv", sep="/"))
head(SI)

## DEGs per data set
## CT thresholds: basemean > 100, padj < 0.05, abs(logFC) >0.26303
SI_up <- SI[SI$log2FoldChange>0.26303 & SI$padj<0.05 & SI$baseMean>100,]
SI_up <- SI_up[!is.na(SI_up$g_id),]
SI_do <- SI[SI$log2FoldChange<(-0.26303) & SI$padj<0.05 & SI$baseMean>100,]
SI_do <- SI_do[!is.na(SI_do$g_id),]

PC_up <- Pancreas[Pancreas$log2FoldChange>0.26303 & Pancreas$padj<0.05 & Pancreas$baseMean>100,]
PC_up <- PC_up[!is.na(PC_up$g_id),]
PC_do <- Pancreas[Pancreas$log2FoldChange<(-0.26303) & Pancreas$padj<0.05 & Pancreas$baseMean>100,]
PC_do <- PC_do[!is.na(PC_do$g_id),]

#### STRING analysis ###
# create a new STRING_db object
string_db <- STRINGdb$new()

SIup_map <- string_db$map( SI_up, "g_symbol", removeUnmappedRows = TRUE )

hits = SIup_map$STRING_id[1:200]

string_db$plot_network(hits)

# get clusters
clustersList = string_db$get_clusters(SIup_map$STRING_id[1:600])

# plot first 4 clusters
par(mfrow=c(2,2))
for(i in seq(1:4)){
  string_db$plot_network(clustersList[[i]])
}
