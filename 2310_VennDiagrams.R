## GO enrichment Plots for Carolyn Thompson
## Venn diagrams and overlapping DEGs between the two organs
## Author: Line Wulff
## Date created: 2023.10.24

#### Libraries ####
library(ggplot2)
library(tidyverse)
library(viridis)
library(stringr)
library(VennDiagram)

#### General variables ####
dir <- "/Users/linewulff/Documents/work/projects/2309_CT_GOPlot"
setwd(dir)
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)

#### Read in data from CT ####
Pancreas <- read.csv(paste(dir,"inputdata","Panc_all_genes.csv", sep="/"))
head(Pancreas)

SI <- read.csv(paste(dir,"inputdata","SI_all_genes.csv", sep="/"))
head(SI)

### DEGs per data set ####
# CT thresholds: basemean > 100, padj < 0.05, abs(logFC) >0.26303
# Gene symbols
SI_up <- SI[SI$log2FoldChange>0.26303 & SI$padj<0.05 & SI$baseMean>100,]$g_symbol
SI_up <- SI_up[!is.na(SI_up)]
SI_do <- SI[SI$log2FoldChange<(-0.26303) & SI$padj<0.05 & SI$baseMean>100 ,]$g_symbol
SI_do <- SI_do[!is.na(SI_do)]

PC_up <- Pancreas[Pancreas$log2FoldChange>0.26303 & Pancreas$padj<0.05 & Pancreas$baseMean>100,]$g_symbol
PC_up <- PC_up[!is.na(PC_up)]
PC_do <- Pancreas[Pancreas$log2FoldChange<(-0.26303) & Pancreas$padj<0.05 & Pancreas$baseMean>100,]$g_symbol
PC_do <- PC_do[!is.na(PC_do)]

# Gene ENSEMBL IDs
# SI_up <- SI[SI$log2FoldChange>0.26303 & SI$padj<0.05 & SI$baseMean>100,]$g_id
# SI_up <- SI_up[!is.na(SI_up)]
# SI_do <- SI[SI$log2FoldChange<(-0.26303) & SI$padj<0.05 & SI$baseMean>100 ,]$g_id
# SI_do <- SI_do[!is.na(SI_do)]
# 
# PC_up <- Pancreas[Pancreas$log2FoldChange>0.26303 & Pancreas$padj<0.05 & Pancreas$baseMean>100,]$g_id
# PC_up <- PC_up[!is.na(PC_up)]
# PC_do <- Pancreas[Pancreas$log2FoldChange<(-0.26303) & Pancreas$padj<0.05 & Pancreas$baseMean>100,]$g_id
# PC_do <- PC_do[!is.na(PC_do)]

#### Venn diagrams ####
# All
venn.diagram(
  x = list(PC_up, PC_do, SI_up, SI_do),
  category.names = c("PC up" ,"PC down" , "SI up", "SI down"),
  output=F,
  
  ## image
  # Output features
  imagetype="tiff",
  filename = 'SIPancreasAllDEGs_VennDiagramm.tiff',
  height = 300, width = 300, resolution = 150)

# Up 
venn.diagram(
  x = list(PC_up, SI_up),
  category.names = c("PC up", "SI up"),
  output=F,
  
  ## image
  # Output features
  imagetype="tiff",
  filename = 'SIPancreasUpDEGs_VennDiagramm.tiff',
  height = 300, width = 300, resolution = 150)

# Down
venn.diagram(
  x = list(PC_do, SI_do),
  category.names = c("PC down", "SI down"),
  output=F,
  
  ## image
  # Output features
  imagetype="tiff",
  filename = 'SIPancreasDownDEGs_VennDiagramm.tiff',
  height = 300, width = 300, resolution = 150)

#### Overlaps ####
com_up <- intersect(SI_up, PC_up)
com_down <- intersect(SI_do, PC_do)

comm <- vectorlist(list(UpReg = com_up, DownReg = com_down))
# ENSEMBL IDs translated
# comm <- vectorlist(list(UpReg_gID = com_up, UpReg_symb =  Pancreas[Pancreas$g_id %in% com_up,]$g_symbol, DownReg_gID = com_down, DownReg_symb = Pancreas[Pancreas$g_id %in% com_down,]$g_symbol))


write.csv(comm, file = paste(dato,"overlappingDEGs_PancreasSI.csv",sep = "_"))
