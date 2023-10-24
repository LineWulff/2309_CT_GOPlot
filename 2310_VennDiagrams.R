## GO enrichment Plots for Carolyn Thompson
## VErsion 2 - Venn diagrams 
## Author: Line Wulff
## Date created: 2023.10.24

#### Libraries ####
library(ggplot2)
library(tidyverse)
library(viridis)
library(stringr)
library(VennDiagram)


#### General varaibles ####
dir <- "/Users/linewulff/Documents/work/projects/2309_CT_GOPlot"
setwd(dir)
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)

#### Read in data from CT ####
Pancreas <- read.csv(paste(dir,"inputdata","Panc_DEG_100+.csv", sep="/"))
head(Pancreas)

SI <- read.csv(paste(dir,"inputdata","SI_DEG_100+.csv", sep="/"))
head(SI)

### DEGs per data set ####
SI_up <- SI[SI$log2FoldChange>0.26303 & SI$padj<0.05,]$g_symbol
SI_do <- SI[SI$log2FoldChange<(-0.26303) & SI$padj<0.05,]$g_symbol

PC_up <- Pancreas[Pancreas$log2FoldChange>0.26303 & Pancreas$padj<0.05,]$g_symbol
PC_do <- Pancreas[Pancreas$log2FoldChange<(-0.26303) & Pancreas$padj<0.05,]$g_symbol

#### Venn diagrams ####
# All
venn.diagram(
  x = list(PC_up, PC_do, SI_up, SI_do),
  category.names = c("Pancreas up" ,"Pancreas down" , "SI up", "SI down"),
  output=F,
  
  ## image
  # Output features
  imagetype="tiff",
  filename = 'SIPancreasAllDEGs_VennDiagramm.tiff',
  height = 300, width = 300, resolution = 150)

# Up 
venn.diagram(
  x = list(PC_up, SI_up),
  category.names = c("Pancreas up", "SI up"),
  output=F,
  
  ## image
  # Output features
  imagetype="tiff",
  filename = 'SIPancreasUpDEGs_VennDiagramm.tiff',
  height = 300, width = 300, resolution = 150)

# Down
venn.diagram(
  x = list(PC_down, SI_down),
  category.names = c("Pancreas down", "SI down"),
  output=F,
  
  ## image
  # Output features
  imagetype="tiff",
  filename = 'SIPancreasUpDEGs_VennDiagramm.tiff',
  height = 300, width = 300, resolution = 150)

#### Overlaps ####
com_up <- intersect(SI_up, PC_up)
com_down <- intersect(SI_do, PC_do)

comm <- vectorlist(list(UpReg = com_up, DownReg = com_down))

write.csv(comm, file = paste(dato,"overlappingDEGs_PancreasSI.csv",sep = "_"))
