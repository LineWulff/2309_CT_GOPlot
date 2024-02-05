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

#### Biplot ####
Panc <- Pancreas[,c("g_symbol","log2FoldChange","padj")]
Panc$tissue <- "Pancreas"
Panc <- Panc[!is.na(Panc$log2FoldChange),]
Panc <- Panc[!is.na(Panc$padj),]
Panc <- Panc[!is.na(Panc$g_symbol),]
head(Panc)
dim(Panc)
Panc <- Panc[Panc$padj<0.05,]
dim(Panc)
colnames(Panc)[2] <- "log2FC_Panc"
Panc <- Panc[!Panc$g_symbol %in% Panc$g_symbol[duplicated(Panc$g_symbol)],]


SI <- SI[,c("g_symbol","log2FoldChange","padj")]
SI$tissue <- "SI"
SI <- SI[!is.na(SI$log2FoldChange),]
SI <- SI[!is.na(SI$padj),]
SI <- SI[!is.na(SI$g_symbol),]
head(SI)
dim(SI)
SI <- SI[SI$padj<0.05,]
dim(SI)
colnames(SI)[2] <- "log2FC_SI"
#SI <- SI[!SI$g_symbol %in% SI$g_symbol[duplicated(SI$g_symbol)],]

# make df from SI and add overlapping genes from Pancreas
biplot_df <- SI
biplot_df$log2FC_Panc <- 0
biplot_df[biplot_df$g_symbol %in% intersect(Panc$g_symbol, SI$g_symbol),]$log2FC_Panc <- Panc[Panc$g_symbol %in% SI$g_symbol,]$log2FC_Panc
# now add pancreas genes not in SI 
biplot_df2 <- Panc[!Panc$g_symbol %in% SI$g_symbol,]
biplot_df2$log2FC_SI <- 0
# concatenate together
biplot_df <- rbind(biplot_df,biplot_df2)
# make colour scheme
biplot_df$shared <- "Not shared"
biplot_df[biplot_df$log2FC_SI>0 & biplot_df$log2FC_Panc>0,]$shared <- "SI up, PC up"
biplot_df[biplot_df$log2FC_SI<0 & biplot_df$log2FC_Panc<0,]$shared <- "SI down, PC down"
biplot_df[biplot_df$log2FC_SI<0 & biplot_df$log2FC_Panc>0,]$shared <- "SI down, PC up"
biplot_df[biplot_df$log2FC_SI>0 & biplot_df$log2FC_Panc<0,]$shared <- "SI up, PC down"

top10 <- biplot_df %>% arrange(log2FC_Panc) %>% top_n(n=10)
top10 <- biplot_df %>% group_by(shared) %>% top_n(n=10)
biplot_df[order(biplot_df$log2FC_Panc),]

cols_shared <- c("SI up, PC up"="#35978f","SI down, PC down"="#dfc27d",
                 "SI down, PC up"="#c2a5cf","SI up, PC down"="#74a9cf",
                 "Not shared"="lightgrey")

pdf(paste(dir,"/output/",dato,"_Biplot_DEGs_PancreasSI_v2.pdf",sep=""), height = 4, width = 5)
ggplot(biplot_df[biplot_df$shared!="Not shared",], 
       aes(x = log2FC_SI, y = log2FC_Panc, colour = shared, ))+ #label = g_symbol))+
  geom_point()+
 # geom_text(data = biplot_df[1:10,], colour = "black")+
  scale_colour_manual(values = cols_shared)+
  labs(colour = "Shared DEGs", x = "log2FC SI", y = "log2FC Pancreas")+
  theme_minimal()
dev.off()

# Save csv of biplot wo. "not shared" for Carolyn
write.csv(biplot_df[biplot_df$shared!="Not shared",], 
          file = paste(dir,"/output/",dato,"_biplotdata_DEGsPancreasSI.csv", sep=""))

write.csv(biplot_df[biplot_df$shared=="Not shared",], 
          file = paste(dir,"/output/",dato,"_biplotdata_DEGsPancreasSI_notshared.csv", sep=""))
