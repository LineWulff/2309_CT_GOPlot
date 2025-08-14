## GO enrichment Plots for Carolyn Thompson
## Venn diagrams and overlapping DEGs between the two organs
## Author: Line Wulff
## Date created: 2023.10.24

#### Libraries ####
library(ggplot2)
library(tidyverse)
library(viridis)
library(stringr)
library(ggrepel)
library(scales)
library(ComplexHeatmap)
library(circlize)

#### General variables ####
dir <- "/Users/linewulff/Documents/work/projects/2309_CT_GOPlot"
setwd(dir)
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)
## colouring
myColorRamp <- function(colors, values) {
  v <- (values - min(values))/diff(range(values))
  x <- colorRamp(colors)(v)
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}
mycols <- colorRamp2(c(-2, 0, 2),c("#762a83","black","goldenrod1"))
mycols

samp_val_col <- c("HA"="royalblue3","GF"="goldenrod1")

#### Read in data from CT ####
data1 <- read.csv("/Users/linewulff/Documents/work/projects/2309_CT_GOPlot/inputdata/2508_Fetalmyeloid.csv",header=T,row.names = 1)
head(data1)
data2 <- read.csv("/Users/linewulff/Documents/work/projects/2309_CT_GOPlot/inputdata/2508_25-01mTOR_forHM.csv",header = T, row.names = 1)
head(data2)

# Differentially expressed
Pancreas <- read.csv(paste(dir,"inputdata","Panc_all_genes.csv", sep="/"))
head(Pancreas)
Panc$sign <- "not sign."
Panc[Panc$log2FC_Panc > 0.26303 & Panc$padj < 0.05,]$sign <- 'up' 
Panc[Panc$log2FC_Panc < -0.26303 & Panc$padj < 0.05,]$sign <- 'down' 

# Significant to include
data1_sig <- rownames(data1)[rownames(data1) %in% Panc[Panc$sign!="not sign.",]$g_symbol]
#118
data2_sig <- rownames(data2)[rownames(data2) %in% Panc[Panc$sign!="not sign.",]$g_symbol]
#33

#### Heatmap of significant transcripts ####
anno <- vsd@colData[,c("sc_cluster","replicate")]
colnames(anno) <- c("sc_cluster","replicate")

data_ha = HeatmapAnnotation(colonization=c("GF","GF","GF","GF","HA","HA","HA","HA"),
                             col = list(colonization = c("HA"="royalblue3","GF"="goldenrod1")))
# data1
pdf(paste(dir,"/output/",dato,"_Heatmap_Fetalmyeloid_PancreasDEGs_v1.pdf",sep=""), height = 16, width=5)
Heatmap(t(scale(t(data1[data1_sig,]))),
        col = mycols,
        name="Exp. level",
        show_column_names = FALSE,
        cluster_columns = FALSE,
        show_row_dend = FALSE,
        top_annotation = data_ha)
dev.off()

# data1
pdf(paste(dir,"/output/",dato,"_Heatmap_Fetalmyeloid_PancreasDEGs_v2_wDendro.pdf",sep=""), height = 16, width=7)
Heatmap(t(scale(t(data1[data1_sig,]))),
        col = mycols,
        name="Exp. level",
        show_column_names = FALSE,
        cluster_columns = FALSE,
        #show_row_dend = FALSE,
        row_dend_width = unit(3, "cm"),
        top_annotation = data_ha)
dev.off()


#data2
pdf(paste(dir,"/output/",dato,"_Heatmap_25-01mTOR_PancreasDEGs_v1.pdf",sep=""), height = 6, width=5)
Heatmap(t(scale(t(data2[data2_sig,]))),
        col = mycols,
        name="Exp. level",
        show_column_names = FALSE,
        cluster_columns = FALSE,
        #show_row_dend = FALSE,
        top_annotation = data_ha)
dev.off()

pdf(paste(dir,"/output/",dato,"_Heatmap_25-01mTOR_PancreasDEGs_v2_wDendro.pdf",sep=""), height = 6, width=7)
Heatmap(t(scale(t(data2[data2_sig,]))),
        col = mycols,
        name="Exp. level",
        show_column_names = FALSE,
        cluster_columns = FALSE,
        #show_row_dend = FALSE,
        row_dend_width = unit(3, "cm"),
        top_annotation = data_ha)
dev.off()

