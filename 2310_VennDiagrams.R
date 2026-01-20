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
library(ggrepel)

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

# same as above but with genes marked like plot 2 below
markDEGS <- read.csv("/Users/linewulff/Documents/work/projects/2309_CT_GOPlot/inputdata/Metabolism_shared.csv", header = F)
head(markDEGS)
markDEGS <- markDEGS$V1
length(markDEGS)

pdf(paste(dir,"/output/",dato,"_Biplot_DEGs_PancreasSI_v2_markedDEGs.pdf",sep=""), height = 4, width = 5)
ggplot(biplot_df[biplot_df$shared!="Not shared",], 
       aes(x = log2FC_SI, y = log2FC_Panc, colour = shared, ))+ 
  geom_point()+
  geom_point(data = biplot_df[biplot_df$g_symbol %in% markDEGS,], shape = 21, colour = "black",fill = "white")+
  scale_colour_manual(values = cols_shared)+
  labs(colour = "Shared DEGs", x = "log2FC SI", y = "log2FC Pancreas")+
  theme_minimal()
dev.off()


# Save csv of biplot wo. "not shared" for Carolyn
write.csv(biplot_df[biplot_df$shared!="Not shared",], 
          file = paste(dir,"/output/",dato,"_biplotdata_DEGsPancreasSI.csv", sep=""))

write.csv(biplot_df[biplot_df$shared=="Not shared",], 
          file = paste(dir,"/output/",dato,"_biplotdata_DEGsPancreasSI_notshared.csv", sep=""))

#### Update with volcano plots for Carolyn's ICMI talk ####
#  Panc_all_genes.csv - (Log2(FC) cut-off > 0.26303 and Padj cut-off <0.05).
## plot 1 - panc genes with cut-offs
Panc$sign <- "not sign."
Panc[Panc$log2FC_Panc > 0.26303 & Panc$padj < 0.05,]$sign <- 'up' 
Panc[Panc$log2FC_Panc < -0.26303 & Panc$padj < 0.05,]$sign <- 'down' 
Panc_cols <- c("not sign."="lightgrey", "up"="royalblue3","down"="goldenrod1")

pdf(paste(dir,"/output/",dato,"_Biplot_PancreasDEGs_v1_onlycols.pdf",sep=""), height = 4, width = 5)
ggplot(Panc, aes(x = log2FC_Panc, y = -log10(padj), colour = sign, ))+ 
  geom_point()+
  # geom_text(data = biplot_df[1:10,], colour = "black")+
  geom_vline(xintercept = c(-0.26303,0.26303), linetype = 'dashed')+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  scale_colour_manual(values = Panc_cols)+
  labs(colour = "DEGs", x = "log2FC Panc", y = "-log10(adj. p-value)")+
  theme_minimal()
dev.off()

## Plot 2 - now with specific genes marked by colour
markDEGS <- read.csv("/Users/linewulff/Documents/work/projects/2309_CT_GOPlot/inputdata/Lipidmetabolism_Panc4volc.csv", header = F)
head(markDEGS)
markDEGS <- markDEGS$V1
length(markDEGS)

pdf(paste(dir,"/output/",dato,"_Biplot_PancreasDEGs_v2_markedDEGs.pdf",sep=""), height = 4, width = 5)
ggplot(Panc, aes(x = log2FC_Panc, y = -log10(padj), colour = sign, ))+ 
  geom_point()+
  geom_point(data = Panc[Panc$g_symbol %in% markDEGS,], shape = 21, colour = "black",fill = "white")+
  #geom_text(data = Panc[Panc$g_symbol %in% markDEGS,], colour = "black")+
  geom_vline(xintercept = c(-0.26303,0.26303), linetype = 'dashed')+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  scale_colour_manual(values = Panc_cols)+
  labs(colour = "DEGs", x = "log2FC Panc", y = "-log10(adj. p-value)")+
  theme_minimal()
dev.off()

## Plot 3 - biplot with same colouring
cols_shared_new <- c("SI up, PC up"="royalblue3","SI down, PC down"="goldenrod1",
                 "SI down, PC up"="lightgrey","SI up, PC down"="lightgrey",
                 "Not shared"="lightgrey")

pdf(paste(dir,"/output/",dato,"_Biplot_DEGs_PancreasSI_v3.pdf",sep=""), height = 4, width = 5)
ggplot(biplot_df[biplot_df$shared!="Not shared",], 
       aes(x = log2FC_SI, y = log2FC_Panc, colour = shared, ))+ #label = g_symbol))+
  geom_point()+
  # geom_text(data = biplot_df[1:10,], colour = "black")+
  scale_colour_manual(values = cols_shared_new)+
  labs(colour = "Shared DEGs", x = "log2FC SI", y = "log2FC Pancreas")+
  theme_minimal()
dev.off()

markDEGS <- read.csv("/Users/linewulff/Documents/work/projects/2309_CT_GOPlot/inputdata/Metabolism_shared.csv", header = F)
head(markDEGS)
markDEGS <- markDEGS$V1
length(markDEGS)

pdf(paste(dir,"/output/",dato,"_Biplot_DEGs_PancreasSI_v3_markedDEGs.pdf",sep=""), height = 4, width = 5)
ggplot(biplot_df[biplot_df$shared!="Not shared",], 
       aes(x = log2FC_SI, y = log2FC_Panc, colour = shared, ))+ 
  geom_point()+
  geom_point(data = biplot_df[biplot_df$g_symbol %in% markDEGS,], shape = 21, colour = "black",fill = "white")+
  scale_colour_manual(values = cols_shared_new)+
  labs(colour = "Shared DEGs", x = "log2FC SI", y = "log2FC Pancreas")+
  theme_minimal()
dev.off()


## Plot 4 - same as plot 3 but with DEGs marked
markDEGS_sh <- read.csv("/Users/linewulff/Documents/work/projects/2309_CT_GOPlot/inputdata/Lipidmetabolism_shared4biplot.csv", header = F)
head(markDEGS_sh)
markDEGS_sh <- markDEGS_sh$V1
length(markDEGS_sh)

pdf(paste(dir,"/output/",dato,"_Biplot_DEGs_PancreasSI_markedDEGS_v4.pdf",sep=""), height = 4, width = 5)
ggplot(biplot_df[biplot_df$shared!="Not shared",], 
       aes(x = log2FC_SI, y = log2FC_Panc, colour = shared, ))+ #label = g_symbol))+
  geom_point()+
  geom_point(data = biplot_df[biplot_df$g_symbol %in% markDEGS_sh,], shape = 21, colour = "black", fill = "white")+
  # geom_text(data = biplot_df[1:10,], colour = "black")+
  scale_colour_manual(values = cols_shared_new)+
  labs(colour = "Shared DEGs", x = "log2FC SI", y = "log2FC Pancreas")+
  theme_minimal()
dev.off()

## Plot 5 - same as above but points labelled
pdf(paste(dir,"/output/",dato,"_Biplot_DEGs_PancreasSI_markedandlabsDEGs_v5.pdf",sep=""), height = 4, width = 5)
ggplot(biplot_df[biplot_df$shared!="Not shared",], 
       aes(x = log2FC_SI, y = log2FC_Panc, colour = shared, label = g_symbol))+
  geom_point()+
  geom_point(data = biplot_df[biplot_df$g_symbol %in% markDEGS_sh,], shape = 21, colour = "black", fill = "white")+
  geom_label_repel(data = biplot_df[biplot_df$g_symbol %in% markDEGS_sh,], colour = "black", max.overlaps = 35)+
  scale_colour_manual(values = cols_shared_new)+
  labs(colour = "Shared DEGs", x = "log2FC SI", y = "log2FC Pancreas")+
  theme_minimal()
dev.off()

#### Same as for ICMI talk but SI data - only plot 1 ####
#  SI_all_genes.csv - (Log2(FC) cut-off > 0.26303 and Padj cut-off <0.05).
## plot 1 - panc genes with cut-offs
SI$sign <- "not sign."
SI[SI$log2FC_SI > 0.26303 & SI$padj < 0.05,]$sign <- 'up' 
SI[SI$log2FC_SI < -0.26303 & SI$padj < 0.05,]$sign <- 'down' 
SI_cols <- c("not sign."="lightgrey", "up"="royalblue3","down"="goldenrod1")

pdf(paste(dir,"/output/",dato,"_Biplot_SIDEGs_v1_onlycols.pdf",sep=""), height = 4, width = 5)
ggplot(SI, aes(x = log2FC_SI, y = -log10(padj), colour = sign, ))+ 
  geom_point()+
  # geom_text(data = biplot_df[1:10,], colour = "black")+
  geom_vline(xintercept = c(-0.26303,0.26303), linetype = 'dashed')+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  scale_colour_manual(values = SI_cols)+
  labs(colour = "DEGs", x = "log2FC SI", y = "-log10(adj. p-value)")+
  theme_minimal()
dev.off()

#### Update 13/8-25 - Pancreas plot edits
## Remove Col6a5
## Add labels if passing significance thresholds
Panc$sign <- "not sign."
Panc[Panc$log2FC_Panc > 0.26303 & Panc$padj < 0.05,]$sign <- 'up' 
Panc[Panc$log2FC_Panc < -0.26303 & Panc$padj < 0.05,]$sign <- 'down' 
Panc_cols <- c("not sign."="lightgrey", "up"="royalblue3","down"="goldenrod1")

pdf(paste(dir,"/output/",dato,"_Biplot_PancreasDEGs_v2_Col6a5rem_onlycols.pdf",sep=""), height = 4, width = 5)
ggplot(Panc[Panc$g_symbol!="Col6a5",], aes(x = log2FC_Panc, y = -log10(padj), colour = sign, ))+ 
  geom_point()+
  # geom_text(data = biplot_df[1:10,], colour = "black")+
  geom_vline(xintercept = c(-0.26303,0.26303), linetype = 'dashed')+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  scale_colour_manual(values = Panc_cols)+
  labs(colour = "DEGs", x = "log2FC Panc", y = "-log10(adj. p-value)")+
  theme_minimal()
dev.off()

markDEGS_v2 <- read.csv("/Users/linewulff/Documents/work/projects/2309_CT_GOPlot/inputdata/Betacellgenes_forvolcanoplot.csv",header = F)
head(markDEGS_v2)
markDEGS_v2 <- markDEGS_v2$V1
length(markDEGS_v2) #39
dim(Panc[Panc$g_symbol %in% markDEGS_v2 & Panc$sign!="not sign.",]) #22 significant

pdf(paste(dir,"/output/",dato,"_Biplot_PancreasDEGs_v2_Col6a5rem_selectDEGlabels.pdf",sep=""), height = 4, width = 5)
ggplot(Panc[Panc$g_symbol!="Col6a5",], aes(x = log2FC_Panc, y = -log10(padj), colour = sign, label = g_symbol))+ 
  geom_point()+
  # geom_text(data = biplot_df[1:10,], colour = "black")+
  geom_vline(xintercept = c(-0.26303,0.26303), linetype = 'dashed')+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_label_repel(data = Panc[Panc$g_symbol %in% markDEGS_v2 & Panc$sign!="not sign.",], colour = "black", max.overlaps = 50)+
  scale_colour_manual(values = Panc_cols)+
  labs(colour = "DEGs", x = "log2FC Panc", y = "-log10(adj. p-value)")+
  theme_minimal()
dev.off()

#### 26-01-20 - Volcano plots for beta cell genes ####
markDEGS_v3 <- read.csv("/Users/linewulff/Documents/work/projects/2309_CT_GOPlot/inputdata/Slcforvolcanoplot.csv",header = F)
markDEGS_v2 <- read.csv("/Users/linewulff/Documents/work/projects/2309_CT_GOPlot/inputdata/Betacellgenes_forvolcanoplot.csv",header = F)
markDEGS_v2 <- markDEGS_v2$V1; markDEGS_v3 <- markDEGS_v3$V1

pdf(paste(dir,"/output/",dato,"_volcplot_PancreasDEGs_v2_Col6a5rem_slcandbetacellsgenes.pdf",sep=""), height = 4, width = 5)
ggplot(Panc[Panc$g_symbol!="Col6a5",], aes(x = log2FC_Panc, y = -log10(padj), colour = sign, label = g_symbol))+ 
  geom_point()+
  geom_point(data = Panc[Panc$g_symbol %in% markDEGS_v2 & Panc$sign!="not sign.",], fill = "darkolivegreen",colour="black", shape = 21)+
  geom_point(data = Panc[Panc$g_symbol %in% markDEGS_v3 & Panc$sign!="not sign.",], fill = "firebrick",colour="black", shape =21)+
  # geom_text(data = biplot_df[1:10,], colour = "black")+
  geom_vline(xintercept = c(-0.26303,0.26303), linetype = 'dashed')+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_label_repel(data = Panc[Panc$g_symbol %in% markDEGS_v2 & Panc$sign!="not sign.",], colour = "darkolivegreen", max.overlaps = 100)+
  geom_label_repel(data = Panc[Panc$g_symbol %in% markDEGS_v3 & Panc$sign!="not sign.",], colour = "firebrick", max.overlaps = 100)+
  scale_colour_manual(values = Panc_cols)+
  labs(colour = "DEGs", x = "log2FC Panc", y = "-log10(adj. p-value)")+
  theme_minimal()
dev.off()

pdf(paste(dir,"/output/",dato,"_volcplot_PancreasDEGs_v2_Col6a5rem_slcgenes.pdf",sep=""), height = 4, width = 5)
ggplot(Panc[Panc$g_symbol!="Col6a5",], aes(x = log2FC_Panc, y = -log10(padj), colour = sign, label = g_symbol))+ 
  geom_point()+
  geom_point(data = Panc[Panc$g_symbol %in% markDEGS_v3 & Panc$sign!="not sign.",],colour="black", shape =21)+
  # geom_text(data = biplot_df[1:10,], colour = "black")+
  geom_vline(xintercept = c(-0.26303,0.26303), linetype = 'dashed')+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_label_repel(data = Panc[Panc$g_symbol %in% markDEGS_v3 & Panc$sign!="not sign.",], colour = "black", max.overlaps = 100)+
  scale_colour_manual(values = Panc_cols)+
  labs(colour = "DEGs", x = "log2FC Panc", y = "-log10(adj. p-value)")+
  theme_minimal()
dev.off()

pdf(paste(dir,"/output/",dato,"_volcplot_PancreasDEGs_v2_Col6a5rem_betacellgenes.pdf",sep=""), height = 4, width = 5)
ggplot(Panc[Panc$g_symbol!="Col6a5",], aes(x = log2FC_Panc, y = -log10(padj), colour = sign, label = g_symbol))+ 
  geom_point()+
  geom_point(data = Panc[Panc$g_symbol %in% markDEGS_v2 & Panc$sign!="not sign.",],colour="black", shape =21)+
  # geom_text(data = biplot_df[1:10,], colour = "black")+
  geom_vline(xintercept = c(-0.26303,0.26303), linetype = 'dashed')+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_label_repel(data = Panc[Panc$g_symbol %in% markDEGS_v2 & Panc$sign!="not sign.",], colour = "black", max.overlaps = 100)+
  scale_colour_manual(values = Panc_cols)+
  labs(colour = "DEGs", x = "log2FC Panc", y = "-log10(adj. p-value)")+
  theme_minimal()
dev.off()





