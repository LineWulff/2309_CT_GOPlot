## GO enrichment Plots for CT
## Creation of radarplot based selected terms by CT.
## Author: Line Wulff
## Date created: 2023.10.26

#### Libraries ####
library(ggplot2)
library(tidyverse)
library(tidyr)
library(viridis)
library(stringr)
library(ggpubr)
library(ggradar)

#### General varaibles ####
dir <- "/Users/linewulff/Documents/work/projects/2309_CT_GOPlot"
setwd(dir)
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)

#### Read in data from CT ####
GO_df_SI <- read.csv(paste(dir,"inputdata","Metabolismpathways_quicklook_sh2.csv", sep="/"))
GO_df_SI <- GO_df_SI[1:20,c("Pathway","Strength","Ratio","p.adj")]
GO_df_SI$pop <- "SI"
head(GO_df_SI)

#GO_df_SI$Pathway <- factor(GO_df_SI$Pathway, levels = c(GO_df_SI[order(GO_df_SI$Strength),]$Pathway))



terms <- c(GO_df_SI$Pathway, GO_df_PC$Pathway)

## DF for rdarplot
GO_rad <- GO_df_SI[,c("pop","p.adj","Pathway")] #GO_df_SI$Pathway %in% Pathways
missing_t <- terms[!terms %in% GO_df_SI$Pathway]
for (item in missing_t){
  GO_rad <- rbind(GO_rad, c("SI",1,item))}

GO_rad <- rbind(GO_rad, GO_df_PC[,c("pop","p.adj","Pathway")])
missing_t <- terms[!terms %in% GO_df_PC$Pathway]
for (item in missing_t){
  GO_rad <- rbind(GO_rad, c("PC",1,item))}

## Adjusting dataframe
head(GO_rad)
GO_rad$p.adj <- as.numeric(GO_rad$p.adj)
rownames(GO_rad) <- NULL

#each term needs to be a column
GO_rad <- spread(GO_rad, key = Pathway, value = p.adj)

head(GO_rad)
colmax <- ncol(GO_rad) #2:ncol value in following 
GO_rad[,2:colmax] <- -log10(GO_rad[,2:colmax])
min(GO_rad[,2:colmax])
max(GO_rad[,2:colmax])

# ordering radarplot
ord <- c("pop",
  # one side
  "Antigen processing-Cross presentation" ,                                      
  
  "Activation of NF-kappaB in B cells",
  "TRIF(TICAM1)-mediated TLR4 signaling", 
  "MyD88 dependent cascade initiated on endosome", 
  "CLEC7A (Dectin-1) signaling", 
  "MyD88 cascade initiated on plasma membrane", 
  "Dectin-1 mediated noncanonical NF-kB signaling",                              
  "ER-Phagosome pathway",

  "FCERI mediated NF-kB activation"  ,                                           
  "MyD88:MAL(TIRAP) cascade initiated on plasma membrane" , 
  "MAPK targets/ Nuclear events mediated by MAP kinases",
  "FGFR2 alternative splicing",
  "DDX58/IFIH1-mediated induction of interferon-alpha/beta",
  
  "TNFR1-induced NFkappaB signaling pathway" ,
  "NEW" ,
  "Antigen Presentation: Folding, assembly and peptide loading of class I MHC",
  "MAP kinase activation" , 
  "TRAF6 mediated induction of NFkB and MAP kinases upon TLR7/8 or 9 activation",
  "ROS and RNS production in phagocytes",
  "Noncanonical NF-kB signaling" ,                                               
  "Toll Like Receptor 3 (TLR3) Cascade"
)

ggradar(GO_rad[,ord],
        # Plot area - change according to highest and lowest values
        grid.min = 0, grid.mid = 1.30103, grid.max = 4, #0.1padj=1, 0.05padj.=1.30103
        # Polygons
        group.line.width = 1, 
        group.point.size = 0,
        group.colours = c("blue3","coral"), 
        # Background and grid lines
        background.circle.colour = "white",
        gridline.mid.colour = "black", #if set to logFC threshold set  to different colour
        # Text and labels
        #base.size = 5,
        legend.position = "bottom",
        label.centre.y = FALSE,
        label.gridline.min = FALSE,
        label.gridline.mid = FALSE,
        label.gridline.max = FALSE,
        #values.radar = c("1", "p.adj=0.05", "lowest"),
        axis.label.size = 4,
        axis.label.offset =1.05)
