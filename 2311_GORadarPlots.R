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

#### Energy generation ####
## Read in data from CT ##
groupvar <- "EnergyGeneration"
GO_df <- read.csv(paste(dir,"inputdata","Energygeneration_top7_FDR.csv", sep="/"))
colnames(GO_df) <- c("term.ID","Pathway","Pancreas","SI")
head(GO_df)

#each term needs to be a column
GO_rad <- as.data.frame(t(GO_df[,c("Pancreas","SI")]))
colnames(GO_rad) <- GO_df$Pathway

head(GO_rad)
colmax <- ncol(GO_rad) #2:ncol value in following 
GO_rad[,1:colmax] <- -log10(GO_rad[,1:colmax])
min(GO_rad[,1:colmax])
max(GO_rad[,1:colmax])

GO_rad <- cbind(GO_rad, group = rownames(GO_rad))

# ordering radarplot
ord <- c("group",
  # one side
  "Mitochondrion_organization",
  "Generation_of_precursor_metabolites_and_energy",
  "Energy_derivation_by_oxidation_of_organic_compounds",
  "Oxidative_phosphorylation",
  "Cellular_respiration",
  "Proton_motive_force-driven_ATP_synthesis",

  "Mitochondrial_respiratory_chain_complex_assembly"
)

pdf(paste(dir,"/output/",dato,"_GORadar_",groupvar,".pdf",sep = ""), height = 4, width = 4)
ggradar(GO_rad[,ord],
        # Plot area - change according to highest and lowest values
        grid.min = 0, grid.mid = 1.30103, grid.max = 6, #0.1padj=1, 0.05padj.=1.30103
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
dev.off()

#### Immune processes ####
## Read in data from CT ##
groupvar <- "ImmuneProcesses"
GO_df <- read.csv(paste(dir,"inputdata","Immuneprocesses_top7_FDR_1000down.csv", sep="/"))
colnames(GO_df) <- c("term.ID","Pathway","Pancreas","SI")
head(GO_df)

#each term needs to be a column
GO_rad <- as.data.frame(t(GO_df[,c("Pancreas","SI")]))
colnames(GO_rad) <- GO_df$Pathway

head(GO_rad)
colmax <- ncol(GO_rad) #2:ncol value in following 
GO_rad[,1:colmax] <- -log10(GO_rad[,1:colmax])
min(GO_rad[,1:colmax])
max(GO_rad[,1:colmax])

GO_rad <- cbind(GO_rad, group = rownames(GO_rad))

# ordering radarplot
ord <- c("group",
         # one side
         "Regulation_of_MAPK_cascade",

         "Negative_regulation_of_MAPK_cascade",
         "Immune_system_development",
         "Positive_regulation_of_protein_serine/threonine_kinase_activity",
         "Innate_immune_response",
         "Positive_regulation_of_MAPK_cascade",
         "Response_to_cytokine"
)

pdf(paste(dir,"/output/",dato,"_GORadar_",groupvar,".pdf",sep = ""), height = 4, width = 4)
ggradar(GO_rad[,ord],
        # Plot area - change according to highest and lowest values
        grid.min = 0, grid.mid = 1.30103, grid.max = 5, #0.1padj=1, 0.05padj.=1.30103
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
dev.off()
#### Metabolic processes ####
## Read in data from CT ##
groupvar <- "MetabolicProcesses"
GO_df <- read.csv(paste(dir,"inputdata","Metabolicprocesses_top7_FDR_1000down.csv", sep="/"))
colnames(GO_df) <- c("term.ID","Pathway","Pancreas","SI")
head(GO_df)

#each term needs to be a column
GO_rad <- as.data.frame(t(GO_df[,c("Pancreas","SI")]))
colnames(GO_rad) <- GO_df$Pathway

head(GO_rad)
colmax <- ncol(GO_rad) #2:ncol value in following 
GO_rad[,1:colmax] <- -log10(GO_rad[,1:colmax])
min(GO_rad[,1:colmax])
max(GO_rad[,1:colmax])

GO_rad <- cbind(GO_rad, group = rownames(GO_rad))

# ordering radarplot
ord <- c("group",
         # one side
         "Carbohydrate_derivative_metabolic_process",
         "Peptidyl-amino_acid_modification",
         "Cellular_amide_metabolic_process",
         "Carbohydrate_derivative_biosynthetic_process",
         "Positive_regulation_of_proteolysis",
         "Lipid_biosynthetic_process",
         "Positive_regulation_of_phosphate_metabolic_process"
)

pdf(paste(dir,"/output/",dato,"_GORadar_",groupvar,".pdf",sep = ""), height = 4, width = 4)
ggradar(GO_rad[,ord],
        # Plot area - change according to highest and lowest values
        grid.min = 0, grid.mid = 1.30103, grid.max = 14, #0.1padj=1, 0.05padj.=1.30103
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
dev.off()
#### Growth Factor ####
## Read in data from CT ##
groupvar <- "GrowthFactor"
GO_df <- read.csv(paste(dir,"inputdata","Growthfactor_top7_FDR.csv", sep="/"))
colnames(GO_df) <- c("term.ID","Pathway","Pancreas","SI")
head(GO_df)

#each term needs to be a column
GO_rad <- as.data.frame(t(GO_df[,c("Pancreas","SI")]))
colnames(GO_rad) <- GO_df$Pathway

head(GO_rad)
colmax <- ncol(GO_rad) #2:ncol value in following 
GO_rad[,1:colmax] <- -log10(GO_rad[,1:colmax])
min(GO_rad[,1:colmax])
max(GO_rad[,1:colmax])

GO_rad <- cbind(GO_rad, group = rownames(GO_rad))

# ordering radarplot
ord <- c("group",
         # one side
         "Response_to_transforming_growth_factor_beta",
         "Cellular_response_to_growth_factor_stimulus",
         "Transforming_growth_factor_beta_receptor_signaling_pathway",
         "Cellular_response_to_transforming_growth_factor_beta_stimulus",
         "Regulation_of_cellular_response_to_growth_factor_stimulus",
         "Regulation_of_transforming_growth_factor_beta_receptor_signaling_pathway",
         "Response_to_growth_factor"

)

pdf(paste(dir,"/output/",dato,"_GORadar_",groupvar,".pdf",sep = ""), height = 4, width = 4)
ggradar(GO_rad[,ord],
        # Plot area - change according to highest and lowest values
        grid.min = 0, grid.mid = 1.30103, grid.max = 5, #0.1padj=1, 0.05padj.=1.30103
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
dev.off()

#### Cell cycle ####
## Read in data from CT ##
groupvar <- "CellCycle"
GO_df <- read.csv(paste(dir,"inputdata","Cellcycle_top7_FDR_1000down.csv", sep="/"))
colnames(GO_df) <- c("term.ID","Pathway","Pancreas","SI")
head(GO_df)
GO_df <- GO_df[1:7,]

#each term needs to be a column
GO_rad <- as.data.frame(t(GO_df[,c("Pancreas","SI")]))
colnames(GO_rad) <- GO_df$Pathway

head(GO_rad)
colmax <- ncol(GO_rad) #2:ncol value in following 
GO_rad[,1:colmax] <- -log10(GO_rad[,1:colmax])
min(GO_rad[,1:colmax])
max(GO_rad[,1:colmax])

GO_rad <- cbind(GO_rad, group = rownames(GO_rad))

# ordering radarplot
ord <- c("group",
         # one side
         "Regulation of cell cycle phase transition",
         "Regulation of mitotic cell cycle",
         "Cell cycle process",
         "Regulation of cell cycle process",
         "Cell division",
         "Mitotic cell cycle",
         "Mitotic cell cycle process"

)

pdf(paste(dir,"/output/",dato,"_GORadar_",groupvar,".pdf",sep = ""), height = 4, width = 4)
ggradar(GO_rad[,ord],
        # Plot area - change according to highest and lowest values
        grid.min = 0, grid.mid = 1.30103, grid.max = 15, #0.1padj=1, 0.05padj.=1.30103
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
dev.off()

#### Barrier Function ####
## Read in data from CT ##
groupvar <- "BarrierFunction"
GO_df <- read.csv(paste(dir,"inputdata","Barrier_function_top7_FDR_SIbias.csv", sep="/"))
colnames(GO_df) <- c("term.ID","Pathway","Pancreas","SI")
head(GO_df)
#GO_df <- GO_df[1:7,]

#each term needs to be a column
GO_rad <- as.data.frame(t(GO_df[,c("Pancreas","SI")]))
colnames(GO_rad) <- GO_df$Pathway

head(GO_rad)
colmax <- ncol(GO_rad) #2:ncol value in following 
GO_rad[,1:colmax] <- -log10(GO_rad[,1:colmax])
min(GO_rad[,1:colmax])
max(GO_rad[,1:colmax])

GO_rad <- cbind(GO_rad, group = rownames(GO_rad))

# ordering radarplot
ord <- c("group",
         # one side
         "Regulation of Wnt signaling pathway",
         "Epithelial cell development",
         "Regulation of canonical Wnt signaling pathway",
         "Morphogenesis of an epithelium",
         "Regulation of epithelial cell migration",
         "Positive regulation of epithelial cell migration",
         "Epithelial tube morphogenesis"
)

pdf(paste(dir,"/output/",dato,"_GORadar_",groupvar,".pdf",sep = ""), height = 4, width = 4)
ggradar(GO_rad[,ord],
        # Plot area - change according to highest and lowest values
        grid.min = 0, grid.mid = 1.30103, grid.max = 5, #0.1padj=1, 0.05padj.=1.30103
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
dev.off()
