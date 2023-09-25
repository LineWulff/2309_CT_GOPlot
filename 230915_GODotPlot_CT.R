## GO enrichment Plots for Carolyn Thompson
## Author: Line Wulff
## Date created: 2023.09.15

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


#### Plot data ####
pdf(paste(dir,"/output/",dato,"GO_DotPlot_sh2.pdf",sep = ""), height = 6,width = 8)
ggplot(data = GO_df, aes(x = Strength, y = Pathway, 
                        color = Ratio, size = -log10(p.adj))) + 
  geom_point() +
  scale_color_viridis_c()+
  theme_bw() + 
  ylab("") + 
  xlab("Strength") 
dev.off()

# p adj conversion
# 3
10^(-3)
# 5
10^(-5)
# 7
10^(-7)