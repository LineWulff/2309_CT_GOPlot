## GO enrichment Plots for Carolyn Thompson
## Author: Line Wulff
## Date created: 2023.09.15

#Example seem online
# create fake data
set.seed(1024) # keep reproducibility
go <- paste0("GO", sample(1000:2000, 5))
data <- data.frame("GOs" = rep(go, 2), 
                   "Condition" = rep(c("A", "B"), each = 5),
                   "GeneRatio" = 1 / sample(10, 10), 
                   "p.adjust" = 0.05 / sample(10, 10))

# plot: dot plot
ggplot(data = data, aes(x = Condition, y = GOs, 
                        color = `p.adjust`, size = GeneRatio)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GO enrichment analysis")

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
GO_df <- read.csv(paste(dir,"inputdata","Metabolismpathways_quicklook.csv", sep="/"))
GO_df <- GO_df[1:20,c("Pathway","Strength","Ratio","p.adj")]
head(GO_df)

GO_df$Pathway <- factor(GO_df$Pathway, levels = c(GO_df[order(GO_df$Strength),]$Pathway))


#### Plot data ####
pdf(paste(dir,"/output/",dato,"GO_DotPlot.pdf",sep = ""), height = 6,width = 8)
ggplot(data = GO_df, aes(x = Strength, y = Pathway, 
                        color = Ratio, size = -log10(p.adj))) + 
  geom_point() +
  scale_color_viridis_c()+
  theme_bw() + 
  ylab("") + 
  xlab("Strength") 
dev.off()
