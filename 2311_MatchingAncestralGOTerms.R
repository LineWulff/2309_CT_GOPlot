# GO ancestral terms
## Author: Line Wulff
## Date created: 2023.11.01

#### Libraries ####
library(ggplot2)
library(tidyverse)
library(tidyr)
library(viridis)
library(stringr)
library(dplyr)
library(GOSim)

#### General varaibles ####
dir <- "/Users/linewulff/Documents/work/projects/2309_CT_GOPlot"
setwd(dir)
dato <- str_sub(str_replace_all(Sys.Date(),"-","_"), 3, -1)

#### read in data and data base ####
# Ancestors picked from GOSim database
Anc_ref <- getAncestors()
# Matching of ancestors by GO:id to text terms done with https://www.informatics.jax.org/downloads/reports/go_terms.mgi
GO_terms <- read.csv(paste(dir,"inputdata/BiologicalProcesses.txt",sep="/"), header = F, sep = "\t");colnames(GO_terms) <- c("BP","GOid","term")

head(GO_terms)

# if terms have only GO:0008150 as ancestor - they are the right layer! (one down from Biological processes)
Anc_ref[["GO:0046034"]] #example normal
Anc_ref[["GO:0044848"]] #example top layer

# find all top layer terms
Anc_top <- c()
for (term in names(Anc_ref)){
  if (length(Anc_ref[[term]])==2 & Anc_ref[[term]][1]=="GO:0008150"){
    Anc_top <- c(Anc_top,term)
  }}
length(Anc_top)
# data frame of top layer terms
Anc_top_df <- GO_terms[GO_terms$GOid %in% Anc_top,]
Anc_top_df
dim(Anc_top_df)[1] # 5 missing
# add missing manually
missing <- c("cellular process","multicellular organismal process","biological phase","biological process involved in intraspecies interaction between organisms","biological regulation")
names(missing) <- Anc_top[!Anc_top %in% GO_terms$GOid] 

for (GO in names(missing)){
  Anc_top_df <- Anc_top_df %>% rbind(c("Biological Process",GO,missing[GO]))
}
tail(Anc_top_df)

#### Adding top layers to CT sheets ####
Panc <- read.csv(paste(dir,"inputdata/Panc_100+_1.2_GO.txt",sep="/"), header = T, sep = "\t")
head(Panc)

Panc$Ancestral_Term_GOid <- NA
Panc$Ancestral_Term <- NA
for (term in Panc$X.term.ID){
   Ancestors <- Anc_top_df$GOid[Anc_top_df$GOid %in% Anc_ref[[term]]]
   Panc[Panc$X.term.ID==term,]$Ancestral_Term_GOid <- toString(Ancestors) #allows multiple
   Panc[Panc$X.term.ID==term,]$Ancestral_Term <- toString(Anc_top_df[Anc_top_df$GOid %in% Ancestors,]$term) # allows multiple
}
write.csv(Panc,paste(dir,"output/Panc_100+_1.2_GO_wAncestTerms.csv",sep = "/"))

SI <- read.csv(paste(dir,"inputdata/SI_100+_1.2_GO.txt",sep="/"), header = T, sep = "\t")
head(SI)

SI$Ancestral_Term_GOid <- NA
SI$Ancestral_Term <- NA
for (term in SI$X.term.ID){
  Ancestors <- Anc_top_df$GOid[Anc_top_df$GOid %in% Anc_ref[[term]]]
  SI[SI$X.term.ID==term,]$Ancestral_Term_GOid <- toString(Ancestors) #allows multiple
  SI[SI$X.term.ID==term,]$Ancestral_Term <- toString(Anc_top_df[Anc_top_df$GOid %in% Ancestors,]$term) # allows multiple
}
write.csv(SI,paste(dir,"output/SI_100+_1.2_GO_wAncestTerms.csv",sep = "/"))
