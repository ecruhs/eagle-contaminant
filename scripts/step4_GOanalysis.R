#DEseq and GO analysis
library(BiocManager)
library(ggplot2)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

install.packages("ggraph") #no from compilation
library(ggraph)
library(clusterProfiler)

#set directory and load data
setwd("/Users/emilyruhs/Desktop/UChi_Brook_Lab/GitHub_repos/eagle-contaminant")

#upload data
DEseq_data <- read.csv('output-data/resSig_p005.csv', header = TRUE, sep = ",")

upreg.dat <- subset(DEseq_data, log2FoldChange>=1.95) #only choosing those with a log2FoldChange >=2
downreg.dat <- subset(DEseq_data, log2FoldChange<0)

upreg.list <- select(upreg.dat, c(2))
downreg.list <- select(downreg.dat, c(2))
head(downreg.list)

#okay there is no reference annotation for bald eagles so we need to make one
library(RSQLite)
library(AnnotationForge)
available.db0pkgs() #chicken is the closest available

makeOrgPackageFromNCBI(version = "0.1",
  author = "Some One <so@someplace.org>",
  maintainer = "Some One <so@someplace.org>",
  outputDir = ".",
  tax_id = "52644",
  genus = "Haliaeetus",
  species = "leucocephalus")

install.packages("./org.Hleucocephalus.eg.db", repos=NULL)

## Makes an organism package for eagle data.frames:
eagleFile <- system.file("extdata","eagle_info.txt",
                         package="AnnotationForge")
eagle <- read.table(eagleFile,sep="\t")
eagle <- read.csv("working-data/eagle_info.csv")

## Now prepare some data.frames
fSym <- finch[,c(2,3,9)]
fSym <- fSym[fSym[,2]!="-",]
fSym <- fSym[fSym[,3]!="-",]
colnames(fSym) <- c("GID","SYMBOL","GENENAME")









#now the GO analysis
#downreg list first since it's smaller
ggo <- groupGO(gene     = downreg.list,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

head(ggo)