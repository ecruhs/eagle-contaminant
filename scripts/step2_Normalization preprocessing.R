
############################################################
#Quality assessment
#From: https://bioinformatics-core-shared-training.github.io/RNAseq_September_2019/html/02_Preprocessing_Data.html#quality-assessment
library(DESeq2)
library(tidyverse)
library(tibble)

setwd("/Users/emilyruhs/Desktop/UChi_Brook_Lab/GitHub_repos/eagle-contaminant")

#read in meta data
sampleinfo = read.csv("working-data/eagle_metadata.csv")
sampleinfo

#reading in the count data from featureCounts
seqdata = read.csv("working-data/RUHS_allcounts_edit.csv",comment="#")
seqdata
seqdata2 <- select(seqdata, c(1,7:30))
# In the seqdata object each row represents a gene. The columns contains:
#1 Geneid - Ensembl ID
#2-5 Chr, Start, End, Strand - Genomic locations of exons of the gene
#6 Length - Transcript length of the gene
#7-18 One column for each sample with the count of how many reads were assigned to each gene by featureCounts.

countdata <- seqdata %>%
  column_to_rownames("GeneID") %>% # turn the geneid column into rownames
  as.matrix()


##########################################################################
#NORMALIZATION
#convert counts to DESeqDataSet object
# first lets check that our rows and columns match
all(sampleinfo$animal == colnames(countdata))
#should come back as true

#a simple model where would just be concerned with contrasting the differential expression between treatments (null and lps). 
# create the design formula
#design <- as.formula(~ treatment)
# create the DESeqDataSet object
ddsObj <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = sampleinfo,
                                 design = ~treatment)

#The estimateSizeFactors command of DESeq2 applies the “median ratio method” of normalisation(Huber 2010). 
#This generates a set of normalization factors and adds them to the ddsObj in the colData slot.

# Apply normalisation to DDS object
ddsObj <- estimateSizeFactors(ddsObj)
sizeFactors(ddsObj)
#A normalization factor below one indicates that the library size will be scaled down, as there is more suppression (i.e., composition bias)
#in that library relative to the other libraries. This is also equivalent to scaling the counts upwards in that sample. 
#Conversely, a factor above one scales up the library size and is equivalent to downscaling the counts.
normalizedCounts=counts(ddsObj,normalized=TRUE)
logNormalizedCounts=log2(normalizedCounts+1)



save(countdata, sampleinfo, ddsObj,normalizedCounts,logNormalizedCounts, file="preprocessing.RData")
write.csv(normalizedCounts, file="output-data/normalizedcounts.csv")
write.csv(logNormalizedCounts, file="output-data/lognormalizedcounts.csv")
