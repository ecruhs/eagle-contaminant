
############################################################
#Quality assessment
#From: https://bioinformatics-core-shared-training.github.io/RNAseq_September_2019/html/02_Preprocessing_Data.html#quality-assessment
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(ggfortify)

setwd("/Users/emilyruhs/Desktop/UChi_Brook_Lab/GitHub_repos/eagle-contaminant")

#read in meta data
sampleinfo = read.csv("working-data/eagle_metadata.csv")
sampleinfo

#reading in the count data from featureCounts
seqdata = read.csv("working-data/RUHS_allcounts_edit.csv",comment="#")
seqdata
# In the seqdata object each row represents a gene. The columns contains:
#1 Geneid - Ensembl ID
#2-5 Chr, Start, End, Strand - Genomic locations of exons of the gene
#6 Length - Transcript length of the gene
#7-18 One column for each sample with the count of how many reads were assigned to each gene by featureCounts.

#subset and create new data tables using dplyr
# newTable = sampleinfo
# str(newTable)
# newTable$treatment <- as.factor(newTable$treatment)
# newTable = filter(newTable, treatment=="lps")
# newTable = filter


#now remove the columns that aren't geneid or sample counts
seqdata2 <- select(seqdata, c(1,7:30))

#remove eagle 6 until we know whats going on with the low count values
#seqdata2 <- select(seqdata2, c(1:5, 8:25))

#format the data into a suitable format for DESeq2
#DESeq2 requires a simple object containing only the count data, we’ll keep the gene ID by setting them as the row names
#Let’s create new counts data object, countdata, that contains only the counts for the samples.
#Our animal object contains a column with the sample names. We should adjust the column names of our matrix to match them
#- we just need to remove the .bam suffix.
#It is also critical to ensure that the samples in the columns are in the same order as the rows of sampleinfo. 
#When we load these objects into DESeq2 for the analysis it will not guess which row of the sampleinfo belongs to which 
#column of the counts matrix, it will assume the same order.
#We’ll use to new dplyr and tibble commands:
     #columns_to_rownames to set the rownames using a named column
     #rename_all which allows to rename all the columns using a string function

library(dplyr)
library(tibble)
countdata <- seqdata2 %>%
  column_to_rownames("GeneID") %>% # turn the geneid column into rownames
  as.matrix()

head(countdata)

#For many analysis methods it is advisable to filter out as many genes as possible prior to starting the 
#analysis in order to decrease the impact on false discovery rates when applying multiple testing correction.
#This is normally done by filtering out genes with low numbers of reads, which are likely to be uninformative.
#With DESeq this is not necessary as it applies a process it calls independent filtering during the analysis process. 
#On the other hand, some filtering for genes that are very lowly expressed does reduce the size of the data matrix, 
#meaning that less memory is required and processing steps are carried out faster.
#We will keep all genes where the total number of reads across all samples is greater than 5.

#str(countdata)
dim(countdata)
keep <- rowSums(countdata) > 5
countdata <- countdata[keep,]
dim(countdata)

#assess the quality of the data
#plot how many reads we have per sample
librarySizes <- colSums(countdata)
dev.new()
barplot(librarySizes, 
        names=names(librarySizes), 
        las=2, 
        main="Barplot of library sizes")
abline(h=20e6, lty=2)


#count distribution boxplot
#Count data is not normally distributed, so if we want to examine the distributions of the raw counts it is helpful to 
#transform the data on to a log scale. 

#Typically we use a log2 transformation, however, because the data is count data and will contain many 0s we need to 
#add a count of 1 to every value in order to prevent attempting log2(0) from creating errors.
logcounts <- log2(countdata + 1)

#We’ll check the distribution of read counts using a boxplot and add some colour to see if there is any difference 
#between sample groups.
# make a colour vector
statusCol <- match(sampleinfo$treatment, c("control", "lps")) + 1

#statusCol <- match(sampleinfo$county, c("WINNEBAGO", "OUTAGAMIE", "BROWN", "BAYFIELD", "DOUGLAS")) + 1
# Check distributions of samples using boxplots
boxplot(logcounts, 
        xlab="", 
        ylab="Log2(Counts)",
        las=2,
        col=statusCol)
# Let's add a blue horizontal line that corresponds to the median
abline(h=median(as.matrix(logcounts)), col="black")
#eagle 006_lps is quite low - ask Zarema about this

#A principal component analysis (PCA) is an example of an unsupervised analysis, where we don’t specify the grouping of
#the samples. If your experiment is well controlled and has worked well, we should find that replicate samples cluster 
#closely, whilst the greatest sources of variation in the data should be between treatments/sample groups. 
#It is also an incredibly useful tool for checking for outliers and batch effects.
#To run the PCA we should first normalise our data for library size and transform to a log scale. 
#DESeq2 provides two commands that can be used to do this, here we will use the command rlog. rlog performs a log2 scale 
#transformation in a way that compensates for differences between samples for genes with low read count and also normalizes 
#between samples for library size.
#You can read more about rlog, its alternative vst and the comparison between the two here:
#http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-data-transformations
#To plot the PCA results we will use the autoplot function from the ggfortify package (Tang, Horikoshi, and Li 2016). 
#ggfortify is built on top of ggplot2 and is able to recognise 
#common statistical objects such as PCA results or linear model results and automatically generate summary plot of the 
#results in an appropriate manner.


#make sure deseq2 and ggfortify libraries are loaded
#vstcounts <- vst(countdata)
rlogcounts <- rlog(countdata)
#rlog data "logcounts" has already been normalised for both library size and composition bias


# run PCA
pcDat <- prcomp(t(rlogcounts))
# plot PCA
autoplot(pcDat)

# We can use colour and shape to identify the Type and the Species of each
# sample
autoplot(pcDat,
         data = sampleinfo, 
         colour="treatment", 
         shape="species",
         size=5)

# setting shape to FALSE causes the plot to default to using the labels
library(ggrepel)
#dev.new()
autoplot(pcDat,
         data = sampleinfo, 
         colour="treatment", 
         size=5) + geom_text_repel(aes(x=PC1, y=PC2, label=animal), box.padding = 0.8)

autoplot(pcDat,
         data = sampleinfo, 
         colour="county", 
         size=5) + geom_text_repel(aes(x=PC1, y=PC2, label=animal), box.padding = 0.8)




#create heatmap of the count matrix
# library("pheatmap")
# library("vsn")
# ntd <- normTransform(countdata)
# select <- order(rowMeans(counts(countdata,normalized=TRUE)),
#                 decreasing=TRUE)[1:20]
# df <- as.data.frame(colData(dds)[,c("county","treatment")])
# pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#          cluster_cols=FALSE, annotation_col=df)
