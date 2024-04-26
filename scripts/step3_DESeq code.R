#eagle project - DESeq analysis

########################################## if you haven't installed packages yet ##############
install.packages("htmltools")
library(htmltools)
install.packages("BiocManager")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

#https://lashlock.github.io/compbio/R_presentation.html

########################################## if packages are installed #########################
library(Matrix)
library(BiocManager)
library("DESeq2")
library(ggplot2)

setwd("/Users/emilyruhs/Desktop/UChi_Brook_Lab/GitHub_repos/eagle-contaminant")

#upload data
seqdata <- read.csv('working-data/RUHS_allcounts_edit.csv', header = TRUE, sep = ",")
seqdata2 <- select(countData, c(1,7:30))

head(countData)
str(countData)
#need to convert ensgene to a factor
countData$GeneID= as.factor(countData$GeneID)

sampleinfo = read.csv("working-data/eagle_metadata.csv")
sampleinfo
#-------------------------------------------------------------------------

#Construct DESeqDataSet Object
#this was done in "step2_Normalization preprocessing" Rscript

#Design specifies how the counts from each gene depend on our variables in the metadata
#For this dataset the factor we care about is our treatment status (treatment)

ddsObj

#run DESEQ function
dds <- DESeq(ddsObj)
resultsNames(dds) #"intercept" and "treatment_lps_vs_control" #always puts the second in the alphabet first
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing

#estimateSizeFactors
#This calculates the relative library depth of each sample 
#estimateDispersions
#estimates the dispersion of counts for each gene 
#nbinomWaldTest
#calculates the significance of coefficients in a Negative Binomial GLM using the size and dispersion outputs
#results table
res <- results(dds)
head(results(dds, tidy=TRUE))
# or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="treatment_Null_vs_LPS", type="apeglm")

#summary of differential gene expression
summary(res)

#summary list by p-value
res <- res[order(res$padj),]
head(res)

#plotCounts
#we can use plotCounts fxn to compare the normalized counts
#between treated and control groups for our top 6 genes
dev.new()
par(mfrow=c(2,3))

#rename with genes
dds$treatment <- factor(dds$treatment, levels=c("control", "lps")) #put null before lps

plotCounts(dds, gene="LOC104828915", intgroup="treatment")
plotCounts(dds, gene="LOC104836757", intgroup="treatment")
plotCounts(dds, gene="IRG1", intgroup="treatment")
plotCounts(dds, gene="RSAD2", intgroup="treatment")
plotCounts(dds, gene="LOC104836340", intgroup="treatment")
plotCounts(dds, gene="FMNL3", intgroup="treatment")

#Volcano Plot
#reset par
#dev.new()
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-12,8)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#PCA
#look at how our samples group by treatment
#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation
vsdata <- vst(dds, blind=FALSE)

#dev.new()
par(mfrow=c(1,1))
plotPCA(vsdata, intgroup="treatment") #using the DESEQ2 plotPCA fxn we can




#download and save DEseq data - but only the significant genes
res <- res[order(res$padj),]
resSig <- subset(res, res$padj < 0.1) #setting p-value higher
resSig<-resSig[ order( resSig$log2FoldChange ),]
resSig_p01<-resSig[1:485,] #all the results - 485 genes is the length
write.csv(resSig_p01,"output-data/resSig_p01.csv")


resSig <- subset(res, res$padj < 0.05) #setting p-value higher
resSig<-resSig[ order( resSig$log2FoldChange ),]
resSig_p005<-resSig[1:369,] #all the results - 369 genes is the length
write.csv(resSig_p005,"output-data/resSig_p005.csv")
