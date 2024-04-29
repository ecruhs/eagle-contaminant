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
seqdata <- read.csv('working-data/RUHS_lpscounts_edit.csv', header = TRUE, sep = ",")
seqdata2 <- select(countData, c(1,7:30))

head(countdata)
str(countdata)
#need to convert ensgene to a factor
countdata$GeneID= as.factor(countdata$GeneID)

sampleinfo = read.csv("working-data/eagle_metadata_lps.csv")
sampleinfo
#-------------------------------------------------------------------------

#Construct DESeqDataSet Object
#this was done in "step2_Normalization preprocessing" Rscript

#Design specifies how the counts from each gene depend on our variables in the metadata
#For this dataset the factor we care about is our treatment status (treatment)

ddsObj

#run DESEQ function
dds <- DESeq(ddsObj)
resultsNames(dds) #"intercept" and "LC_cat_b_vs_a" #always puts the second in the alphabet first
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
dds$LC_cat <- factor(dds$LC_cat, levels=c("a", "b")) #put null before lps

plotCounts(dds, gene="HPGDS", intgroup="LC_cat")
plotCounts(dds, gene="DCN", intgroup="LC_cat")
plotCounts(dds, gene="PLOD2", intgroup="LC_cat")
plotCounts(dds, gene="FLT1", intgroup="LC_cat")
plotCounts(dds, gene="IGSF9B", intgroup="LC_cat")
plotCounts(dds, gene="CARHSP1", intgroup="LC_cat")

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
plotPCA(vsdata, intgroup="LC_cat") #using the DESEQ2 plotPCA fxn we can




#download and save DEseq data - but only the significant genes
res <- res[order(res$padj),]
resSig <- subset(res, res$padj < 0.1) #setting p-value higher
resSig<-resSig[ order( resSig$log2FoldChange ),]
resSig_p01<-resSig[1:1822,] #all the results - 1822 genes is the length
write.csv(resSig_p01,"output-data/resSig_p01.csv")


resSig <- subset(res, res$padj < 0.05) #setting p-value higher
resSig<-resSig[ order( resSig$log2FoldChange ),]
resSig_p005<-resSig[1:975,] #all the results - 975 genes is the length
write.csv(resSig_p005,"output-data/resSig_p005.csv")
