# Setup
## Bioconductor and CRAN libraries used
library(ggplot2)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)

setwd("C:/Users/durbe/Documents/repositories/StemAway/Training/DEanalysis")
## load data

data <- read.table("./data/Mov10_full_counts.txt", header=T, row.names=1)
meta <- read.table("./meta/Mov10_full_meta.txt", header=T, row.names=1)

### check classes of data we just brought in
class(data)
class(meta)

### check that sample names match in both files
all(names(data) %in% rownames(meta))
all(names(data) == rownames(meta))

# DESeq2 - high sensitivity and precision while controlling the false positive rate
## create DESeqDataSet object

dds <- DESeqDataSetFromMatrix(countData=data, colData=meta, design=~sampletype)

View(counts(dds))

# normalization of counts
dds <- estimateSizeFactors(dds) # computed based on median of ratios method - accounts for gene length and sequencing depth
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="./data/normalized_counts.txt", sep="\t", quote=F, col.names=NA)
sizeFactors(dds) #should be identical to those generated in estimateSizeFactors(dds)
estimateSizeFactors(dds) #?
colSums(counts(dds)) #total number of reads for each sample
colSums(counts(dds, normalized=T))

# transformation of counts
rld <- rlog(dds, blind=TRUE)

#quality assessment and exploratory analysis
## PCA
plotPCA(rld, intgroup="sampletype")

### Exercise: Plot the PCA using all of the genes in your original count matrix. Hint: you can use nrow() to help get the total number of genes

## Hierarchical Clustering
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
pheatmap(rld_cor, color=brewer.pal(6,"Blues"))

#running DESeq2
dds <- DESeq(dds)

### specify coefficient or contrast we want to build results table for
res <- results(dds, name="sampletype_MOV10_overexpression_vs_control")
res <- results(dds, contrast=c("sampletype","MOV10_overexpression","control"))

### log fold change shrinkage for visualization and ranking
resLFC <- lfcShrink(dds, coef="sampletype_MOV10_overexpression_vs_control", type="apeglm")

### p-values and adjusted p-values
resOrdered <- res[order(res$pvalue),]
summary(res)

res05 <- results(dds, alpha=0.05)
summary(res05)

# Exploring and Exporting Results
## MA-plot
plotMA(resLFC, ylim=c(-2,2))
?identify # interactively detect row number of individual genes by clicking on the plot


## Volcano Plot
library(EnhancedVolcano)

### build data frame containing logfold change and padjval
df <- as.data.frame(res@listData)
rownames(df) <- rownames(data)
df <- df[complete.cases(df),] #remove nas

### make plot
volcano_plot <- EnhancedVolcano(df,
                               lab=rownames(df),
                               pCutoff=0.01,
                               FCcutoff=0.5,
                               x="log2FoldChange",
                               y="padj",
                               pointSize=1,
                               legendLabSize=10,
                               xlim=c(-6.5,6.5),
                               labSize=3.0)

volcano_plot
