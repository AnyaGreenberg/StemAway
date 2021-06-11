setwd('/Volumes/SamsungUSB/stemaway/bi-mentoring/level1/2021summer')

library(hgu133plus2.db)
library(WGCNA)
library(limma)
library(EnhancedVolcano)
library(ggplot2)
library(pheatmap)

### REMOVE OUTLIERS
rma <- read.csv("data/rma.csv", row.names=1)
paper <- read.csv("data/selected-rma.csv", row.names=1)



### ANNOTATION
symbols <- AnnotationDbi::select(hgu133plus2.db, keys=rownames(rma), columns=c("SYMBOL"))
symbols <- symbols[!duplicated(symbols$PROBEID),]

if (all(row.names(rma) == symbols$PROBEID)) {
  rma$SYMBOL <- symbols$SYMBOL
}

rma <- na.omit(rma)
rma <- collapseRows(rma[-121], rowGroup=rma$SYMBOL, rowID=rownames(rma))

###

symbols.selected <- AnnotationDbi::select(hgu133plus2.db, keys=rownames(paper), columns=c("SYMBOL"))
symbols.selected <- symbols.selected[!duplicated(symbols.selected$PROBEID),]

if (all(row.names(paper) == symbols.selected$PROBEID)) {
  paper$SYMBOL <- symbols.selected$SYMBOL
}

paper <- na.omit(paper)
paper <- collapseRows(paper[-63], rowGroup=paper$SYMBOL, rowID=rownames(paper))



### GENE FILTERING
means <- rowMeans(rma$datETcollapsed)
perc2 <- as.numeric(quantile(means, probs=0.02, na.rm=T))
filt <- rma$datETcollapsed[which(means > perc2),]
# write.csv(filt, "dge/filt.csv")
filt <- read.csv("dge/filt.csv", row.names=1)

###

means <- rowMeans(paper$datETcollapsed)
perc2 <- as.numeric(quantile(means, probs=0.02, na.rm=T))
filt.selected <- paper$datETcollapsed[which(means > perc2),]
# write.csv(filt.selected, "dge/filt-selected.csv")
filt.selected <- read.csv("dge/filt-selected.csv", row.names=1)


### LIMMA
meta <- read.csv("metadata/metadata.txt", row.names=1)
Tissue <- factor(meta$Tissue, levels = c('NORMAL','CANCER'), ordered = F)
row.names(meta) <- meta$Sample
design <- model.matrix(~Tissue, meta)

lm <- lmFit(filt, design)
# contrast <- makeContrasts(Tissuecancer-Tissuenormal, levels=design)
# fit <- contrasts.fit(lm, contrast)
fit <- eBayes(lm)


tT <- topTable(fit, p.value=0.05, adjust.method="fdr", sort.by="P", genelist=row.names(filt), number=length(row.names(filt)))
# write.table(tT, "dge/rma-tt.txt")
tT <- read.table("dge/rma-tt.txt", row.names=1)

###

meta.selected <- read.csv("metadata/selected-metadata.csv", row.names=1)
Tissue.selected <- factor(meta.selected$Tissue, levels = c('NORMAL','CANCER'), ordered = F)
row.names(meta.selected) <- meta.selected$Sample
design.selected <- model.matrix(~Tissue, meta.selected)

lm.selected <- lmFit(filt.selected, design.selected)
fit.selected <- eBayes(lm.selected)


tT.selected <- topTable(fit.selected, p.value=0.05, adjust.method="fdr", sort.by="P", genelist=row.names(filt.selected), number=length(row.names(filt.selected)))
# write.table(tT.selected, "dge/rma-tt-selected.txt")
tT.selected <- read.table("dge/rma-tt-selected.txt", row.names=1)



## VOLCANO PLOT
EnhancedVolcano(tT, lab=tT$ID, x="logFC", y="adj.P.Val", 
                pointSize=1, legendLabSize=10, labSize=3.0,
                title="Volcano Plot", subtitle="RMA Normalization")

###

EnhancedVolcano(tT.selected, lab=tT.selected$ID, x="logFC", y="adj.P.Val", 
                pointSize=1, legendLabSize=10, labSize=3.0,
                title="Volcano Plot", subtitle="RMA Normalization (selected samples)")



## HEATMAP
tT50 <- topTable(fit, p.value=0.05, adjust.method="fdr", sort.by="P", genelist=row.names(filt), number=50)
input <- filt[(row.names(filt) %in% row.names(tT50)),]
group <- data.frame(Tissue=meta$Tissue)
row.names(group) <- meta$Sample
pheatmap(input, annotation_col=group, cluster_rows=T, main="Top 50 DEG", 
         annotation_colors=list(Tissue=c(CANCER="hotpink1", NORMAL="navy")))

###
tT50.selected <- topTable(fit.selected, p.value=0.05, adjust.method="fdr", sort.by="P", genelist=row.names(filt.selected), number=50)
input.selected <- filt.selected[(row.names(filt.selected) %in% row.names(tT50.selected)),]
group.selected <- data.frame(Tissue=meta.selected$Tissue)
row.names(group.selected) <- meta.selected$Sample
pheatmap(input.selected, annotation_col=group.selected, cluster_rows=T, main="Top 50 DEG (selected samples)", 
         annotation_colors=list(Tissue=c(CANCER="hotpink1", NORMAL="navy")))
