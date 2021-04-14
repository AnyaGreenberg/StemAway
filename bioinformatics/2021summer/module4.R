setwd('/Volumes/SamsungUSB/stemaway/bi-mentoring/level1/2021summer')

# library(hgu133plus2.db)
library(biomaRt)
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
rma <- collapseRows(rma[-96], rowGroup=rma$SYMBOL, rowID=rownames(rma))

##### BIOMART
mart <- useEnsembl(biomart='genes', dataset='hsapiens_gene_ensembl')
probe_ids <- rownames(rma)

gn <- getBM(attributes=c('affy_hg_u133_plus_2', 'external_gene_name'),
            filters='affy_hg_u133_plus_2',
            values=probe_ids,
            mart=mart)



### GENE FILTERING
means <- rowMeans(rma$datETcollapsed)
perc2 <- as.numeric(quantile(means, probs=0.02, na.rm=T))
filt <- rma$datETcollapsed[which(means > perc2),]
# write.csv(filt, "DGE/filt.csv")
# filt <- read.csv("DGE/filt.csv", row.names=1)



### LIMMA
meta <- read.csv("Metadata/metadata.txt", row.names=1)
meta <- meta[-c(34,47,77),]
Tissue <- factor(meta$Tissue, levels = c('Normal','Cancer'), ordered = F)
row.names(meta) <- meta$Sample
design <- model.matrix(~Tissue, meta)

lm <- lmFit(filt, design)
# contrast <- makeContrasts(Tissuecancer-Tissuenormal, levels=design)
# fit <- contrasts.fit(lm, contrast)
fit <- eBayes(lm)


tT <- topTable(fit, p.value=0.05, adjust.method="fdr", sort.by="P", genelist=row.names(filt), number=length(row.names(filt)))
# write.table(tT, "DGE/rma.txt")
tT <- read.table("DGE/rma.txt", row.names=1)


## VOLCANO PLOT
EnhancedVolcano(tT, lab=tT$ID, x="logFC", y="adj.P.Val", 
                pointSize=1, legendLabSize=10, labSize=3.0,
                title="Volcano Plot", subtitle="RMA Normalization")



## HEATMAP
tT50 <- topTable(fit, p.value=0.05, adjust.method="fdr", sort.by="P", genelist=row.names(filt), number=50)
input <- filt[(row.names(filt) %in% row.names(tT50)),]
group <- data.frame(Tissue=meta$Group)
row.names(group) <- meta$Sample
pheatmap(input, annotation_col=group, cluster_rows=T, main="Top 50 DEG", 
         annotation_colors=list(Tissue=c(B1_Normal="hotpink1", B1_Cancer="navy", B2_Normal="goldenrod2", B2_Cancer="aquamarine4")))
