setwd("C:/Users/durbe/Documents/repositories/StemAway/Deliverables")
library(hgu133plus2.db)
library(limma)
library(EnhancedVolcano)

norm <- read.csv("./data/norm.csv")
colnames(norm)[1] <- "PROBEID"
norm <- norm[,-44]
norm <- norm[,-14]

# annotation
keys <- norm$PROBEID
ps <- select(hgu133plus2.db, keys=keys, columns=c("SYMBOL"), keytype="PROBEID")

## create new column that contains 1 symbol for each probeset id
### - the first probset id is assigned to each symbol
symbol <- merge(ps, norm, by="PROBEID", all.x=TRUE)
symbol <- symbol[!duplicated(symbol$PROBEID),]
symbol <- symbol[!duplicated(symbol$SYMBOL),]

# gene filtering
means <- rowMeans(symbol[,3:64])
perc4 <- as.numeric(quantile(means, probs=0.04, na.rm=TRUE))
filt <- symbol[which(means > perc4),]
row.names(filt) <- filt$PROBEID
filt <- filt[-1]
write.csv(filt, "./data/filt.csv")


# limma -- https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE8671
genelist <- filt$SYMBOL
filt <- filt[-1]
type <- factor(paste(rep(c("control", "cancer"), each=31)))

## compute pooled variance
design <- model.matrix(~0+type)
write.csv(design, "./data/model_matrix.csv")

corfit <- duplicateCorrelation(filt, design)

lm <- lmFit(filt, design, cor=corfit$consensus.correlation)

## create coefficient matrix for contrasts
contrast <- makeContrasts(typecancer-typecontrol, levels=design)
colnames(contrast) <- c("Ca-Co")

## compute estimated contrasts
fitC <- contrasts.fit(lm, contrast)

## compute moderated contrast t-test
efitC <- eBayes(fitC)

## plot p-values
hist(efitC$p.value[,1], main="Cancer vs. Control")

## compute gene list
tT05 <- topTable(efitC, adjust.method="fdr", p.value=0.05, genelist=genelist, n=100)
write.csv(tT05, "./data/topT_05.csv")

tT1e5 <- topTable(efitC, adjust.method="fdr", p.value=1e-5, genelist=genelist, n=100)
write.csv(tT1e5, "./data/topT_1e-5.csv")

## volcano plot -- made with all data
vp <- topTable(efitC, adjust.method="fdr", genelist=genelist, n=52488)

volcano_plot <- EnhancedVolcano(vp,
                                lab=vp$ID,
                                x="logFC",
                                y="adj.P.Val",
                                pCutoff=0.01,
                                FCcutoff=1.0,
                                pointSize=1,
                                legendLabSize=10,
                                labSize=3.0,
                                subtitle="alpha = 0.01")

volcano_plot
