setwd("C:/Users/durbe/Documents/repositories/StemAway/Deliverables")
library(GEOquery)
library(affy)

# 1. download series matrix file from GEO 
gse <- getGEO(filename="./meta/GSE8671_series_matrix.txt.gz")

# 2. copy phenotype data and transpose into a new sheet
pheno <- gse@phenoData@data
pheno <- pheno[-c(13,43),]

# 3. remove non-relevant columns/rows
pheno <- pheno[,c(42,44,45)]

# 4. rename rest of columns with meaningful names (avoid spaces)
colnames(pheno) <- c("location","size", "tissue")

# 5. remove unwanted strings from cells (replace empty strings)
for (i in 1:length(pheno$size)) {
  temp <- strsplit(pheno$size[i], split=" ")
  pheno$size[i] <- as.numeric(temp[[1]][1])
}

# 6. save sheet as tab-delimited (.txt)
## write.csv(pheno, "./meta/pheno.csv")
write.table(pheno, "./meta/pheno.txt", sep="\t")

# 7. read file into R using read.table() function
pheno <- read.table("./meta/pheno.txt", sep="\t", header=TRUE)

# 8. have a look into the object - summarize(phenodata)
s <- summary(pheno)

# 9. read in the .CEL files, introduce phenotypic data into the created expression set using the argument pData(expressionSet) <- new("AnnotatedDataFrame", data=object)
filt <- read.csv("./data/filt.csv")
row.names(filt) <-  filt$X
filt <- filt[c(-1,-2)]
filt <- as.matrix(filt)

all(rownames(pheno)==colnames(filt)) 

### METHOD 1
metadata <- data.frame(labelDescription=c("Location of sample", "Size of sample", "Case/Control Status"), row.names=c("location", "size", "score"))
phenoData <- new("AnnotatedDataFrame", data=pheno, varMetadata=metadata)
exprSet <- ExpressionSet(assayData=filt, phenoData=phenoData, annotation="hgu133plus2")

### METHOD 2
minSet <- ExpressionSet(assayData=filt, annotation="hgu133plus2")
pData(minSet) <- pheno

  
# 10. save the expressionSet object (including both expression matrix and phenodata) using saveRDS() function
saveRDS(exprSet, file="./data/phenoData.rds")
saveRDS(minSet, file="./data/phenoData2.rds")

# ### test
# met1 <- readRDS("./data/phenoData.rds")
# met2 <- readRDS("./data/phenoData2.rds")
