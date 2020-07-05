setwd("C:/Users/durbe/Documents/repositories/StemAway/Deliverables")
library(org.Hs.eg.db) #why are we using a different database?
library(topGO)
library(clusterProfiler)
library(pathview)
library(magrittr)
library(tidyr)

#5.1
data <- read.csv("./data/vp.csv")
genes <- data$X[abs(data$logFC) > .2] #threshold = 2
id <- bitr(genes, fromType="SYMBOL", OrgDb=org.Hs.eg.db, toType=c("ENSEMBL","ENTREZID"))


#5.2
## groupGO()
ggo1 <- groupGO(id$ENTREZID, org.Hs.eg.db, readable=TRUE, level=1)
ggo5 <- groupGO(id$ENTREZID, org.Hs.eg.db, readable=TRUE, level=5)
### -- level increases result in more data (file size went from kB to MB)
barplot(ggo1, main="groupGO with 1 Level")
barplot(ggo5, main="groupGO with 5 Levels")
### -- level increases increase the number of subcategories therefore increasing the number of categories for the bar plot

## enrichGO()
ego <- enrichGO(id$ENTREZID, org.Hs.eg.db, readable=TRUE)
write.csv(as.data.frame(ego), "./data/summaryEnrichGO.csv") #are we supposed to write the results or a summary of the results?
### -- warning message
# In summary(ego) :
#   summary method to convert the object to data.frame is deprecated, please use as.data.frame instead.
egoF <- enrichGO(id$ENTREZID, org.Hs.eg.db, readable=FALSE)
egoFsr <- setReadable(egoF, org.Hs.eg.db, keyType="ENTREZID")
dotplot(egoFsr, title="enrichGO setReadable dotplot")
barplot(egoFsr, main="Top 20 Terms from enrichGO") #only have 8?
plotGOgraph(egoFsr)

#5.3