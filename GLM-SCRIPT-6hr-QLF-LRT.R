library(statmod)
library(edgeR)
library(gplots)
library(RColorBrewer)
library(foreign)
library(bigPint)
library(tidyverse)
library(dendextend)
library(reshape2)
library(mixOmics)
library(HTSFilter)

readCount06 <- read.table("gm_6h_B1-B34.csv", header=TRUE, sep=",", row.names=1)
Treat06 <- factor(substring(colnames(readCount06),1,2))
Treat06<-relevel(Treat06,ref = "B.")

y06<-DGEList(counts = readCount06, group=Treat06)
keep06 <- filterByExpr(y06)
y06 <- y06[keep06, ,keep.lib.sizes=FALSE]
y06 <- calcNormFactors(y06)


design06 <- model.matrix(~0+Treat06)


rownames(design06)<-colnames(y06)
y06 <- estimateDisp(y06, design06, robust=TRUE)

fit06L<- glmFit(y06, design06)
CONTRASTS06L <- makeContrasts( BvsP6 = Treat06P.-Treat06B.,
                               BvsH6 = Treat06H.-Treat06B.,
                               HvsP6 = Treat06P.-Treat06H.,
                               levels = design06)

GLM06L<- paste(colnames(CONTRASTS06L),sep="")
for (i in 1:ncol(CONTRASTS06L)){
  current.glmQLFTest <- glmLRT(fit06L, contrast = CONTRASTS06L[,i])
  outputFile06 <- paste0((colnames(CONTRASTS06L)[i]), "06hr.csv")
  write.table(current.glmQLFTest, file=outputFile06, sep=',', quote=FALSE)
  assign(GLM06L[i], current.glmQLFTest)}

topBvsP6<- topTags(BvsP6, n=nrow(BvsP6$table), adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
topBvsH6<- topTags(BvsH6, n=nrow(BvsH6$table), adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
topHvsP6<- topTags(HvsP6, n=nrow(HvsP6$table), adjust.method = "BH", sort.by = "PValue", p.value = 0.05)


topBvsP6np<- topTags(BvsP6, n=nrow(BvsP6$table), adjust.method = "BH", sort.by = "PValue")
topBvsH6np<- topTags(BvsH6, n=nrow(BvsH6$table), adjust.method = "BH", sort.by = "PValue")
topHvsP6np<- topTags(HvsP6, n=nrow(HvsP6$table), adjust.method = "BH", sort.by = "PValue")


library(readr)
CombiMarkers <- read_table2("CombiMarkers.csv",
                            col_names = FALSE)
View(CombiMarkers)
CombiMarkers$X1<-paste(CombiMarkers$X1,"v2" ,sep=".")


CombiMarkersBP6<-subset(topBvsP6np$table, rownames(topBvsP6np$table) %in% CombiMarkers$X1 )
CombiMarkersBP6L<-CombiMarkersBP6[,1]
CombiMarkersBH6<-subset(topBvsH6np$table, rownames(topBvsH6np$table) %in% CombiMarkers$X1 )
CombiMarkersBH6L<-CombiMarkersBH6[,1]
CombiMarkersHP6<-subset(topHvsP6np$table, rownames(topHvsP6np$table) %in% CombiMarkers$X1 )
CombiMarkersHP6L<-CombiMarkersHP6[,1]



CombiMarkersBP1L<-data.matrix(CombiMarkersBP1L)
rownames(CombiMarkersBP1L)<-rownames(CombiMarkersBP1)

#combine LogFCs in excel and read in as heatmap 6 
heatmap6<-as.data.frame(heatmap6)
rownames(heatmap6)<-heatmap6[,1]
heatmap6$SID<-NULL


pheatmap::pheatmap(heatmap6, cluster_rows = TRUE, cluster_cols = FALSE,clustering_method ='ward.D', show_rownames = FALSE, cellwidth = 15)

MiAMP_SpolyBP05<-subset(topBvsP05np$table, rownames(topBvsP05np$table) %in% MiAMP_Spoly$`Gene Name` )
MiAMP_SpolyBP05L<-MiAMP_SpolyBP05[,1]
MiAMP_SpolyBP05L<-data.matrix(MiAMP_SpolyBP05L)
rownames(MiAMP_SpolyBP05L)<-rownames(MiAMP_SpolyBP05)





write.table(topBvsP6$table, file="topTags-BvsP-6hr-Batch1-gm-LRT.csv", sep=',', quote=FALSE)
write.table(topBvsH6$table, file="topTags-BvsH-6hr-Batch1-gm-LRT.csv", sep=',', quote=FALSE)
write.table(topHvsP6$table, file="topTags-HvsP-6hr-Batch1-gm-LRT.csv", sep=',', quote=FALSE)

write.table(topBvsP6np$table, file="topTags-BvsP-6hr-Batch1-gm-LRTnp.csv", sep=',', quote=FALSE)
write.table(topBvsH6np$table, file="topTags-BvsH-6hr-Batch1-gm-LRTnp.csv", sep=',', quote=FALSE)
write.table(topHvsP6np$table, file="topTags-HvsP-6hr-Batch1-gm-LRTnp.csv", sep=',', quote=FALSE)
