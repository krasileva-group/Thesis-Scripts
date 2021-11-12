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

readCount1 <- read.table("gm_1h-BP-B19-B20-P21.csv", header=TRUE, sep=",", row.names=1)
Treat1 <- factor(substring(colnames(readCount1),1,2))
Treat1<-relevel(Treat1,ref = "B.")

y1<-DGEList(counts = readCount1, group=Treat1)
keep1 <- filterByExpr(y1)
y1 <- y1[keep1, ,keep.lib.sizes=FALSE]
y1 <- calcNormFactors(y1)


design1 <- model.matrix(~0+Treat1)


rownames(design1)<-colnames(y1)
y1 <- estimateDisp(y1, design1, robust=TRUE)

fit1L<- glmFit(y1, design1)
CONTRASTS1L <- makeContrasts(BvsP = Treat1P.-Treat1B.,
                             BvsH = Treat1H.-Treat1B.,
                             HvsP = Treat1P.-Treat1H.,
                             levels = design1)

GLM1L<- paste(colnames(CONTRASTS1L),sep="")
for (i in 1:ncol(CONTRASTS1L)){
  current.glmQLFTest <- glmLRT(fit1L, contrast = CONTRASTS1L[,i])
  outputFile1 <- paste0((colnames(CONTRASTS1L)[i]), "1hr.csv")
  write.table(current.glmQLFTest, file=outputFile1, sep=',', quote=FALSE)
  assign(GLM1L[i], current.glmQLFTest)}

topBvsP1<- topTags(BvsP, n=nrow(BvsP$table), adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
topBvsH1<- topTags(BvsH, n=nrow(BvsH$table), adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
topHvsP1<- topTags(HvsP, n=nrow(HvsP$table), adjust.method = "BH", sort.by = "PValue", p.value = 0.05)


topBvsP1np<- topTags(BvsP, n=nrow(BvsP$table), adjust.method = "BH", sort.by = "PValue")
topBvsH1np<- topTags(BvsH, n=nrow(BvsH$table), adjust.method = "BH", sort.by = "PValue")
topHvsP1np<- topTags(HvsP, n=nrow(HvsP$table), adjust.method = "BH", sort.by = "PValue")

library(readr)
CombiMarkers <- read_table2("CombiMarkers.csv",
                            col_names = FALSE)
View(CombiMarkers)
CombiMarkers$X1<-paste(CombiMarkers$X1,"v2" ,sep=".")


CombiMarkersBP1<-subset(topBvsP1np$table, rownames(topBvsP1np$table) %in% CombiMarkers$X1 )
CombiMarkersBP1L<-CombiMarkersBP1[,1]
CombiMarkersBH1<-subset(topBvsH1np$table, rownames(topBvsH1np$table) %in% CombiMarkers$X1 )
CombiMarkersBH1L<-CombiMarkersBH1[,1]
CombiMarkersHP1<-subset(topHvsP1np$table, rownames(topHvsP1np$table) %in% CombiMarkers$X1 )
CombiMarkersHP1L<-CombiMarkersHP1[,1]



CombiMarkersBP1L<-data.matrix(CombiMarkersBP1L)
rownames(CombiMarkersBP1L)<-rownames(CombiMarkersBP1)


pheatmap::pheatmap(CombiMarkersBP1L, cluster_rows = TRUE, cluster_cols = FALSE,clustering_method ='ward.D', show_rownames = FALSE, cellwidth = 15)

MiAMP_SpolyBP05<-subset(topBvsP05np$table, rownames(topBvsP05np$table) %in% MiAMP_Spoly$`Gene Name` )
MiAMP_SpolyBP05L<-MiAMP_SpolyBP05[,1]
MiAMP_SpolyBP05L<-data.matrix(MiAMP_SpolyBP05L)
rownames(MiAMP_SpolyBP05L)<-rownames(MiAMP_SpolyBP05)




write.table(topBvsP1, file="topTags-BvsP-1hr-Batch1-gm-LRT.csv", sep=',', quote=FALSE)
write.table(topBvsH1, file="topTags-BvsH-1hr-Batch1-gm-LRT.csv", sep=',', quote=FALSE)
write.table(topHvsP1, file="topTags-HvsP-1hr-Batch1-gm-LRT.csv", sep=',', quote=FALSE)
write.table(topBvsP1np, file="topTags-BvsP-1hr-Batch1-gm-LRTnp.csv", sep=',', quote=FALSE)
write.table(topBvsH1np, file="topTags-BvsH-1hr-Batch1-gm-LRTnp.csv", sep=',', quote=FALSE)
write.table(topHvsP1np, file="topTags-HvsP-1hr-Batch1-gm-LRTnp.csv", sep=',', quote=FALSE)

            
                        
