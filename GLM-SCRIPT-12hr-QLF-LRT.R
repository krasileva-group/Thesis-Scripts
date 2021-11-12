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

readCount12 <- read.table("gm_12h_B1.csv", header=TRUE, sep=",", row.names=1)
Treat12 <- factor(substring(colnames(readCount12),1,2))
Treat12<-relevel(Treat12,ref = "B.")
     

y12<-DGEList(counts = readCount12, group=Treat12)
keep12 <- filterByExpr(y12)
y12 <- y12[keep12, ,keep.lib.sizes=FALSE]
y12 <- calcNormFactors(y12)


design12 <- model.matrix(~0+Treat12)


rownames(design12)<-colnames(y12)
y12 <- estimateDisp(y12, design12, robust=TRUE)


fit12L<- glmFit(y12, design12)
CONTRASTS12 <- makeContrasts( BvsP12 = Treat12P.-Treat12B.,
                              BvsH12 = Treat12H.-Treat12B.,
                              HvsP12 = Treat12P.-Treat12H.,
                              levels = design12)

GLM12L<- paste(colnames(CONTRASTS12),sep="")
for (i in 1:ncol(CONTRASTS12)){
  current.glmQLFTest <- glmLRT(fit12L, contrast = CONTRASTS12[,i])
  outputFile12 <- paste0((colnames(CONTRASTS12)[i]), "12hr.csv")
  write.table(current.glmQLFTest, file=outputFile12, sep=',', quote=FALSE)
  assign(GLM12L[i], current.glmQLFTest)}

topBvsP12<- topTags(BvsP12, n=nrow(BvsP12$table), adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
topBvsH12<- topTags(BvsH12, n=nrow(BvsH12$table), adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
topHvsP12<- topTags(HvsP12, n=nrow(HvsP12$table), adjust.method = "BH", sort.by = "PValue", p.value = 0.05)


topBvsP12np<- topTags(BvsP12, n=nrow(BvsP12$table), adjust.method = "BH", sort.by = "PValue")
topBvsH12np<- topTags(BvsH12, n=nrow(BvsH12$table), adjust.method = "BH", sort.by = "PValue")
topHvsP12np<- topTags(HvsP12, n=nrow(HvsP12$table), adjust.method = "BH", sort.by = "PValue")


write.table(topBvsP12$table, file="topTags-BvsP-12hr-Batch1-gm-LRT.csv", sep=',', quote=FALSE)
write.table(topBvsH12$table, file="topTags-BvsH-12hr-Batch1-gm-LRT.csv", sep=',', quote=FALSE)
write.table(topHvsP12$table, file="topTags-HvsP-12hr-Batch1-gm-LRT.csv", sep=',', quote=FALSE)

write.table(topBvsP12np$table, file="topTags-BvsP-12hr-Batch1-gm-LRTnp.csv", sep=',', quote=FALSE)
write.table(topBvsH12np$table, file="topTags-BvsH-12hr-Batch1-gm-LRTnp.csv", sep=',', quote=FALSE)
write.table(topHvsP12np$table, file="topTags-HvsP-12hr-Batch1-gm-LRTnp.csv", sep=',', quote=FALSE)





BHdown12 <- topBvsH12$table$FDR <= 0.05 & topBvsH12$table$logFC <= -1
BHup12 <- topBvsH12$table$FDR <= 0.05 & topBvsH12$table$logFC >= 1
BPup12 <- topBvsP12$table$FDR <= 0.05 & topBvsP12$table$logFC >= 1
BPdown12 <- topBvsP12$table$FDR <= 0.05 & topBvsP12$table$logFC <= -1
HPdown12 <- topHvsP12$table$FDR <= 0.05 & topHvsP12$table$logFC <= -1
HPup12 <- topHvsP12$table$FDR <= 0.05 & topHvsP12$table$logFC >= 1
downBH12<-topBvsH12$table[BHdown12,]
upBH12<-topBvsH12$table[BHup12,]
upBP12<-topBvsP12$table[BPup12,]
downBP12<-topBvsP12$table[BPdown12,]
downHP12<-topHvsP12$table[HPdown12,]
upHP12<-topHvsP12$table[HPup12,]
write.table(upBH12, file="topTags-BH12-UP-Lfc1-FDR0-5-H116.csv", sep=',', quote=FALSE)
write.table(downBH12, file="topTags-BH12-DOWN-Lfc1-FDR0-5-H116.csv", sep=',', quote=FALSE)
write.table(downBP12, file="topTags-BP12-DOWN-Lfc1-FDR0-5-H116.csv", sep=',', quote=FALSE)
write.table(upBP12, file="topTags-BP12-UP-Lfc1-FDR0-5-H116.csv", sep=',', quote=FALSE)
write.table(upHP12, file="topTags-HP12-UP-Lfc1-FDR0-5-H116.csv", sep=',', quote=FALSE)
write.table(downHP12, file="topTags-HP12-DOWN-Lfc1-FDR0-5-H116.csv", sep=',', quote=FALSE)


volcanoData12 <- cbind(topBvsP12np$table$logFC, -log10(topBvsP12np$table$FDR))
colnames(volcanoData12) <- c("logFC", "negLogPval")
