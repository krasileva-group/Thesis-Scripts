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

readCount1 <- read.table("../ReadCountAdjustedTreatTime-1hr.csv", header=TRUE, sep=",", row.names=1)
Treat1 <- factor(substring(colnames(readCount1),1,2))
Treat1<-relevel(Treat1,ref = "B.")
batch1<-factor(c(rep(1,5),rep(2,7),rep(1,10),rep(2,2)))

y1<-DGEList(counts = readCount1, group=Treat1)
keep1 <- filterByExpr(y1)
y1 <- y1[keep1, keep.lib.sizes=FALSE]
y1 <- calcNormFactors(y1)


design1 <- model.matrix(~0+Treat1+batch1)


rownames(design1)<-colnames(y1)
y1 <- estimateDisp(y1, design1, robust=TRUE)

fit1 <- glmQLFit(y1, design1)

CONTRASTS1 <- makeContrasts( BvsPQ = Treat1P.-Treat1B.,
                            BvsHQ = Treat1H.-Treat1B.,
                            HvsPQ = Treat1P.-Treat1H.,
                            levels = design1)

GLM1<- paste(colnames(CONTRASTS1),sep="")


for (i in 1:ncol(CONTRASTS1)){
current.glmQLFTest <- glmQLFTest(fit1, contrast = CONTRASTS1[,i])
outputFile1 <- paste0((colnames(CONTRASTS1)[i]), "1hr.csv")
write.table(current.glmQLFTest, file=outputFile1, sep=',', quote=FALSE)
assign(GLM1[i], current.glmQLFTest)}

topBvsP1Q<- topTags(BvsPQ, n=nrow(keep1), adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
topBvsH1Q<- topTags(BvsHQ, n=nrow(keep1), adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
topHvsP1Q<- topTags(HvsPQ, n=nrow(keep1), adjust.method = "BH", sort.by = "PValue", p.value = 0.05)




fit1L<- glmFit(y1, design1)
CONTRASTS1L <- makeContrasts( BvsPL = Treat1P.-Treat1B.,
                             BvsHL = Treat1H.-Treat1B.,
                             HvsPL = Treat1P.-Treat1H.,
                             levels = design1)

GLM1L<- paste(colnames(CONTRASTS1L),sep="")
for (i in 1:ncol(CONTRASTS1)){
  current.glmQLFTest <- glmLRT(fit1L, contrast = CONTRASTS1L[,i])
  outputFile1 <- paste0((colnames(CONTRASTS1L)[i]), "1hr.csv")
  write.table(current.glmQLFTest, file=outputFile1, sep=',', quote=FALSE)
  assign(GLM1L[i], current.glmQLFTest)}

topBvsP1L<- topTags(BvsPL, n=nrow(keep1), adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
topBvsH1L<- topTags(BvsHL, n=nrow(keep1), adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
topHvsP1L<- topTags(HvsPL, n=nrow(keep1), adjust.method = "BH", sort.by = "PValue", p.value = 0.05)


topBvsP1Lnp<- topTags(BvsPL, n=nrow(keep1), adjust.method = "BH", sort.by = "PValue")
topBvsH1Lnp<- topTags(BvsHL, n=nrow(keep1), adjust.method = "BH", sort.by = "PValue")
topHvsP1Lnp<- topTags(HvsPL, n=nrow(keep1), adjust.method = "BH", sort.by = "PValue")

write.table(topBvsP1L, file="topTags-BvsP-1hr-Batch-glmLRT.csv", sep=',', quote=FALSE)
write.table(topBvsH1L, file="topTags-BvsH-1hr-Batch-glmLRT.csv", sep=',', quote=FALSE)
write.table(topHvsP1L, file="topTags-HvsP-1hr-Batch-glmLRT.csv", sep=',', quote=FALSE)
write.table(topBvsP1Q, file="topTags-BvsP-1hr-Batch-glmQLF.csv", sep=',', quote=FALSE)
write.table(topBvsH1Q, file="topTags-BvsH-1hr-Batch-glmQLF.csv", sep=',', quote=FALSE)
write.table(topHvsP1Q, file="topTags-HvsP-1hr-Batch-glmQLF.csv", sep=',', quote=FALSE)
            
                        
