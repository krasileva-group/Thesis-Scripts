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

readCount05 <- read.table("../ReadCountAdjustedTreatTime-0-5hr-B1-noHrcC.csv", header=TRUE, sep=",", row.names=1)
Treat05 <- factor(substring(colnames(readCount05),1,2))
Treat05<-relevel(Treat05,ref = "B.")

y05<-DGEList(counts = readCount05, group=Treat05)
keep05 <- filterByExpr(y05)
y05 <- y05[keep05, keep.lib.sizes=FALSE]
y05 <- calcNormFactors(y05)


design05 <- model.matrix(~0+Treat05)


rownames(design05)<-colnames(y05)
y05 <- estimateDisp(y05, design05, robust=TRUE)

fit05 <- glmFit(y05, design05)

CONTRASTS05 <- makeContrasts( BvsP05 = Treat05P.-Treat05B.,
                              BvsH05 = Treat05H.-Treat05B.,
                              HvsP05 = Treat05P.-Treat05H.,
                              levels = design05)

GLM05<- paste(colnames(CONTRASTS05),sep="")


for (i in 1:ncol(CONTRASTS05)){
  current.glmQLFTest <- glmLRT(fit05, contrast = CONTRASTS05[,i])
  outputFile05 <- paste0((colnames(CONTRASTS05)[i]), "05hr.csv")
  write.table(current.glmQLFTest, file=outputFile05, sep=',', quote=FALSE)
  assign(GLM05[i], current.glmQLFTest)}

topBvsP05<- topTags(BvsP05, n=nrow(keep05), adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
topBvsH05<- topTags(BvsH05, n=nrow(keep05), adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
topHvsP05<- topTags(HvsP05, n=nrow(keep05), adjust.method = "BH", sort.by = "PValue", p.value = 0.05)

topBvsP05np<- topTags(BvsP05, n=nrow(keep05), adjust.method = "BH", sort.by = "PValue")
topBvsH05np<- topTags(BvsH05, n=nrow(keep05), adjust.method = "BH", sort.by = "PValue")
topHvsP05np<- topTags(HvsP05, n=nrow(keep05), adjust.method = "BH", sort.by = "PValue")