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

readCount05 <- read.table("../ReadCountAdjustedTreatTime-0-5hr-B3P9.csv", header=TRUE, sep=",", row.names=1)
Treat05 <- factor(substring(colnames(readCount05),1,2))
Treat05<-relevel(Treat05,ref = "B.")
batch05<-factor(c(rep(1,3),rep(2,2),rep(1,2),rep(2,5),rep(1,1),rep(2,1),rep(1,1)))

y05<-DGEList(counts = readCount05, group=Treat05)
keep05 <- filterByExpr(y05)
y05 <- y05[keep05, keep.lib.sizes=FALSE]
y05 <- calcNormFactors(y05)


design05 <- model.matrix(~0+Treat05+batch05)


rownames(design05)<-colnames(y05)
y05 <- estimateDisp(y05, design05, robust=TRUE)

fit05 <- glmQLFit(y05, design05)

CONTRASTS05 <- makeContrasts( BvsP05Q = Treat05P.-Treat05B.,
                              BvsH05Q = Treat05H.-Treat05B.,
                              HvsP05Q = Treat05P.-Treat05H.,
                              levels = design05)

GLM05<- paste(colnames(CONTRASTS05),sep="")


for (i in 1:ncol(CONTRASTS05)){
  current.glmQLFTest <- glmQLFTest(fit05, contrast = CONTRASTS05[,i])
  outputFile05 <- paste0((colnames(CONTRASTS05)[i]), "05hr.csv")
  write.table(current.glmQLFTest, file=outputFile05, sep=',', quote=FALSE)
  assign(GLM05[i], current.glmQLFTest)}

topBvsP05Q<- topTags(BvsP05Q, n=nrow(keep05), adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
topBvsH05Q<- topTags(BvsH05Q, n=nrow(keep05), adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
topHvsP05Q<- topTags(HvsP05Q, n=nrow(keep05), adjust.method = "BH", sort.by = "PValue", p.value = 0.05)




fit05L<- glmFit(y05, design05)
CONTRASTS05L <- makeContrasts( BvsPL05 = Treat05P.-Treat05B.,
                               BvsHL05 = Treat05H.-Treat05B.,
                               HvsPL05 = Treat05P.-Treat05H.,
                               levels = design05)

GLM05L<- paste(colnames(CONTRASTS05L),sep="")
for (i in 1:ncol(CONTRASTS05)){
  current.glmQLFTest <- glmLRT(fit05L, contrast = CONTRASTS05L[,i])
  outputFile05 <- paste0((colnames(CONTRASTS05L)[i]), "05hr.csv")
  write.table(current.glmQLFTest, file=outputFile05, sep=',', quote=FALSE)
  assign(GLM05L[i], current.glmQLFTest)}

topBvsP05L<- topTags(BvsPL05, n=nrow(keep05), adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
topBvsH05L<- topTags(BvsHL05, n=nrow(keep05), adjust.method = "BH", sort.by = "PValue", p.value = 0.05)
topHvsP05L<- topTags(HvsPL05, n=nrow(keep05), adjust.method = "BH", sort.by = "PValue", p.value = 0.05)


topBvsP05Lnp<- topTags(BvsPL05, n=nrow(keep05), adjust.method = "BH", sort.by = "PValue")
topBvsH05Lnp<- topTags(BvsHL05, n=nrow(keep05), adjust.method = "BH", sort.by = "PValue")
topHvsP05Lnp<- topTags(HvsPL05, n=nrow(keep05), adjust.method = "BH", sort.by = "PValue")