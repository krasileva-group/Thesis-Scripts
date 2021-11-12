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
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)

readCount05 <- read.table("gm_05h_B1.csv", header=TRUE, sep=",", row.names=1)
Treat05 <- factor(substring(colnames(readCount05),1,2))
Treat05<-relevel(Treat05,ref = "B.")

y05<-DGEList(counts = readCount05, group=Treat05)
keep05 <- filterByExpr(y05)
y05 <- y05[keep05, keep.lib.sizes=FALSE]
y05 <- calcNormFactors(y05)


design05 <- model.matrix(~0+Treat05)


rownames(design05)<-colnames(y05)
y05 <- estimateDisp(y05, design05, robust=TRUE)

et<-exactTest(y05, pair=c("B.","P."))
topBvsP05np<- topTags(et, n=nrow(y05$counts), adjust.method = "BH", sort.by = "PValue")
topBvsP05<- topTags(et, n=nrow(y05$counts), adjust.method = "BH", sort.by = "PValue", p.value = 0.05)

topBP05F<-topBvsP05$table[keeps,]


topBP05FL<-topBP05F[,1]
topBP05FL<-data.matrix(topBP05FL)
rownames(topBP05FL)<-rownames(topBP05F$table)
breaksList = seq(-2, 2, by = 0.5)

pheatmap::pheatmap(topBP05FL, cluster_rows = TRUE, cluster_cols = FALSE,clustering_method ='ward.D', color=viridis::cividis(10), cellwidth = 15, breaks = breaksList)


library(readr)
MiAMP_Spoly <- read_csv("/Volumes/Orthofinder-KVK/RNAseq-Duckweed/analysis/EdgeR/batch-corrected/MiAMP-Spoly.csv")
marker <- read_table2("Marker1.csv", col_names = FALSE)
celldeath <- read_table2("celldeath.csv", col_names = FALSE)
Zipfel_marker <- read_table2("Zipfel-marker.csv", col_names = FALSE)
Marta_core_PTI_Genes <- read_table2("Marta-core-PTI-Genes.csv",col_names = FALSE)

celldeath$X1<-paste(celldeath$X1,"v2" ,sep=".")
Zipfel_marker$X1<-paste(Zipfel_marker$X1,"v2" ,sep=".")
Marta_core_PTI_Genes$X1<-paste(Marta_core_PTI_Genes$X1,"v2" ,sep=".")
MiAMP_Spoly$`Gene Name`<-paste(MiAMP_Spoly$`Gene Name`, 'v2', sep = ".")
marker$X1<-paste(marker$X1,"v2" ,sep=".")

MiAMP_SpolyBP05<-subset(topBvsP05np$table, rownames(topBvsP05np$table) %in% MiAMP_Spoly$`Gene Name` )
MiAMP_SpolyBP05L<-MiAMP_SpolyBP05[,1]
MiAMP_SpolyBP05L<-data.matrix(MiAMP_SpolyBP05L)
rownames(MiAMP_SpolyBP05L)<-rownames(MiAMP_SpolyBP05)

MarkerBP05<-subset(topBvsP05np$table, rownames(topBvsP05np$table) %in% dataset$X1 )
MarkerBP05L<-MarkerBP05[,1]
MarkerBP05L<-data.matrix(MarkerBP05L)
rownames(MarkerBP05L)<-rownames(MarkerBP05)



celldeathBP05<-subset(topBvsP05np$table, rownames(topBvsP05np$table) %in% celldeath$X1)
Zipfel_markerBP05<-subset(topBvsP05np$table, rownames(topBvsP05np$table) %in% Zipfel_marker$X1)
corePTIBP05<-subset(topBvsP05np$table, rownames(topBvsP05np$table) %in% Marta_core_PTI_Genes$X1)

corePTIBP05L<-corePTIBP05[,1]
corePTIBP05L<-data.matrix(corePTIBP05L)
rownames(corePTIBP05L)<-rownames(corePTIBP05)

celldeathBP05L<-celldeathBP05[,1]
celldeathBP05L<-data.matrix(celldeathBP05L)
rownames(celldeathBP05L)<-rownames(celldeathBP05)

Zipfel_markerBP05L<-Zipfel_markerBP05[,1]
Zipfel_markerBP05L<-data.matrix(Zipfel_markerBP05L)
rownames(Zipfel_markerBP05L)<-rownames(Zipfel_markerBP05)

library(readr)
CombiMarkers <- read_table2("CombiMarkers.csv",
col_names = FALSE)
View(CombiMarkers)
CombiMarkers$X1<-paste(CombiMarkers$X1,"v2" ,sep=".")
CombiMarkersBP05<-subset(topBvsP05np$table, rownames(topBvsP05np$table) %in% CombiMarkers$X1 )
CombiMarkersBP05L<-CombiMarkersBP05[,1]
CombiMarkersBP05L<-data.matrix(CombiMarkersBP05L)
rownames(CombiMarkersBP05L)<-rownames(CombiMarkersBP05)

pheatmap::pheatmap(corePTIBP05L, cluster_rows = TRUE, cluster_cols = FALSE, color=viridis::cividis(length(breaksList)), cellwidth = 15, breaks = breaksList)
pheatmap::pheatmap(celldeathBP05L, cluster_rows = TRUE, cluster_cols = FALSE,clustering_method ='ward.D', color=viridis::cividis(length(breaksList)), cellwidth = 15, breaks = breaksList)
pheatmap::pheatmap(Zipfel_markerBP05L, cluster_rows = TRUE, cluster_cols = FALSE,clustering_method ='ward.D', color=viridis::cividis(length(breaksList)), cellwidth = 15, breaks = breaksList)
pheatmap::pheatmap(MarkerBP05L, cluster_rows = TRUE, cluster_cols = FALSE,clustering_method ='ward.D', color=viridis::cividis(length(breaksList)), cellwidth = 15, breaks = breaksList)
pheatmap::pheatmap(CombiMarkersBP05L, cluster_rows = TRUE, cluster_cols = FALSE,clustering_method ='ward.D', show_rownames = FALSE, cellwidth = 15)




write.table(topBvsP05np$table, file="topTags-HvsP-05hr-Batch-exactTest-noP-filter.csv", sep=',', quote=FALSE)
write.table(topBP05F, file="topTags-HvsP-05hr-Batch-exactTest-FDR05-LogFC1.csv", sep=',', quote=FALSE)
write.table(MiAMP_SpolyBP05, file="topTags-HvsP-05hr-Batch-exactTest-MiAMP.csv", sep=',', quote=FALSE)
write.table(Zipfel_markerBP05, file="topTags-HvsP-05hr-Batch-exactTest-Zipfel.csv", sep=',', quote=FALSE)
write.table(celldeathBP05, file="topTags-HvsP-05hr-Batch-exactTest-celldeath.csv", sep=',', quote=FALSE)
write.table(corePTIBP05, file="topTags-HvsP-05hr-Batch-exactTest-corePTI.csv", sep=',', quote=FALSE)
write.table(CombiMarkersBP05, file="topTags-HvsP-05hr-Batch-exactTest-combiMarkers.csv", sep=',', quote=FALSE)
#####
keeps <- topBvsP05np$table$FDR <= 0.05 & abs(topBvsP05np$table$logFC) >= 1


logCounts05 <- cpm(y05,log=TRUE)
logCounts05cPTI<-logCounts[corePTIBP05,]

