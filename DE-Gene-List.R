jaz<-read.table('jazgid.csv', header = TRUE, sep=',')
wrky<-read.table('WRKYGID.csv', header = FALSE, sep=',')
miamp<-read.table('MIAMPGID.csv', header = FALSE, sep=',')
nlr<-read.table('NLRGID.csv', header = FALSE, sep=',')

wrky$V1<-paste(wrky$V1,'v2',sep='.')
miamp$V1<-paste(miamp$V1,'v2',sep='.')
nlr$V1<-paste(nlr$V1,'v2',sep='.')
jaz$JAZ.GID<-paste(jaz$JAZ.GID, 'v2', sep = ".")

JAZBH12<-subset(topBvsH12np$table, rownames(topBvsH12np$table) %in% jaz$JAZ.GID )
JAZBP12<-subset(topBvsP12np$table, rownames(topBvsP12np$table) %in% jaz$JAZ.GID )
JAZHP12<-subset(topHvsP12np$table, rownames(topHvsP12np$table) %in% jaz$JAZ.GID )

JAZBH6<-subset(topBvsH6np$table, rownames(topBvsH6np$table) %in% jaz$JAZ.GID )
JAZBP6<-subset(topBvsP6np$table, rownames(topBvsP6np$table) %in% jaz$JAZ.GID )
JAZHP6<-subset(topHvsP6np$table, rownames(topHvsP6np$table) %in% jaz$JAZ.GID )

JAZBH1<-subset(topBvsH1np$table, rownames(topBvsH1np$table) %in% jaz$JAZ.GID )
JAZBP1<-subset(topBvsP1np$table, rownames(topBvsP1np$table) %in% jaz$JAZ.GID )
JAZHP1<-subset(topHvsP1np$table, rownames(topHvsP1np$table) %in% jaz$JAZ.GID )

JAZBP05<-subset(topBvsP05np$table, rownames(topBvsP05np$table) %in% jaz$JAZ.GID )


WRKBH12<-subset(topBvsH12np$table, rownames(topBvsH12np$table) %in% wrky$V1)
WRKBP12<-subset(topBvsP12np$table, rownames(topBvsP12np$table) %in% wrky$V1)
WRKHP12<-subset(topHvsP12np$table, rownames(topHvsP12np$table) %in% wrky$V1)

WRKBH6<-subset(topBvsH6np$table, rownames(topBvsH6np$table) %in% wrky$V1)
WRKBP6<-subset(topBvsP6np$table, rownames(topBvsP6np$table) %in% wrky$V1)
WRKHP6<-subset(topHvsP6np$table, rownames(topHvsP6np$table) %in% wrky$V1)

WRKBH1<-subset(topBvsH1np$table, rownames(topBvsH1np$table) %in% wrky$V1)
WRKBP1<-subset(topBvsP1np$table, rownames(topBvsP1np$table) %in% wrky$V1)
WRKHP1<-subset(topHvsP1np$table, rownames(topHvsP1np$table) %in% wrky$V1)

WRKBP05<-subset(topBvsP05np$table, rownames(topBvsP05np$table) %in% wrky$V1 )


NLRBH12<-subset(topBvsH12np$table, rownames(topBvsH12np$table) %in% nlr$V1)
NLRBP12<-subset(topBvsP12np$table, rownames(topBvsP12np$table) %in% nlr$V1)
NLRHP12<-subset(topHvsP12np$table, rownames(topHvsP12np$table) %in% nlr$V1)

NLRBH6<-subset(topBvsH6np$table, rownames(topBvsH6np$table) %in% nlr$V1)
NLRBP6<-subset(topBvsP6np$table, rownames(topBvsP6np$table) %in% nlr$V1)
NLRHP6<-subset(topHvsP6np$table, rownames(topHvsP6np$table) %in% nlr$V1)

NLRBH1<-subset(topBvsH1np$table, rownames(topBvsH1np$table) %in% nlr$V1)
NLRBP1<-subset(topBvsP1np$table, rownames(topBvsP1np$table) %in% nlr$V1)
NLRHP1<-subset(topHvsP1np$table, rownames(topHvsP1np$table) %in% nlr$V1)

NLRBP05<-subset(topBvsP05np$table, rownames(topBvsP05np$table) %in% nlr$V1)



AMPBH12<-subset(topBvsH12np$table, rownames(topBvsH12np$table) %in% miamp$V1)
AMPBP12<-subset(topBvsP12np$table, rownames(topBvsP12np$table) %in% miamp$V1)
AMPHP12<-subset(topHvsP12np$table, rownames(topHvsP12np$table) %in% miamp$V1)

AMPBH6<-subset(topBvsH6np$table, rownames(topBvsH6np$table) %in% miamp$V1)
AMPBP6<-subset(topBvsP6np$table, rownames(topBvsP6np$table) %in% miamp$V1)
AMPHP6<-subset(topHvsP6np$table, rownames(topHvsP6np$table) %in% miamp$V1)

AMPBH1<-subset(topBvsH1np$table, rownames(topBvsH1np$table) %in% miamp$V1)
AMPBP1<-subset(topBvsP1np$table, rownames(topBvsP1np$table) %in% miamp$V1)
AMPHP1<-subset(topHvsP1np$table, rownames(topHvsP1np$table) %in% miamp$V1)

AMPBP05<-subset(topBvsP05np$table, rownames(topBvsP05np$table) %in% miamp$V1)


MARBH12<-subset(topBvsH12np$table, rownames(topBvsH12np$table) %in% CombiMarkers$X1)
MARBP12<-subset(topBvsP12np$table, rownames(topBvsP12np$table) %in% CombiMarkers$X1)
MARHP12<-subset(topHvsP12np$table, rownames(topHvsP12np$table) %in% CombiMarkers$X1)

MARBH6<-subset(topBvsH6np$table, rownames(topBvsH6np$table) %in% CombiMarkers$X1)
MARBP6<-subset(topBvsP6np$table, rownames(topBvsP6np$table) %in% CombiMarkers$X1)
MARHP6<-subset(topHvsP6np$table, rownames(topHvsP6np$table) %in% CombiMarkers$X1)

MARBH1<-subset(topBvsH1np$table, rownames(topBvsH1np$table) %in% CombiMarkers$X1)
MARBP1<-subset(topBvsP1np$table, rownames(topBvsP1np$table) %in% CombiMarkers$X1)
MARHP1<-subset(topHvsP1np$table, rownames(topHvsP1np$table) %in% CombiMarkers$X1)

MARBP05<-subset(topBvsP05np$table, rownames(topBvsP05np$table) %in% CombiMarkers$X1)

MARBP05p <- MARBP05$FDR <= 0.05
MARBP05p <- MARBP05[MARBP05p,]

MARBP1p <- MARBP1$FDR <= 0.05
MARBP1p <- MARBP1[MARBP1p,]
MARBH1p <- MARBH1$FDR <= 0.05
MARBH1p <- MARBH1[MARBH1p,]
MARHP1p <- MARHP1$FDR <= 0.05
MARHP1p <- MARHP1[MARHP1p,]


MARBP6p <- MARBP6$FDR <= 0.05
MARBP6p <- MARBP6[MARBP6p,]
MARBH6p <- MARBH6$FDR <= 0.05
MARBH6p <- MARBH6[MARBH6p,]
MARHP6p <- MARHP6$FDR <= 0.05
MARHP6p <- MARHP6[MARHP6p,]

MARBP12p <- MARBP12$FDR <= 0.05
MARBP12p <- MARBP12[MARBP12p,]
MARBH12p <- MARBH12$FDR <= 0.05
MARBH12p <- MARBH12[MARBH12p,]
MARHP12p <- MARHP12$FDR <= 0.05
MARHP12p <- MARHP12[MARHP12p,]

GISpo$V1<-paste(GISpo$V1,'v2', sep='.')
GIBH12<-subset(topBvsH12np$table, rownames(topBvsH12np$table) %in% GISpo$V1)
GIBP12<-subset(topBvsP12np$table, rownames(topBvsP12np$table) %in% GISpo$V1)
GIHP12<-subset(topHvsP12np$table, rownames(topHvsP12np$table) %in% GISpo$V1)

GIBH6<-subset(topBvsH6np$table, rownames(topBvsH6np$table) %in% GISpo$V1)
GIBP6<-subset(topBvsP6np$table, rownames(topBvsP6np$table) %in% GISpo$V1)
GIHP6<-subset(topHvsP6np$table, rownames(topHvsP6np$table) %in% GISpo$V1)

GIBH1<-subset(topBvsH1np$table, rownames(topBvsH1np$table) %in% GISpo$V1)
GIBP1<-subset(topBvsP1np$table, rownames(topBvsP1np$table) %in% GISpo$V1)
GIHP1<-subset(topHvsP1np$table, rownames(topHvsP1np$table) %in% GISpo$V1)

GIBP05<-subset(topBvsP05np$table, rownames(topBvsP05np$table) %in% GISpo$V1)

BP5 <- topBvsP05np$table$FDR <= 0.05 & abs(topBvsP05np$table$logFC) >= 1

BP1 <- topBvsP1np$table$FDR <= 0.05 & abs(topBvsP1np$table$logFC) >= 1
BH1 <- topBvsH1np$table$FDR <= 0.05 & abs(topBvsH1np$table$logFC) >= 1
HP1 <- topHvsP1np$table$FDR <= 0.05 & abs(topHvsP1np$table$logFC) >= 1

BP6 <- topBvsP6np$table$FDR <= 0.05 & abs(topBvsP6np$table$logFC) >= 1
BH6 <- topBvsH6np$table$FDR <= 0.05 & abs(topBvsH6np$table$logFC) >= 1
HP6 <- topHvsP6np$table$FDR <= 0.05 & abs(topHvsP6np$table$logFC) >= 1

BP12 <- topBvsP12np$table$FDR <= 0.05 & abs(topBvsP12np$table$logFC) >= 1
BH12 <- topBvsH12np$table$FDR <= 0.05 & abs(topBvsH12np$table$logFC) >= 1
HP12 <- topHvsP12np$table$FDR <= 0.05 & abs(topHvsP12np$table$logFC) >= 1







write.table(MARBH12, file="MARBH12.csv", sep=',', quote=FALSE)
write.table(MARBP12, file="MARBP12.csv", sep=',', quote=FALSE)
write.table(MARHP12, file="MARHP12.csv", sep=',', quote=FALSE)
write.table(MARBH6, file="MARBH6.csv", sep=',', quote=FALSE)
write.table(MARBP6, file="MARBP6.csv", sep=',', quote=FALSE)
write.table(MARHP6, file="MARHP6.csv", sep=',', quote=FALSE)
write.table(MARBH1, file="MARBH1.csv", sep=',', quote=FALSE)
write.table(MARBP1, file="MARBP1.csv", sep=',', quote=FALSE)
write.table(MARHP1, file="MARHP1.csv", sep=',', quote=FALSE)
write.table(MARBP05, file="MARBP05.csv", sep=',', quote=FALSE)

write.table(AMPBH12, file="AMPBH12.csv", sep=',', quote=FALSE)
write.table(AMPBP12, file="AMPBP12.csv", sep=',', quote=FALSE)
write.table(AMPHP12, file="AMPHP12.csv", sep=',', quote=FALSE)
write.table(AMPBH6, file="AMPBH6.csv", sep=',', quote=FALSE)
write.table(AMPBP6, file="AMPBP6.csv", sep=',', quote=FALSE)
write.table(AMPHP6, file="AMPHP6.csv", sep=',', quote=FALSE)
write.table(AMPBH1, file="AMPBH1.csv", sep=',', quote=FALSE)
write.table(AMPBP1, file="AMPBP1.csv", sep=',', quote=FALSE)
write.table(AMPHP1, file="AMPHP1.csv", sep=',', quote=FALSE)
write.table(AMPBP05, file="AMPBP05.csv", sep=',', quote=FALSE)

write.table(NLRBP05, file="NLRBP05.csv", sep=',', quote=FALSE)
write.table(NLRBP1, file="NLRBP1.csv", sep=',', quote=FALSE)
write.table(NLRBH1, file="NLRBH1.csv", sep=',', quote=FALSE)
write.table(NLRHP1, file="NLRHP1.csv", sep=',', quote=FALSE)
write.table(NLRBP6, file="NLRBP6.csv", sep=',', quote=FALSE)
write.table(NLRBH6, file="NLRBH6.csv", sep=',', quote=FALSE)
write.table(NLRHP6, file="NLRHP6.csv", sep=',', quote=FALSE)
write.table(NLRBP12, file="NLRBP12.csv", sep=',', quote=FALSE)
write.table(NLRBH12, file="NLRBH12.csv", sep=',', quote=FALSE)
write.table(NLRHP12, file="NLRHP12.csv", sep=',', quote=FALSE)

write.table(WRKHP12, file="WRKHP12.csv", sep=',', quote=FALSE)
write.table(WRKBP12, file="WRKBP12.csv", sep=',', quote=FALSE)
write.table(WRKBH12, file="WRKBH12.csv", sep=',', quote=FALSE)
write.table(WRKHP1, file="WRKHP1.csv", sep=',', quote=FALSE)
write.table(WRKBP1, file="WRKBP1.csv", sep=',', quote=FALSE)
write.table(WRKBH1, file="WRKBH1.csv", sep=',', quote=FALSE)
write.table(WRKHP6, file="WRKHP6.csv", sep=',', quote=FALSE)
write.table(WRKBP6, file="WRKBP6.csv", sep=',', quote=FALSE)
write.table(WRKBH6, file="WRKBH6.csv", sep=',', quote=FALSE)
write.table(WRKBP05, file="WRKBP05.csv", sep=',', quote=FALSE)

write.table(JAZBP05, file="JAZBP05.csv", sep=',', quote=FALSE)
write.table(JAZBP1, file="JAZBP1.csv", sep=',', quote=FALSE)
write.table(JAZBH1, file="JAZBH1.csv", sep=',', quote=FALSE)
write.table(JAZHP1, file="JAZHP1.csv", sep=',', quote=FALSE)
write.table(JAZBP6, file="JAZBP6.csv", sep=',', quote=FALSE)
write.table(JAZBH6, file="JAZBH6.csv", sep=',', quote=FALSE)
write.table(JAZHP6, file="JAZHP6.csv", sep=',', quote=FALSE)
write.table(JAZBP12, file="JAZBP12.csv", sep=',', quote=FALSE)
write.table(JAZBH12, file="JAZBH12.csv", sep=',', quote=FALSE)
write.table(JAZHP12, file="JAZHP12.csv", sep=',', quote=FALSE)

