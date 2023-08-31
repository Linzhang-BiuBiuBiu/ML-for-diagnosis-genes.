inputFile="normalize_train.txt"      
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))   
source("Cibersort.R")       

outTab=CIBERSORT("ref.txt", inputFile, perm=1000, QN=TRUE)

outTab=outTab[outTab[,"P-value"]<0.05,]
outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])
outTab=rbind(id=colnames(outTab),outTab)
write.table(outTab, file="CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)
