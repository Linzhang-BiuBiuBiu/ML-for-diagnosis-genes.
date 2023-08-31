library(pROC)                  
expFile="diffGeneExp.txt"      
geneFile="interGenes.txt"      
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))   

rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
y=gsub("(.*)\\_(.*)", "\\2", colnames(rt))
y=ifelse(y=="con", 0, 1)

geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)

train_roc <- data.frame()

for(x in as.vector(geneRT[,1])){
	roc1=roc(y, as.numeric(rt[x,]))
	ci1=ci.auc(roc1, method="bootstrap")
	ciVec=as.numeric(ci1)
	pdf(file=paste0("train_ROC.",x,".pdf"), width=5, height=5)
	plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=x)
	text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"−",sprintf("%.03f",ciVec[3])), col="red")
	dev.off()
	train_r <- cbind(variable=x,value=ciVec[2])
	train_roc <- rbind(train_r,train_roc)
}
write.table(train_roc, file="train_roc.txt", sep="\t", quote=F, row.names=T, col.names=T)