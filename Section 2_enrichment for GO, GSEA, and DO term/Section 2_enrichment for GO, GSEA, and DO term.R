
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("DOSE")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pvalueFilter=0.05       
qvalueFilter=0.05       


colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))      
rt=read.table("diff.txt", header=T, sep="\t", check.names=F)     

genes=as.vector(rt[,1])

entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]      

#GO term
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)

# DO term
kk=enrichDO(gene=gene, ont="DO", pvalueCutoff=1, qvalueCutoff=1, readable=TRUE)
DO=as.data.frame(kk)
DO=DO[(DO$pvalue<pvalueFilter & DO$qvalue<qvalueFilter),]
write.table(DO, file="DO.txt", sep="\t", quote=F, row.names = F)

#GSEA term
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

inputFile="all.txt"        
gmtFile="GSEA.gmt"      
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))      #???ù???Ŀ¼

rt=read.table(inputFile, header=T, sep="\t", check.names=F)
logFC=as.vector(rt[,2])
names(logFC)=as.vector(rt[,1])
logFC=sort(logFC, decreasing=T)

gmt=read.gmt(gmtFile)

#GSEA term
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result.txt",sep="\t",quote=F,row.names = F)

