#### download the GEO data
rm(list = ls())

#working condition

options( 'download.file.method.GEOquery' = 'libcurl' )
options('GEOquery.inmemory.gpl'=FALSE)
options(stringsAsFactors = F)

library(GEOquery)
library(Biobase)
library(limma)
library(stringi)
library(stringr)
library(openxlsx)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

rt <- data.frame("GSE59867","GSE60993",
                    "GSE62646","GSE48060")

 for(i in rt[,1]){
  options(digits=12)
  gset <- getGEO(i,destdir = ".",
                 AnnotGPL = T,
                 getGPL = T)

  expression_data <- exprs(gset[[1]])
  clinical_data <- pData(gset[[1]])
  
  index_data <- c("ID","Gene symbol","Gene Symbol","GENE_SYMBOL","Symbol")
  f_Data <- fData(gset[[1]])
  intersect_data <- intersect(index_data,colnames(fData(gset[[1]])))
  
  if(length(intersect_data)>1){
  index_prob <- colnames(fData(gset[[1]]))
  probe_data <- fData(gset[[1]])[,intersect_data]

  expression_prob <- cbind(ID = rownames(expression_data),expression_data)
  
  prob_gene_expression <- merge(probe_data,expression_prob,
                                by="ID")
  
  prob_gene_expression[prob_gene_expression==""] <- NA
  
  prob_gene_expression <- na.omit(prob_gene_expression)
  Gene_names <- gsub("(.*)//(.*)","\\2",prob_gene_expression[,2])
  
  prob_gene_expression[,2] <- Gene_names
  
  prob_gene_expression <- avereps(prob_gene_expression, ID=prob_gene_expression[,2])
  
  rownames(prob_gene_expression) <- prob_gene_expression[,2]
  gene_expression <- prob_gene_expression[,c(-1,-2)]
  gene_expression=avereps(gene_expression)
  
  clinialNames <- paste0(i,"_clinical.xlsx") 
  expNames <- paste0(i,"_expression.txt")

  write.xlsx(clinical_data,file = clinialNames)
  write.table(gene_expression,file=expNames,sep = "\t",col.names = T,
              row.names = T,quote=F)
  }
  
  }

####integration of data 
library(limma)
library(sva)
outFile="merge.txt"       

files=grep("expression.txt$", dir(), value=T)
geneList=list()

for(file in files){
    if(file==outFile){next}
    rt=read.table(file, header=T, sep="\t", check.names=F) 
    rt <- rt[-1,]
    ROWnames <- gsub("(.*)//(.*)","\\2",rownames(rt))
    data <- cbind(ID=ROWnames,rt)
    data=avereps(data,ID=data[,1])
    rownames(data) <- data[,1]
    data <- data[,-1]
    geneNames=rownames(data)     
    uniqGene=unique(geneNames)      
    header=unlist(strsplit(file, "\\.|\\-"))
    geneList[[header[1]]]=uniqGene
}

interGenes=Reduce(intersect, geneList)

allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
    inputFile=files[i]
    header=unlist(strsplit(inputFile, "\\.|\\-"))

        rt=read.table(inputFile, header=T, sep="\t", check.names=F)
    rt <- rt[-1,]
    ROWnames <- gsub("(.*)//(.*)","\\2",rownames(rt))
    data <- cbind(ID=ROWnames,rt)
    data=avereps(data,ID=data[,1])
    rownames(data) <- data[,1]
    data <- data[,-1]
    exp=data[,1:ncol(data)]
    dimnames=list(rownames(exp),colnames(exp))
    data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
    rt=avereps(data)
    
    qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
    if(LogC){
        rt[rt<0]=0
        rt=log2(rt+1)}
    rt=normalizeBetweenArrays(rt)
    
    if(i==1){
        intersectGENESS <- intersect(interGenes,rownames(rt))
        allTab=rt[intersectGENESS,]
    }else{
        intersectGENESS <- intersect(interGenes,rownames(rt))
        allTab=cbind(allTab, rt[intersectGENESS,])
    }
    batchType=c(batchType, rep(i,ncol(rt)))
}

#merge data and eliminate the branching effect
outTab=ComBat(allTab, batchType, par.prior=TRUE)
outTab=rbind(geneNames=colnames(outTab), outTab)
rownames(outTab) <- gsub("-","_",rownames(outTab))
write.table(outTab, file="merge.txt", sep="\t", quote=F, col.names=)

####explore the DEGs
library(limma)
library(pheatmap)
library(tidyverse)

inputFile="merge.txt"      
logFCfilter=0.7
adj.P.Val.Filter=0.05
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))      

rt=read.table(inputFile, header=T, sep="\t", check.names=T)
exp <- rt[rownames(rt)!="geneNames",]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

##define the colnames of control 
sampleName1=c()
files=dir()
files=grep("s1.txt$", files, value=T)
for(file in files){
    rt=read.table(file, header=F, sep="\t", check.names=F)      
    geneNames=as.vector(rt[,1])      
    uniqGene=unique(geneNames)       
    sampleName1=c(sampleName1, uniqGene)}

##define the colnames of AMI 
sampleName2=c()
files=dir()
files=grep("s2.txt$", files, value=T)
for(file in files){
    rt=read.table(file, header=F, sep="\t", check.names=F)     
    geneNames=as.vector(rt[,1])     
    uniqGene=unique(geneNames)      
    sampleName2=c(sampleName2, uniqGene)
}

##merge the control and AMI
conData=data[,sampleName1]
AMIData=data[,sampleName2]
data=cbind(conData,DCMData)
conNum=ncol(conData)
AMINum=ncol(AMIData)

##difference between AMI and control
Type=c(rep("con",conNum),rep("AMI",DCMNum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","AMI")
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(DCM-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
allDiffOut=rbind(id=colnames(allDiff),allDiff)
write.table(allDiffOut, file="all.txt", sep="\t", quote=F, col.names=F)

##save file with |LogFC|>1 and adj.P.Val<0.01
diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adj.P.Val.Filter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, file="diff.txt", sep="\t", quote=F, col.names=F)

outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData, file="normalize.txt", sep="\t", quote=F, col.names=F)

diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adj.P.Val.Filter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, file="diff.txt", sep="\t", quote=F, col.names=F)

diffGeneExp=data[row.names(diffSig),]
diffGeneExpOut=rbind(id=paste0(colnames(diffGeneExp),"_",Type), diffGeneExp)
write.table(diffGeneExpOut, file="diffGeneExp.txt", sep="\t", quote=F, col.names=F)

#build heatmap.png
geneNum=50
diffSig=diffSig[order(as.numeric(as.vector(diffSig$logFC))),]
diffGeneName=as.vector(rownames(diffSig))
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(2*geneNum)){
    hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
    hmGene=diffGeneName
}
hmExp=data[hmGene,]
Type=c(rep("Con",conNum),rep("DIRI",DILINum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
png(file="heatmap.png", width=6400, height=6400, res=600)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 16,
         fontsize_row=16,
         fontsize_col=16)
dev.off()
