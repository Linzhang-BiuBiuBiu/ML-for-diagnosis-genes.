# LASSO for diagnosis genes
library(tidyverse)
library(caret)
library(glmnet)
set.seed(123)

inputFile="diffGeneExp.txt"      
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))      #???ù???Ŀ¼

rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
rt=t(rt)

x=as.matrix(rt)
y=gsub("(.*)\\_(.*)", "\\2", 
       row.names(rt))
fit=glmnet(x, y, family = "binomial", 
           alpha=1)
cvfit=cv.glmnet(x, y, family="binomial", 
                alpha=1,
                type.measure='deviance',
                nfolds = 10)

modelData <- as.data.frame(cbind(group=y,x))
x <- model.matrix(group~.,modelData)
y <- as.factor(modelData[,1])

pdf(file="cvfit.pdf",width=6,height=5.5)
plot(cvfit)
dev.off()

coef=coef(fit, s = cvfit$lambda.1se)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=lassoGene[-1]
write.table(lassoGene, file="LASSO.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)
write.table(as.matrix(coef), file="LASSO_coef.gene.txt", sep="\t", quote=F, row.names=T, col.names=T)

#SVM for filter diagnosis genes
library(e1071)
library(kernlab)
library(caret)

set.seed(123)
inputFile="diffGeneExp.txt"        
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))      #???ù???Ŀ¼

data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))

#SVM-RFE
Profile=rfe(x=data,
            y=as.numeric(as.factor(group)),
            sizes = c(2,4,6,8, seq(10,40,by=3)),
            rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
            methods="svmRadial")


pdf(file="SVM-RFE.pdf", width=6, height=5.5)
par(las=1)
x = Profile$results$Variables
y = Profile$results$RMSE
plot(x, y, xlab="Variables", ylab="RMSE (Cross-Validation)", col="darkgreen")
lines(x, y, col="darkgreen")
wmin=which.min(y)
wmin.x=x[wmin]
wmin.y=y[wmin]
points(wmin.x, wmin.y, col="blue", pch=16)
text(wmin.x, wmin.y, paste0('N=',wmin.x), pos=2, col=2)
dev.off()

featureGenes=Profile$optVariables
write.table(file="SVM-RFE.gene.txt", featureGenes, sep="\t", quote=F, row.names=F, col.names=F)
write.table(file="SVM-RFE-coef.gene.txt", Profile$variables, sep="\t", quote=F, row.names=T, col.names=T)
