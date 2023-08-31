library(tidyverse)
library(readxl)
library(openxlsx)
library(glmnet)
library(e1071)
library(stringr)
library(randomForest)
library(RSNNS)
library(neuralnet)
library(caret)
library(h2o)
library(NeuralNetTools)
library(ggplot2)
library(rpart)
library(rpart.plot)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(digits = 15)

train_dif <- "diffGeneExp.txt"
rt <- read.table(train_dif, header=T, sep="\t", check.names=T)
rownames(rt) <- rt[,1]
data <- rt[,-1]
y <- gsub("(.*)\\_(.*)","\\2",colnames(data))

group <- as.factor(y)
model_data <- t(data)
model_data_group <- t(rbind(group=group,data))


#RF
set.seed(1234)
rfcla <- randomForest(x=model_data, 
                      y=group, 
                      ntree=500,
                      proximity = TRUE)

rftune <- tuneRF(x=model_data,
                 y=group,stepFactor = 1,
                 ntreeTry = 500)
rftune_min_OBB <- min(which(grepl(min(rftune[,2]),rftune[,2])))
rftune_optimal <- rftune[rftune_min_OBB,1]

#optimal tree number for RF
RandomForest <- randomForest(x=model_data, y=group,
                          ntree=500,mtry=rftune_optimal)

pdf(file="RandomForest.pdf",width = 6,height = 5.5)
plot(RandomForest)
varImpPlot(RandomForest,sort=T,
           main = "The importance of significent genes",
           color = c("blue"))
dev.off()

#NN
y <- gsub("(.*)\\_(.*)", "\\2", colnames(data))
group <- as.factor(y)
model_data <- t(data)
model_data_group <- t(rbind(group=group,data))

#formula for NN
model_data_group[,1] <- factor(model_data_group[,1])
n<- colnames(model_data_group)
form <- as.formula(paste("group~",paste(n[!n %in% "use"],collapse = "+")))
fit <- neuralnet(form,data=model_data_group, hidden = c(1,2,3),
                 err.fct = "sse",
                 linear.output = T)

pdf(file="Neuralnet.pdf",width=7,height = 6)
par(cex = 1)
plotnet(fit,max_sp=1,pad_x=0.7,prune_lty=0.5,circle_cex=3,
        pos_col = "red", neg_col = "grey")
dev.off()

#GBM
Sys.setenv(JAVA_HOME="C:\\Program Files\\Java\\jdk1.8.0_341")
h2o.init(nthreads = 8,max_mem_size = "8G")

#as.h2o formyla
gbm_model_data <- as.h2o(model_data_group)

target <- "group"
predictors <- colnames(model_data_group)[2:ncol(model_data_group)]

#Model for GBM
gbm_grid_cla <- h2o.gbm(x=predictors,
                        y=target,
                        distribution="AUTO",
                        training_frame = gbm_model_data,
                        ntree=100,
                        learn_rate=0.01,
                        sample_rate = 0.8,
                        col_sample_rate = 0.6,
                        seed=1234)
aml <- h2o.automl(y = target,
                  training_frame = gbm_model_data,
                  max_models = 6,
                  seed = 1)
pdf(file="GBM_Variable_importance.PDF",width = 6,height = 5.5)
h2o.varimp_heatmap(aml)
h2o.varimp_plot(gbm_grid_cla,ncol(model_data_group)/8)
dev.off()

#DT
model_data_group <- as.data.frame(model_data_group)
rpart_model <- rpart(group~.,
                     data=model_data_group,
                     method = "class",
                     cp=0.000001)

#figure for DT
pdf(file="Decision trees-1.pdf",width = 4,height = 3)
rpart.plot(rpart_model,type=2,extra = "auto",
           under=T,fallen.leaves=F,cex=0.7,main="Decision trees")
plotcp(rpart_model)
dev.off()
#optimaze the tree number
bestcp <- rpart_model$cptable[which.min(rpart_model$cptable[,"xerror"]),"CP"]
rpart_model.pruned <- prune(rpart_model,cp=bestcp)
par(family="STKaiti")
pdf(file="Decision trees.pdf",width = 10,height = 10)
rpart.plot(rpart_model.pruned,type=2,extra = "auto",
          under=T,fallen.leaves=F,cex=2,main="Optimal Decision trees")
dev.off()

rf_sig <- RandomForest$importance
RandomForest_result <- rf_sig[order(rf_sig[,1],decreasing = T),]
RF_number <- rownames(as.matrix(RandomForest_result))[c(1:50)]
NeuralNet_result <- fit$result.matrix
GBM_result <- summary(gbm_grid_cla)
GBM_number <- GBM_result[,1][c(1:50)]
rpart_number <- as.matrix(rpart_model.pruned$variable.importance)
result_variable <- cbind(RF_number=RF_number,GBM_number=GBM_number,
                         rt_number=rownames(rpart_number)[c(1:50)])

write.table(result_variable,file="result_variable.txt",quote=F,sep = "\t")

rf_imp<- as.data.frame(RandomForest$importance)
rf_imp <- cbind(rownames(rf_imp),rf_imp)

NN_imp <- as.data.frame(fit$result.matrix)
NN_imp <- cbind(rownames(NN_imp),NN_imp)

gbm_imp <- as.data.frame(GBM_result)
gbm_imp <- cbind(rownames(gbm_imp),gbm_imp)

rp_imp <- as.data.frame(rpart_model$variable.importance)
rp_imp <- cbind(rownames(rp_imp),rp_imp)
imp_all <- list(rf_imp,NN_imp,gbm_imp,rp_imp)
write.xlsx(imp_all,file = "imp_all.xlsx")
