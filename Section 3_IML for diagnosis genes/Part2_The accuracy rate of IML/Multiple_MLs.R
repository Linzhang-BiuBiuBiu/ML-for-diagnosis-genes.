library(e1071)
library(kernlab)
library(caret)
library(forecast)
library(h2o)
library(tidyverse)
library(readxl)
library(glmnet)
library(e1071)
library(stringr)
library(randomForest)
library(RSNNS)
library(neuralnet)
library(openxlsx)
library(caret)
library(h2o)
library(NeuralNetTools)
library(ggplot2)
library(rpart)
library(rpart.plot)
library(Metrics)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dif_genes_train <- read.table("train-diff.txt",header=T,
                        sep="\t", 
                        check.names=T)[,1]
dif_genes_test <- read.table("test.normalize.txt",header=T,
                              sep="\t", 
                              check.names=T)[,1]

dif_genes <- as.matrix(read.table("filtered_gege.txt",header=F,
                        sep="\t", 
                        check.names=T))
dif_genes <- intersect(dif_genes,dif_genes_test)
dif_genes <- intersect(dif_genes,dif_genes_train)

rt <- read.table("train.normalize.txt",header=T,
                         sep="\t", 
                         check.names=T)
rownames(rt) <- rt[,1]
train_data <- t(rt[dif_genes,][,-1])
train_data[,1:ncol(train_data)] <- sapply(train_data[,1:ncol(train_data)], as.numeric)

rt <- read.table("test.normalize.txt",header=T,
                         sep="\t", 
                         check.names=T)

rownames(rt) <- rt[,1]
test_data <- t(rt[dif_genes,][,-1])
test_data[,1:ncol(test_data)] <- sapply(test_data[,1:ncol(test_data)], as.numeric)

x <- train_data
y <- gsub("(.*)\\_(.*)", "\\2", rownames(train_data))
#支持向量机训练集模型
set.seed(234)
model <- as.data.frame(cbind(group=as.factor(y),train_data))
model$group <- as.factor(model$group)
y_test <- gsub("(.*)\\_(.*)", "\\2", rownames(test_data))
model_test <- as.data.frame(cbind(group=as.factor(y_test),test_data))
model_test$group <- as.factor(model_test$group)
svm_train = tune.svm(group~.,
                     data=model,
                     kernel="linear",
                     cost=seq(1,200,1),
                     fold=10)
svm_train$best.model

svm_train_data <- as.character(predict(svm_train$best.model,model))
svm_train_accuracy <- accuracy(svm_train_data,model$group)

svm_test <- as.character(predict(svm_train$best.model,test_data))
svm_test_accuracy <- accuracy(svm_test,model_test$group)

#lasso
set.seed(234)
x <- train_data
y <- as.factor(ifelse(y=="con",0,1))
lasso_train <- cv.glmnet(x,y,alpha=1,
                         lambda=seq(1,2000,1),
                         nflod=10,
                         family="binomial",
                         type.measure="class")
lasso_min <- lasso_train$lambda.min
lasso_best <- glmnet(x,y,alpha=1,lambda=lasso_min,family = "binomial")

lasso_test <- predict(lasso_best,test_data)
lasso_test <- as.factor(ifelse(lasso_test>0.5,2,1))
lasso_test_accuracy <- accuracy(as.vector(lasso_test),model_test$group)

lasso_tra <- predict(lasso_best,train_data)
lasso_tra <- as.factor(ifelse(lasso_tra>0.5,1,0))
lasso_train_accuracy <- accuracy(as.vector(lasso_tra),y)

#随机森林
set.seed(234)
rftune <- tuneRF(x=model,y=y,stepFactor = 1,ntreeTry = 500,nflod=10)
rftune_min_OBB <- min(which(grepl(min(rftune[,2]),rftune[,2])))
rftune_optimal <- rftune[rftune_min_OBB,1]

#optimalize_depth
RandomForest <- randomForest(x=model, 
                             y=y,
                             ntree=500,
                             mtry=rftune_optimal,
                             nflod=10)
rf_train_data <- predict(RandomForest,model)
rf_train_data <- ifelse(rf_train_data==0,1,2)
rf_train_accuracy <- accuracy(as.vector(rf_train_data),model$group)

rf_test <- predict(RandomForest,model_test)
rf_test <- ifelse(rf_test==0,1,2)
rf_test_accuracy <- accuracy(as.vector(rf_test),model_test$group)

#神经网络

set.seed(234)
n <- colnames(model)
form <- as.formula(paste("group~",paste(n[!n %in% "use"],collapse = "+")))
model[,1:ncol(model)] <- lapply(model[,1:ncol(model)], as.numeric)

model_test[,1:ncol(model_test)] <- lapply(model_test[,1:ncol(model_test)],
                                          as.numeric)

nn_train <- neuralnet(form,data=model, hidden = c(1,2,3),
                 err.fct = "sse",
                 linear.output = T)

nn_train_data <- predict(nn_train,model)
nn_train_data <- ifelse(nn_train_data>1.5,2,1)
nn_train_accuracy <- accuracy(as.vector(nn_train_data),model$group)

nn_test <- predict(nn_train,model_test)
nn_test <- ifelse(nn_test>1.5,2,1)
nn_test_accuracy <- accuracy(as.vector(nn_test),model_test$group)


#梯度提升机
set.seed(234)
Sys.setenv(JAVA_HOME="C:\\Program Files\\Java\\jdk1.8.0_341")
h2o.init(nthreads = 8,max_mem_size = "8G")

#transform to h2o
model$group <- as.factor(model$group)
gbm_model_data <- as.h2o(model)
gbm_model_test <- as.h2o(model_test)
target <- "group"
predictors <- colnames(model)[2:ncol(model)]

#model
gbm_grid_cla <- h2o.gbm(x=predictors,
                        y=target,
                        distribution="bernoulli",
                        training_frame = gbm_model_data,
                        ntree=100,
                        learn_rate=0.01,
                        sample_rate = 0.8,
                        col_sample_rate = 0.6,
                        seed=1234,nfolds=10)

gbm_train_data <- as.data.frame(h2o.predict(gbm_grid_cla,newdata=gbm_model_data))
gbm_train_accuracy <- accuracy(as.vector(gbm_model_data$group),gbm_train_data$predict)

gbm_test <- as.data.frame(h2o.predict(gbm_grid_cla,newdata=gbm_model_test))
gbm_test_accuracy <- accuracy(as.vector(gbm_model_test$group),gbm_test$predict)

#决策树
set.seed(234)
model_data_group <- as.data.frame(model)
rpart_model <- rpart(group~.,data=model,method = "class",
                     cp=0.000001)
dt_train <- predict(rpart_model,model)
dt_train_accuracy <- accuracy(as.vector(ifelse(dt_train[,1]>0.5,1,2)),model$group)

dt_test <- predict(rpart_model,model_test)
dt_test_accuracy <- accuracy(as.vector(ifelse(dt_test[,1]>0.5,1,2)),model_test$group)

ml_accuracy <- data.frame(lasso_train_accuracy=lasso_train_accuracy,
                          lasso_test_accuracy=lasso_test_accuracy,
                          svm_train_accuracy=svm_train_accuracy,
                          svm_test_accuracy=svm_test_accuracy,
                          rf_train_accuracy=rf_train_accuracy,
                          dt_train_accuracy=dt_train_accuracy,
                          dt_test_accuracy=dt_test_accuracy,
                          gbm_train_accuracy=gbm_train_accuracy,
                          gbm_test_accuracy=gbm_test_accuracy,
                          nn_train_accuracy=nn_train_accuracy,
                          nn_test_accuracy=nn_test_accuracy)

ml_accuracy <- data.frame(train=c(lasso=lasso_train_accuracy,
                          svm=svm_train_accuracy,
                          rf=rf_train_accuracy,
                          dt=dt_train_accuracy,
                          gbm=gbm_train_accuracy,
                          nn=nn_train_accuracy),
                          
                          test=c(lasso_test=lasso_test_accuracy,
                          svm_test=svm_test_accuracy,
                          rf_test=rf_test_accuracy,
                          dt_test=dt_test_accuracy,
                          gbm_test=gbm_test_accuracy,
                          nn_test=nn_test_accuracy))

write.csv(ml_accuracy,file = "ml_accuracy.csv")

