library(nnet)

dataset_path<-"~/ParkinsonBayes/dataset_real2_elastic_loocv_logbyT/dataset_real2_loocv.rds"
result_path<-"~/ParkinsonBayes/dataset_real2_elastic_loocv_logbyT/"
#####################################################################
setwd(result_path)
one_dataset<-readRDS(dataset_path)
dataset <- one_dataset$dataset
otu_names <- one_dataset$otu_names
otu.num<-one_dataset$otu.num
var_interested_name <- "parkinson"
nl <- nlevels(dataset[1, var_interested_name])
n<-one_dataset$n
folds<-one_dataset$folds
folds.index<-one_dataset$folds.index

if (!exists("irep")) irep <- 1
i <- irep #qusb folds=327 times

traindata<-dataset[-folds.index[[i]],]
testdata<-dataset[folds.index[[i]],]
logistic<-glm(parkinson~age+sex,data<-traindata,family = binomial(logit))
testyes<-predict(logistic,testdata,type='response')
trainyes<-predict(logistic,traindata,type='response')

Pr_xgiven_En_log_test<-unname(cbind(log(1-testyes),log(testyes)))
Pr_xgiven_En_log_train<-unname(cbind(log(1-trainyes),log(trainyes)))
case_predicition_non_otu<-which.max(Pr_xgiven_En_log_test)
saveRDS(Pr_xgiven_En_log_test,file = paste0('./fold_',i,'/log_prior_test','.rds'))
saveRDS(Pr_xgiven_En_log_train,file = paste0('./fold_',i,'/log_prior_train','.rds'))
saveRDS(case_predicition_non_otu,file = paste0('./fold_',i,'/case_predicition_non_otu','.rds'))
