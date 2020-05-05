

dataset_path<-"~/ParkinsonBayes/dataset_real2_svml1_loocv_logbyT/dataset_real2_loocv.rds"
result_path<-"~/ParkinsonBayes/dataset_real2_svml1_loocv_logbyT/"
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
##############################
library(glmnet)
trans<-apply(dataset[,1:382],1,function(x) log(x+1)-log(sum(x)+length(x) ))
trans1<-t(trans)
dataset[,1:382]<-trans1

#fitlasso_loocv <- function(X_cov, datafile, notu) {

nsample <- n
#nsample <- 3  
Y_name <- "parkinson"
X_name <- otu_names
X_cov<-c("age","sex")

fit_formula <- as.formula(paste0(Y_name, "~", paste(c(X_cov, X_name), collapse = "+")))
X <- model.matrix(fit_formula,data=dataset)[,-1]

  Y <- dataset$parkinson
  lambda <- readRDS("/home/mac744/ParkinsonBayes/dataset_real2_svml1_loocv_logbyT/_lambda.rds")
  nlambda<-length(lambda)
if (!exists("irep")) irep <- 1
j <- irep #100





SVM_test_point<-vector()
S_L_test_point<-vector()
for (i in 1:folds){
 SVM_test_point[i]<-readRDS(paste0 ("./fold_",i,"/SVM/lambda_",j,"test_point.rds"))
 S_L_test_point[i]<-max.col(readRDS(paste0("./fold_",i,"/S_L/lambda_",j,"Pr_of_PD_test.rds")))
 
 }

SVM_test_error<-1-mean(SVM_test_point==as.character(Y))
S_L_test_error<-1-mean(S_L_test_point==as.numeric(Y))


saveRDS(SVM_test_error,file = paste0("./_test_error/SVM_lambda_",j,".rds"))
saveRDS(SVM_test_point,file = paste0("./_test_error/SVM_lambda_",j,"_point.rds"))

saveRDS(S_L_test_error,file = paste0("./_test_error/S_L_lambda_",j,".rds"))
saveRDS(S_L_test_point,file = paste0("./_test_error/S_L_lambda_",j,"_point.rds"))
