


dataset_path<-"~/ParkinsonBayes/dataset_real2_svmela_loocv_logbyT/dataset_real2_loocv.rds"
result_path<-"~/ParkinsonBayes/dataset_real2_svmela_loocv_logbyT/"
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
library(sparseSVM)
library(nnet)
trans<-apply(dataset[,1:382],1,function(x) log(x+1)-log(sum(x)+length(x) ))
trans1<-t(trans)
dataset[,1:382]<-trans1

#tranfOTU <- function(x) log(x+1)-log(sum(x)+length(x))
#for (i in 1:nrow(dataset))
#{
 # dataset[i,1:382] <- tranfOTU(dataset0[i,1:382])
#}


#fitlasso_loocv <- function(X_cov, datafile, notu) {
   
  nsample <- n
   
  Y_name <- "parkinson"
  X_name <- otu_names
  X_cov<-c("age","sex")
  
 fit_formula <- as.formula(paste0(Y_name, "~", paste(c(X_cov, X_name), collapse = "+")))
  X <- model.matrix(fit_formula,data=dataset)[,-1]
  
  Y <- dataset$parkinson
  #X_c <- sweep(X, 2, colMeans(X), "-")
 # X_s<- sweep(X_c,2, apply(X_c,2,sd), "/")
  lambda <- readRDS("/home/mac744/ParkinsonBayes/dataset_real2_svmela_loocv_logbyT/_lambda.rds")
  #lambda <- scan("/home/xiw378/canola/lassoglmm/realdata/ethnicity/lassoout/lambda.txt")
  
  nlambda<-length(lambda)
  if (!exists("irep")) irep <- 1
  i <- irep #folds=327
  
  #coef<-rep(0,times=384)
  
  
  #Y_svm<-2*as.numeric(Y)-3
    X_tr <- X[-i, ]
    Y_tr <- Y[-i]
    X_ts <- X[i, ]
    Y_ts <-Y[i]
    #Y_tr_label<-Y[-i]
    
    for(j in 1:nlambda){
      ELA<-sparseSVM(X_tr, Y_tr, alpha = 0.5, gamma =1,lambda=lambda[j], preprocess = "standardize")
     ELA_tr_pred<-as.vector(predict(ELA, X_tr))
    ELA_tr_er<-1-mean(ELA_tr_pred==as.character(Y_tr))
    ELA_ts_pred<-as.vector(predict(ELA, X_ts))
    ELA_intercept<-ELA$weights[1]
    ELA_coef<-ELA$weights[-1]
    names(ELA_coef)<-colnames(X_tr)
    ELA_otu_retained<-ELA_coef[-c(1:2)]
    ELA_num_otu_retained <- length(which(ELA_otu_retained!=0))
    ELA_index_otu_retained<-which(ELA_coef[-c(1:2)]!=0)
    ELA_intercept<-ELA$weights[1]
    tr_relative_distance<-X_tr%*%ELA_coef+ELA_intercept
    ts_relative_distance<-X_ts%*%ELA_coef+ELA_intercept
    ##########logitics regression input:relative distance, output:parkinson
    SL_traindata<-data.frame(relative_distance=tr_relative_distance,parkinson=Y_tr)
    SL_testdata<-data.frame(relative_distance=ts_relative_distance,parkinson=Y_ts)
    S_logistic<-glm(parkinson~relative_distance,data=SL_traindata,family = binomial(logit))
    
    trainyes<-predict(S_logistic,SL_traindata,type = "response")
    Pr_of_PD_train<-unname(cbind(1-trainyes,trainyes))
    train_pred_point<-max.col(Pr_of_PD_train)
    S_L_tr_er<-mean(train_pred_point!=as.numeric(Y_tr))
    testyes<-predict(S_logistic,SL_testdata,type = "response")
    Pr_of_PD_test<-unname(cbind(1-testyes,testyes))
    
    saveRDS(ELA_tr_pred,file = paste0("./fold_",i,"/SVM//lambda_",j,"train_point.rds"))
    saveRDS(ELA_ts_pred,file =paste0 ("./fold_",i,"/SVM/lambda_",j,"test_point.rds"))
    saveRDS(ELA_num_otu_retained,file = paste0("./fold_",i,"/SVM/lambda_",j,"num_otu_retained.rds"))
    saveRDS(ELA_index_otu_retained,file = paste0("./fold_",i,"/SVM/lambda_",j,"index_otu_retained.rds"))
    saveRDS(ELA_tr_er,file = paste0("./fold_",i,"/SVM/lambda_",j,"train_error_rate.rds"))
    saveRDS(ELA_otu_retained,file = paste0("./fold_",i,"/SVM/lambda_",j,"beta.rds"))
    saveRDS(ELA_coef,file = paste0("./fold_",i,"/SVM/lambda_",j,"coef.rds"))
    
    saveRDS(ELA_intercept,file = paste0("./fold_",i,"/SVM/lambda_",j,"intercept.rds"))
    saveRDS(tr_relative_distance,file = paste0("./fold_",i,"/SVM/lambda_",j,"tr_relative_distance.rds"))
    saveRDS(ts_relative_distance,file = paste0("./fold_",i,"/SVM/lambda_",j,"ts_relative_distance.rds"))
   ###########logistic result:
     saveRDS(Pr_of_PD_train,file = paste0("./fold_",i,"/S_L/lambda_",j,"Pr_of_PD_train.rds"))
    saveRDS(Pr_of_PD_test,file = paste0("./fold_",i,"/S_L/lambda_",j,"Pr_of_PD_test.rds"))
    saveRDS(S_L_tr_er,file = paste0("./fold_",i,"/S_L/lambda_",j,"train_error.rds"))
    
    
 }
