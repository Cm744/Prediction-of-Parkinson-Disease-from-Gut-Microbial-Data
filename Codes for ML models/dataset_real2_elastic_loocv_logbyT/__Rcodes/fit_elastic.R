
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
##############################
library(glmnet)
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
 # lambda <- scan("/home/xiw378/canola/lassoglmm/realdata/ethnicity/lassoout/lambda.txt")
  lambda<-readRDS("./_lambda.rds")
 alpha_fit<-readRDS("./best_alpha.rds")
  if (!exists("irep")) irep <- 1
  i <- irep #folds=327
  
 
 
    X_tr <- X[-i, ]
    Y_tr <- as.numeric(Y[-i])
    X_ts <- X[i, ]
    Y_ts <- as.numeric(Y[i])
    
    for(j in 1:length(lambda)){
      
      lasso.fit <- glmnet(x=X_tr, y=Y_tr, lambda=lambda[j], family = "binomial",alpha=alpha_fit)
    
    ## calculate the training result
    tr_pred_yes <- predict(lasso.fit, newx = X_tr, type="response")
    tr_pred_no <- 1 - tr_pred_yes
    tr_pred_matrix <- cbind(tr_pred_no , tr_pred_yes)
    tr_prediction <- max.col( tr_pred_matrix)
    tr_er <- 1-mean(tr_prediction==Y_tr)
    
    ## calculate the validation result
    test_pred_yes <- t(predict(lasso.fit, newx = t(X_ts), type = "response"))
    test_pred_no <- 1 - test_pred_yes
    ts_pred_matrix <- cbind(test_pred_no, test_pred_yes)
    ts_prediction <- max.col(ts_pred_matrix)
    num_otu_retained <- length(which(lasso.fit$beta[-c(1:2),]!=0))
    index_otu_retained<-which(lasso.fit$beta[-c(1:2),]!=0)
    beta<-lasso.fit$beta[-c(1:2),]
    #truelab=rep(Y_ts,200)
  
    saveRDS(tr_pred_matrix,file = paste0("./fold_",i,"/lambda_",j,"train_prob.rds"))
    saveRDS(ts_pred_matrix,file = paste0("./fold_",i,"/lambda_",j,"test_prob.rds"))
    saveRDS(tr_prediction,file = paste0("./fold_",i,"/lambda_",j,"train_point.rds"))
    saveRDS(ts_prediction,file =paste0 ("./fold_",i,"/lambda_",j,"test_point.rds"))
    saveRDS(num_otu_retained,file = paste0("./fold_",i,"/lambda_",j,"num_otu_retained.rds"))
    saveRDS(index_otu_retained,file = paste0("./fold_",i,"/lambda_",j,"index_otu_retained.rds"))
    saveRDS(tr_er,file = paste0("./fold_",i,"/lambda_",j,"train_error_rate.rds"))
    saveRDS(beta,file = paste0("./fold_",i,"/lambda_",j,"beta.rds"))
    }
