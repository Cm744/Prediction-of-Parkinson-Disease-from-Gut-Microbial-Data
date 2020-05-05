
dataset_path<-"~/ParkinsonBayes/dataset_real2_lasso_loocv_logbyT/dataset_real2_loocv.rds"
result_path<-"~/ParkinsonBayes/dataset_real2_lasso_loocv_logbyT/"
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
TR<-exp(dataset$logTR)
##############################
library(reshape2)
library(ggplot2)

#dataset_logp <- apply(dataset[,1:382],1,function(x)(log(x+1)-log(sum(x)+length(x))) )
lambda <- scan("/home/xiw378/canola/lassoglmm/realdata/ethnicity/lassoout/lambda.txt")
nlambda<-length(lambda)
genus_names<-one_dataset$genus_names
nsample <- n
########################
trans<-apply(dataset[,1:otu.num],1,function(x) log(x+1)-log(sum(x)+length(x) ))
trans1<-t(trans)
dataset[,1:otu.num]<-trans1
Y_name <- "parkinson"
X_name <- otu_names
X_cov<-c("age","sex")
fit_formula <- as.formula(paste0(Y_name, "~", paste(c(X_cov, X_name), collapse = "+")))
X <- model.matrix(fit_formula,data=dataset)[,-1]


#selected_samples<-Y[-i][sample_order[c(1:20)]]
#group <- sort(selected_samples)


test_error<-vector()

for(m in 1:nlambda){
  test_error[m]<-readRDS(paste0("./_test_error/lambda_",m,".rds"))
  
  
}

min_er_lambda<-which.min(test_error)


if (!exists("irep")) irep <- 1
i <- irep #folds=327



beta<-readRDS(paste0("./fold_",i,"/lambda_",min_er_lambda,"beta.rds"))
beta_rank<-order(-abs(beta))

write.table(beta_rank,file = paste0("~/ParkinsonBayes/dataset_real2_loocv/_beta_rank/fold_",i,".txt"))
