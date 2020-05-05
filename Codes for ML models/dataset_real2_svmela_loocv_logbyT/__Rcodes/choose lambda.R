

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
#X_s<- sweep(X_c,2, apply(X_c,2,sd), "/")
#lambda <- scan("/home/xiw378/canola/lassoglmm/realdata/ethnicity/lassoout/lambda.txt")
# lambda<-c (seq(0.01 ,0.14, 0.03), seq(0.15,0.95, 0.1),seq(1,4.5,0.5))
#lambda<-seq(0.1,1,0.01)
cv.fit <- cv.sparseSVM(X, 2*as.numeric(Y)-3, nfolds = 10, ncores = 2, seed = 1234,alpha=0.5,gamma=1, preprocess = "standardize")

lambda<-cv.fit$lambda
saveRDS(lambda,file="./_lambda.rds")
