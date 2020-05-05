
dataset_path<-"~/ParkinsonBayes/dataset_real2_elastic_loocv_logbyT/dataset_real2_loocv.rds"
result_path<-"~/ParkinsonBayes/dataset_real2_elastic_loocv_logbyT/"
#####################################################################
setwd(result_path)
one_dataset<-readRDS(dataset_path)
dataset <- one_dataset$dataset
otu_names <- one_dataset$otu_names
genus_names <- one_dataset$genus_names

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
# X_s<- sweep(X_c,2, apply(X_c,2,sd), "/")

#lambda <- scan("/home/xiw378/canola/lassoglmm/realdata/ethnicity/lassoout/lambda.txt")
lambda<-readRDS("./_lambda.rds")

nlambda<-length(lambda)
############################
best_lam.n<-which.min(test_error)
beta<-matrix(NA,ncol=382,nrow=327)
for(i in 1:n){
  beta[i,]<-readRDS(paste0("./fold_",i,"/lambda_",best_lam.n,"beta.rds"))
}
colnames(beta)<-names(readRDS(paste0("./fold_",i,"/lambda_",best_lam.n,"beta.rds")))
median_absolute<-apply(beta, 2, function(x) median(abs(x)))
medians<-apply(beta, 2, function(x) median(x))
non_zero_n<-sum(medians!=0)
non_zero_index<-which(medians!=0)
ab_me_order.n<-order(median_absolute,decreasing=TRUE)
ab_me_order<-colnames(beta)[ab_me_order.n]
sign_of_top_beta<-sign(medians)[ab_me_order.n[1:non_zero_n]]
selected_genus<-genus_names[ab_me_order.n[1:non_zero_n]]
top_otus<-vector()
A<-vector()
for(i in 1:length(selected_genus)){a<-unlist(strsplit(selected_genus[i],"f__"))
A[i]<-a[length(a)]
if(is.na(unlist(strsplit(A[i],".g__"))[2])==TRUE){B<-c(unlist(strsplit(A[i],".g__")),"unclassified")}else{B<-unlist(strsplit(A[i],".g__"))}
top_otus[i]<-paste0(B[1],'/',B[2])
}
#####final steps by hands to unify the names, then: top_names<-top_otus[1:25]
saveRDS(top_names,file = "./_top_genera_names.rds")
