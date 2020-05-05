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
#lambda <- scan("/home/xiw378/canola/lassoglmm/realdata/ethnicity/lassoout/lambda.txt")
set.seed(1234)
er<-vector()
alpha<-seq(0,1,0.005)
if (!exists("irep")) irep <- 1
i <- irep #folds=201



c<-cv.glmnet(x=X, y=Y,nfolds =327, family = "binomial",alpha=alpha[i])
g<-glmnet(x=X,y=as.numeric(Y),lambda=c$lambda.min,alpha=alpha[i],family = 'binomial')
pp<-predict(g,newx=X,type="response")
p<-cbind(1-pp,pp)
po<-max.col(p)
er<-1-mean(po==as.numeric(Y))
l<-c$lambda


saveRDS(l,paste0( "./_choose/_lambda_of_alpha_",i,".rds"))
saveRDS(er,paste0( "./_choose/_er_of_alpha_",i,".rds"))
