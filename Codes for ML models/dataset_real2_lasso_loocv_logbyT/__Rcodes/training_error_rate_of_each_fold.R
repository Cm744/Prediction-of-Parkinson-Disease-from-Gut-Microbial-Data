

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
##############################
library(glmnet)

#fitlasso_loocv <- function(X_cov, datafile, notu) {

nsample <- n
#nsample <- 3  
Y_name <- "parkinson"
X_name <- otu_names
X_cov<-c("age","sex")

fit_formula <- as.formula(paste0(Y_name, "~", paste(c(X_cov, X_name), collapse = "+")))
X <- model.matrix(fit_formula,data=dataset)[,-1]

Y <- dataset$parkinson
lambda <- scan("/home/xiw378/canola/lassoglmm/realdata/ethnicity/lassoout/lambda.txt")

if (!exists("irep")) irep <- 1
i <- irep #folds=327



##200 represents 200 lamdas 
#$result <- matrix(NA,200,6))
#for (i in 1:nsample){
X_tr <- X[-i, ]
Y_tr <- as.numeric(Y[-i])
X_ts <- X[i, ]
Y_ts <- as.numeric(Y[i])

train_error_rate<-vector()

for(j in 1:200){

train_pedict_points<-readRDS(paste0('./fold_',i,'/lambda_',j,'train_point.rds'))
train_error_rate[j]<-1-mean(train_pedict_points==Y_tr)
}

#pdf(file =paste0('./fold_',i,'/_training_error_rate_of_current_fold.pdf'))
#par(mfrow = c(1, 1))
plot(log(lambda),train_error_rate, type = "l", lty = 1, col = 1, lwd = 4, xlab = "log(lambda)", ylab =  "  train  Error rate of current fold")
#dev.off()