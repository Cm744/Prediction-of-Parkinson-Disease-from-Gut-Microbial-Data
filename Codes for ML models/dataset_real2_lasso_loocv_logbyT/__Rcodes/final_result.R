
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
library(pROC)
#library (HTLR)
library (HTLR, lib.loc = "/home/longhai/Rdev/HTLR_3.1-1")
#library (htlr, lib.loc = "/home/longhai/Rdev/htlr_3.1-0")
#library (HTLR, lib.loc = "/home/longhai/Rdev/HTLR_3.1-0")
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
lambda <- scan("/home/xiw378/canola/lassoglmm/realdata/ethnicity/lassoout/lambda.txt")



test_error<-rep(0,times=200)

for(i in 1:200){
  test_error[i]<-readRDS(paste0("./_test_error/lambda_",i,".rds"))
  

}


#nlambda <- length (lambda)
#n <- one_dataset$n
#pred_matrix <- array (0, dim = c(n,2,nlambda))

#for (i in 1:n) {
 # for (j in 1:nlambda)
   # pred_matrix[i,,j] <- readRDS(paste0(result_path,"fold_",i,"/lambda_",j,"test_prob.rds"))}



#saveRDS(pred_matrix, file = paste0(result_path, "all_pred_matrix.rds"))
pred_matrix<-readRDS(paste0(result_path, "allpred_matrix.rds"))
allpred_matrix <- readRDS("~/ParkinsonBayes/dataset_real2_lasso_loocv_logbyT/allpred_matrix.rds")

nv_lasso <- matrix(0,n,nlambda)

for (i in 1:n) {
  for (j in 1:nlambda){
    
    nv_lasso[i,j] <- readRDS(paste0(result_path,"fold_",i,"/lambda_",j,"num_otu_retained.rds"))
  }
}

n <- one_dataset$n
pred_matrix_cov <- array (0, dim = c(n,2))

for (i in 1:n) {
    pred_matrix_cov[i,] <- exp(readRDS(paste0(result_path,'fold_',i,'/log_prior_test','.rds')))
}


## find predictive metrics against lambda
y_true <- as.numeric(one_dataset$dataset$parkinson)

auc_lasso <- rep (0, nlambda)
amlp_lasso <- rep(0, nlambda)
er_lasso <- rep(0, nlambda)

for (i in 1:nlambda) {
  eval <- evaluate_pred (pred_matrix[,,i],y_true)
  auc_lasso[i]<-roc(y_true,pred_matrix[,2,i])$auc
  #auc_lasso[i] <- eval$auc
  er_lasso[i] <- eval$er
  amlp_lasso [i] <- eval$amlp
}

which.max(auc_lasso)
which.min(er_lasso)
which.min(amlp_lasso)
#pdf(file ='./final_plot.pdf')

par(mfrow=c(1,3))
plot(log(lambda), auc_lasso)
plot(log(lambda), er_lasso)
plot(log(lambda), amlp_lasso)

PDorder <- order (y_true)
par(mfrow=c(2,2))
eval_lasso_PD_cov <- evaluate_pred(pred_matrix_cov[PDorder,], y_true[PDorder], showplot = TRUE)
eval_lasso_PD <-evaluate_pred(pred_matrix[PDorder,,21], y_true[PDorder], showplot = TRUE)
eval_lasso_PD <-evaluate_pred(pred_matrix[PDorder,,22], y_true[PDorder], showplot = TRUE)
eval_lasso_PD <-evaluate_pred(pred_matrix[PDorder,,34], y_true[PDorder], showplot = TRUE)

par(mfrow = c(1,2), mar=c(4,4,2,1), las = 1)
qtranf <- qunif
plot(qtranf(pred_matrix[PDorder,2,21]), col = y_true[PDorder], 
     pch = c(1,3)[y_true[PDorder]], yaxt='n', 
     ylab = "Predictive Probabilities",
     main = title (main = sprintf ("Error Rate = %4.2f%%, AUC=%5.3f, AMLP = %5.3f", 
                                   er_lasso[21]*100, auc_lasso[21], amlp_lasso[21]), 
                   cex = 0.8, line = 0.5))
prob_labels <- seq (0,1,by=0.1)
axis(side = 2, at = qtranf(prob_labels), labels = prob_labels)
abline(h= qtranf(seq(0,1, by = 0.1)), col="grey", lty=2)

plot(roc(y_true, pred_matrix_cov[,2]), ylim = c(0,1), lty=2)
lines(roc(y_true, pred_matrix[,2,21]), lty=1)


## some code for looking at the predictive results
par(mfrow = c(1,2), mar=c(4,4,2,1), las = 1)

plot(one_dataset$dataset$logTR, qnorm(pred_matrix[,2,21]),  col = y_true)
abline (v = log (10000))

plot(one_dataset$dataset$logTR, qnorm(eval_lasso_PD$probs_at_truelabels),  col = y_true)
abline (v = log (10000))

## looking at the correlation of wrong prediction with total reads

plot(one_dataset$dataset$logTR, qnorm(pred_matrix[,2,21]),  col = y_true)
abline (v = log (10000))

plot(one_dataset$dataset$age, qnorm(pred_matrix[,2,21]),  col = y_true)
plot(jitter(as.numeric(one_dataset$dataset$sex)), jitter(one_dataset$dataset$age), 
     col= heat.colors(100) [cut(pred_matrix[,2,21], breaks=seq(0,1,by=0.01))] 
)

plot(one_dataset$dataset$logTR, jitter(one_dataset$dataset$age), 
     col= heat.colors(100) [cut(eval_lasso_PD$probs_at_truelabels, breaks=seq(0,1,by=0.01))] 
)

par(mfrow = c(1, 1))
plot(log(lambda),test_error, type = "l", lty = 1, col = 1, lwd = 4, xlab = "log(lambda)", ylab =  "  test  Error rate")
abline( h = 0.3455657, col = "2", lty = "dotted", lwd = 5,legend("topright", "age+sex_no_otu 0.3455657", col = "2",
                                                                                  text.col = "black", lty="dotted",lwd=3,box.lty=0))
#dev.off()
which.min(test_error)
#pre<-allpred_matrix[,,22]
yreal<-as.numeric(dataset[,"parkinson"])
#for (i in 1:200){

Auc<-vector()
precision<-vector()
recall<-vector()
F_score<-vector()
for (i in 1:200){
  Auc[i]<-roc(yreal,allpred_matrix[,2,i])$auc
  pred<-data.frame(prob=allpred_matrix[,2,i],obs=yreal)
  tp <- sum(pred$prob > 0.5 & pred$obs == 2)
  fp <- sum(pred$prob > 0.5 & pred$obs == 1)
  tn <- sum(pred$prob < 0.5 & pred$obs == 1)
  fn <- sum(pred$prob < 0.5 & pred$obs == 2)
  precision[i]<-tp/(tp+fp)
  recall[i]<-tp/(tp+fn)
  F_score[i]<-2/(1/precision[i]+1/recall[i])
 
}
par(mfrow=c(2,2))
plot(Auc, ylab = "AUC", col="1")
plot(recall, ylab = "recall", col="2")
plot(precision, ylab = "precision", col="3")
plot(F_score, ylab = "F_score", col="4")
#dev.off()