


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
trans<-apply(dataset[,1:otu.num],1,function(x)log(x+1)-log(sum(x)+length(x) ))
trans1<-t(trans)
dataset[,1:otu.num]<-trans1
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
allpred_matrix <- readRDS("~/ParkinsonBayes/dataset_real2_lasso_loocv_logbyT/allpred_matrix.rds")
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
pdf(file ='./_final.pdf')
#nums_comb<-c(1:otu.num)
########result of rank of amlp and error rate
par(mfrow = c(1, 1))
plot(log(lambda),test_error, type = "l", lty = 1, col = 1, lwd = 4, xlab = "log(lambda)", ylab =  "  test  Error rate")
abline( h =  0.3455657, col = "2", lty = "dotted", lwd = 5,legend("topright", "age+sex_no_otu 0.3455657", col = "2",
                                                                  text.col = "black", lty="dotted",lwd=3,box.lty=0))
dev.off()



par(mfrow = c(1, 2))
plot(log(lambda),test_error, type = "l", lty = 1, col = 1, lwd = 4, xlab = "", ylab =  "  LR_L1 test error")
mtext(side=1, line=2, expression(log(lambda)), cex=1,las=1)

plot(log(lambda),Auc, ylab = "LR_L1 AUC",type = "l",lwd = 4,xlab = "", col="2")
mtext(side=1, line=2, expression(log(lambda)), cex=1,las=1)
par(mfrow = c(1, 3))
plot(log(lambda),recall, ylab = "LR_L1 recall",xlab ="",type = "l",lwd = 4, col="3")
mtext(side=1, line=2, expression(log(lambda)), cex=0.7,las=1)

plot(log(lambda),precision, ylab = "LR_L1 precision",xlab = "",type = "l",lwd = 4, col="4")
mtext(side=1, line=2, expression(log(lambda)), cex=0.7,las=1)

plot(log(lambda),F_score, ylab = "LR_L1 F_score",xlab = "", type = "l",lwd = 4,col="5")
mtext(side=1, line=2, expression(log(lambda)), cex=0.7,las=1)





