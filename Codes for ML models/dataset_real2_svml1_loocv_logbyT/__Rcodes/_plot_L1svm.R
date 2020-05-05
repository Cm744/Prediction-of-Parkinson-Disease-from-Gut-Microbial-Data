


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
y_true<-as.numeric(Y)
lambda <- readRDS("/home/mac744/ParkinsonBayes/dataset_real2_svml1_loocv_logbyT/_lambda.rds")
nlambda<-length(lambda)

##################
non_otu_predmatrix<-readRDS(paste0(result_path,"non_otu_predmatrix.rds"))

non_otu_er<-mean(y_true!=max.col(non_otu_predmatrix))


SVM_test_error<-rep(0,times=nlambda)
for(i in 1:nlambda){
  SVM_test_error[i]<-readRDS(paste0("./_test_error/SVM_lambda_",i,".rds"))

  
}
plot(SVM_test_error)
min(SVM_test_error)

####################
####################
####################

SVM_allpred_point<-matrix(NA,327,ncol=nlambda)

for (l in 1:nlambda) {
  SVM_allpred_point[,l]<-readRDS(paste0("./_test_error/SVM_lambda_",l,"_point.rds"))

  
}
#pred<- readRDS("~/ParkinsonBayes/dataset_real2_lasso_loocv_logbyT/allpred_matrix.rds")
SVM_Auc<-vector()
SVM_precision<-vector()
SVM_recall<-vector()
SVM_F_score<-vector()
for (i in 1:nlambda){
  SVM_allpred_point[,i][which(SVM_allpred_point[,i]=="Yes")]<-2
  SVM_allpred_point[,i][which(SVM_allpred_point[,i]=="No")]<-1
  SVM_Auc[i]<-roc(y_true,as.numeric(SVM_allpred_point[,i]))$auc
  pred<-data.frame(pre=as.numeric(SVM_allpred_point[,i]),obs=y_true)
  tp <- sum(pred$pre ==2 & pred$obs == 2)
  fp <- sum(pred$pre ==2 & pred$obs == 1)
  tn <- sum(pred$pre==1 & pred$obs == 1)
  fn <- sum(pred$pre==1 & pred$obs == 2)
  SVM_precision[i]<-tp/(tp+fp)
  SVM_recall[i]<-tp/(tp+fn)
 SVM_F_score[i]<-2/(1/precision[i]+1/recall[i])
  
}


#pdf(file ='./_final.pdf')
#nums_comb<-c(1:otu.num)
########result of rank of amlp and error rate
par(mfrow = c(2, 2))
plot(lambda,SVM_test_error, type = "l", lty = 1, col = 1, lwd = 4, xlab = "lambda1", ylab =  "  test  Error rate SVM_SVM")

#dev.off()


####################

pdf(file=paste0("./__final_plots/L1SVM/recall_precision_f.pdf"),width=6,height = 6)
par(mar=c(4,4,2,2))
par(mfrow = c(1, 3))

plot(SVM_recall,type="l",col=4,ylab="recall L1SVM")
plot(SVM_precision,type="l",col=5,ylab="precision L1SVM")

plot(SVM_F_score,type="l",col=6,ylab="F1 L1SVM")

dev.off()


pdf(file=paste0("./__final_plots/L1SVM/L1SVM_Error_Rate.pdf"),width=6,height = 6)
par(mar=c(3,4,1,4))

plot(log(lambda), SVM_test_error,col="blue",type="l",lwd=2,ylim=c(0.29,0.50),xlab="", ylab="")
abline(h=non_otu_er,col="darkgreen",lty=2,lwd=2)
abline(h=130/327,col="red",lty=3,lwd=2)
axis(side = 4, at=seq(0,197/327,by=130/327/10), labels=seq(1,-0.5,by=-0.1))
mtext(side=4, line=2, expression(R^2), cex=1.5,las=2)
mtext(side=2, line=2, "ER_L1SVM", cex=1,las=3)
mtext(side=1, line=2, expression(log(lambda)), cex=1)
abline(h=seq(0.29,0.5,by=0.005),v=seq(-12,-4,by=1),col="lightgrey",lty=3)
legend('bottomleft', c('OTU+age+sex', 'age+sex','baseline'),col  =c('blue', 'darkgreen','red'),cex=0.8,lty=c(1,2,3),lwd=2,bg="white")
##

dev.off()
