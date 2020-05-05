

dataset_path<-"~/ParkinsonBayes/dataset_real2_svml1_loocv_logbyT/dataset_real2_loocv.rds"
result_path<-"~/ParkinsonBayes/dataset_real2_svml1_loocv_logbyT/"
#####################################################################
setwd(result_path)
one_dataset<-readRDS(dataset_path)
dataset <- one_dataset$dataset
otu_names <- one_dataset$otu_names
otu.num<-one_dataset$otu.num
genus_names<-one_dataset$genus_names
var_interested_name <- "parkinson"
nl <- nlevels(dataset[1, var_interested_name])
n<-one_dataset$n
folds<-one_dataset$folds
folds.index<-one_dataset$folds.index
##############################
library(glmnet)
library(pROC)
library(ROCR)
library(PRROC)
#library (HTLR)
library (HTLR, lib.loc = "/home/longhai/Rdev/HTLR_3.1-1")
trans<-apply(dataset[,1:382],1,function(x) log(x+1)-log(sum(x)+length(x) ))
trans1<-t(trans)
dataset[,1:382]<-trans1
#fitlasso_loocv <- function(X_cov, datafile, notu) {

nsample <- n

Y_name <- "parkinson"
X_name <- otu_names
X_cov<-c("age","sex")

fit_formula <- as.formula(paste0(Y_name, "~", paste(c(X_cov, X_name), collapse = "+")))
X <- model.matrix(fit_formula,data=dataset)[,-1]

Y <- dataset$parkinson

#lambda <- scan("/home/xiw378/canola/lassoglmm/realdata/ethnicity/lassoout/lambda.txt")
lambda <- readRDS("/home/mac744/ParkinsonBayes/dataset_real2_svml1_loocv_logbyT/_lambda.rds")

nlambda<-length(lambda)

X_c <- sweep(X, 2, colMeans(X), "-")
X_n <- sweep(X_c,2, apply(X_c,2,sd), "/")
exprSet=data.frame(X_n[,-c(1:2)],parkinson=Y)
sample_order<-order(Y)


test_error<-rep(0,times=nlambda)

for(i in 1:nlambda){
  test_error[i]<-readRDS(paste0("./_test_error/S_L_lambda_",i,".rds"))
  
  
}
#non_otu_predmatrix<-matrix(NA,327,2)
#for (i in 1:folds) {
 #non_otu_predmatrix[i,]<- exp(readRDS(paste0(result_path,"fold_",i,'/log_prior_test.rds')))
#}
#saveRDS(non_otu_predmatrix,paste0(result_path,"non_otu_predmatrix.rds"))
non_otu_predmatrix<-readRDS(paste0(result_path,"non_otu_predmatrix.rds"))

#nlambda <- length (lambda)
#n <- one_dataset$n
#pred_matrix <- array (0, dim = c(n,2,nlambda))
#for (i in 1:n) {
#for (j in 1:nlambda)
#pred_matrix[i,,j] <- readRDS(paste0(result_path,"fold_",i,"/S_L/lambda_",j,"Pr_of_PD_test.rds"))}
#saveRDS(pred_matrix, file = paste0(result_path, "all_pred_matrix.rds"))
pred_matrix<-readRDS(paste0(result_path, "all_pred_matrix.rds"))


#nv_lasso <-readRDS("~/ParkinsonBayes/dataset_real2_lasso_loocv_logbyT/_nv.lasso.rds")

n <- one_dataset$n
#pred_matrix_cov <- array (0, dim = c(n,2))



## find predictive metrics against lambda
y_true <- as.numeric(one_dataset$dataset$parkinson)
non_otu_er<-mean(y_true!=max.col(non_otu_predmatrix))
amlp_s_l <- rep(0, nlambda)
er_s_l <- rep(0, nlambda)

for (i in 1:nlambda) {
  eval <- evaluate_pred (pred_matrix[,,i],y_true)
  #auc_lasso[i]<-roc(y_true,pred_matrix[,2,i])$auc
  #auc_lasso[i] <- eval$auc
  er_s_l[i] <- eval$er
  amlp_s_l[i] <- eval$amlp
}
base<-matrix(0.5,nrow=327,ncol=2)
amlp_non_otu<-evaluate_pred (non_otu_predmatrix,y_true)$amlp
amlp_base<-evaluate_pred (base,y_true)$amlp

AUC<-vector()
AU_PR<-vector()
for (i in 1:nlambda){
  AUC[i]<-roc(y_true,pred_matrix[,2,i])$auc
  scores <- data.frame(pre=pred_matrix[,2,i], lab=y_true)
  AU_PR[i] <- pr.curve(scores.class0=scores[scores$lab=="2",]$pre,
                 scores.class1=scores[scores$lab=="1",]$pre,
                 curve=T)$auc.integral}
scores_nonotu <- data.frame(pre=exp(non_otu_predmatrix[,2]), lab=y_true)
AU_PR_non_otu <- pr.curve(scores.class0=scores_nonotu[scores_nonotu$lab=="2",]$pre,
                             scores.class1=scores_nonotu[scores_nonotu$lab=="1",]$pre,
                             curve=T)$auc.integral  

AUC_non_otu<-roc(y_true,non_otu_predmatrix[,2])$auc

which.min(er_s_l)
which.min(amlp_s_l)
best_lam.n<-which.min(er_s_l)

pred <- prediction(pred_matrix[,2,best_lam.n], y_true)
non_otu_pred<-prediction(non_otu_predmatrix[,2],y_true)
perf <- performance(pred,"tpr","fpr")
perf1 <- performance(pred, "prec", "rec")
non_otu_perf<-performance(non_otu_pred, "prec", "rec")

s_l_precision<-vector()
s_l_recall<-vector()
s_l_F_score<-vector()
for (i in 1:nlambda){
  pred<-data.frame(pre=max.col(pred_matrix[,,i]),obs=y_true)
  tp <- sum(pred$pre ==2 & pred$obs == 2)
  fp <- sum(pred$pre ==2 & pred$obs == 1)
  tn <- sum(pred$pre==1 & pred$obs == 1)
  fn <- sum(pred$pre==1 & pred$obs == 2)
  s_l_precision[i]<-tp/(tp+fp)
  s_l_recall[i]<-tp/(tp+fn)
  s_l_F_score[i]<-2/(1/s_l_precision[i]+1/s_l_recall[i])
  
}









###########AUC/AUPR#########
pdf(file=paste0("./__final_plots/L1SVM+logistic/AUC_and_AUC-PR.pdf"),width=6,height = 6)
par(mar=c(4,4,2,1))
par(mar=c(3,4,2,4))
plot(log(lambda), AUC,col="1",type="l",lwd=2,xlab=" ",ylab=" ",lty=1)
par(new=TRUE)
plot(log(lambda), AU_PR,col="2",type="l",lwd=2,axes=F,ylab=" ",xlab="",lty=4)
axis(side=4,col.axis="black")
mtext(side=1, line=2, expression(log(lambda)), cex=1)
mtext(side=4, line=2, "AUPRC", cex=1.2)
mtext(side=2, line=2, "AUC", cex=1.2)
#legend('topleft', c('AUPRC', 'AUC'),col=c('red', 'black'),bty = "o",cex=0.8,lty=1)
abline(h=seq(0.7,0.88,by=0.005),v=seq(-12,-4,by=1),col="lightgrey",lty=3)
legend('topleft', c('AUPRC', 'AUC'),col=c('red', 'black'),lty=c(4,1),bty = "o",cex=0.8)
dev.off()



#########Error raete#######
pdf(file=paste0("./__final_plots/L1SVM+logistic/L1SVM+logistic_Error_Rate.pdf"),width=6,height = 6)
par(mar=c(3,4,1,4))

plot(log(lambda), er_s_l,col="blue",type="l",lwd=2,ylim=c(0.275,0.40),xlab="", ylab="")
abline(h=non_otu_er,col="darkgreen",lty=2,lwd=2)
abline(h=130/327,col="red",lty=3,lwd=2)
axis(side = 4, at=seq(0,130/327,by=130/327/10), labels=seq(1,0,by=-0.1))
mtext(side=4, line=2, expression(R^2), cex=1.5,las=2)
mtext(side=2, line=2, "ER", cex=1.5,las=2)
mtext(side=1, line=2, expression(log(lambda)), cex=1)
abline(h=seq(0.2,0.4,by=0.005),v=seq(-12,-4,by=1),col="lightgrey",lty=3)
legend('bottomleft', c('OTU+age+sex', 'age+sex','baseline'),col  =c('blue', 'darkgreen','red'),cex=0.8,lty=c(1,2,3),lwd=2,bg="white")
##
dev.off()

#########amlp#######
pdf(file=paste0("./__final_plots/L1SVM+logistic/L1SVM+logistic_amlp.pdf"),width=6,height = 6)
par(mar=c(3,4,1,4))

plot(log(lambda), amlp_s_l,col="blue",type="l",lwd=2,xlab="", ylab="")
abline(h=amlp_non_otu,col="darkgreen",lty=2,lwd=2)
abline(h=amlp_base,col="red",lty=3,lwd=2)
#axis(side = 4, at=seq(0,130/327,by=130/327/10), labels=seq(1,0,by=-0.1))
#mtext(side=4, line=2, expression(R^2), cex=1.5,las=2)
mtext(side=2, line=2, "amlp", cex=1.5,las=3)
mtext(side=1, line=2, expression(log(lambda)), cex=1)
abline(h=seq(0.55,0.8,by=0.005),v=seq(-4.5,-1.5,by=1),col="lightgrey",lty=3)
legend('bottomleft', c('OTU+age+sex', 'age+sex','baseline'),col  =c('blue', 'darkgreen','red'),cex=0.8,lty=c(1,2,3),lwd=2,bg="white")
##
dev.off()


#plot(log(lambda), amlp_lasso,col="3",type="l")

### Probability of Parkinson######
#########Probability of Parkinson#######
pdf(file=paste0("./__final_plots/L1SVM+logistic/Probability of Parkinson.pdf"),width=6,height = 6)
par(mar=c(4,4,2,2))
par(mfrow=c(1,1))
PDorder<-order(y_true)
fr<-data.frame(p=pred_matrix[sample_order,,best_lam.n], l=y_true)
#plot(pred_matrix[sample_order,2,22],pch=2*y_true[sample_order]-1,col=-2*y_true[PDorder]+6,ylab="Probability of Parkinson",xlab="Subject ID")
#abline(h=0.5,lty=3,lwd=1)
plot(pred_matrix[sample_order,2,best_lam.n],c(1:327),pch=2*y_true[sample_order]-1,col=-2*y_true[PDorder]+6,xlab="Probability of Parkinson",ylab="Subject ID")
abline(v=0.5,lty=3,lwd=1)
dev.off()

#########ROC#########
pdf(file=paste0("./__final_plots/L1SVM+logistic/roc.pdf"),width=6,height = 6)
par(mar=c(4,4,2,2))
plot(roc(y_true,pred_matrix[,2,best_lam.n]),col="blue",main = "ROC",lty=1,lwd=2)
lines(roc(y_true,exp(non_otu_predmatrix[,2])),col='darkgreen',lty=2,lwd=2)
legend('bottomright', c('OTU+age+sex, AUC=0.78', 'age+sex, AUC=0.66'),col  =c('blue', 'darkgreen'),bty = "o",cex=0.8,lty=c(1,2),lwd=2)
dev.off()
######### P_R curve #########
pdf(file=paste0("./__final_plots/L1SVM+logistic/PRC.pdf"),width=6,height = 6)
par(mar=c(4,4,2,2))
plot(perf1,col="blue",main  = "P-R",ylim=c(0.6,1),lty=1,lwd=2)
lines(non_otu_perf@x.values[[1]],non_otu_perf@y.values[[1]],col="darkgreen",lty=2,lwd=2)
abline(h=197/327,lty=3,col="red",lwd=2)
legend('topright', c('OTU+age+sex, AUPR=0.83', 'age+sex, AUPR=0.72, ','baseline'),col  =c('blue', 'darkgreen','red'),bty = "o",cex=0.8,lty=c(1,2,3),lwd=2)
dev.off()

##############beta################
beta<-matrix(NA,ncol=382,nrow=327)
for(i in 1:n){
beta[i,]<-readRDS(paste0("./fold_",i,"/SVM/lambda_",best_lam.n,"beta.rds"))
}
colnames(beta)<-names(readRDS(paste0("./fold_",i,"/SVM/lambda_",best_lam.n,"beta.rds")))
median_absolute<-apply(beta, 2, function(x) median(abs(x)))
medians<-apply(beta, 2, function(x) median(x))
non_zero_n<-sum(medians!=0)
non_zero_index<-which(medians!=0)
ab_me_order.n<-order(median_absolute,decreasing=TRUE)
ab_me_order<-colnames(beta)[ab_me_order.n]
sign_of_top_beta<-sign(medians)[ab_me_order.n[1:non_zero_n]]
which(sign_of_top_beta==-1)
top_number<-2
rank<-ab_me_order.n
ordered_exprSet <- exprSet[sample_order, c("parkinson",ab_me_order)]
top_genera<-readRDS("./_top_genera_names.rds")
b_data<-data.frame(beta[,rank[1:top_number]],id=rownames(dataset))
beta_data<-melt(b_data,id.vars = c('id') )

pdf(file=paste0("./__final_plots/L1SVM+logistic/top 20 beta_absolute of L1SVM+logistic.pdf"),width=6,height = 6)
boxplot(abs(b_data[1:2]),col = -1*sign_of_top_beta+3,outcol=-1*sign_of_top_beta+3,names=rep("",time=2),par(mar=c(5,5,1,1)),xaxt='n')
labels <- top_genera[c(12,2)]
text(seq(1.25,2.25,by=1), par("usr")[3] - 0.001,  adj = 1, labels =labels, xpd = TRUE,cex=0.8)
mtext(side=2, line=3, "Absolute Values of Coefficients", cex=1)
abline (v=1.5, col = "lightgrey", lty=2,lwd=2)
grid()
legend('topright', c('-positive', '-negative'),fill =c('red', 'blue'),bty = "o",cex=0.8,bg="white")
#
dev.off()


################
#############
ordered_exprSet_add<-matrix(NA,nrow=328,ncol=382)
for (i in 1:382){ordered_exprSet_add[,i]<-c(ordered_exprSet[,i+1],NA)}
for (i in 1:382){ordered_exprSet_add[,i]<-as.numeric(ordered_exprSet_add[,i])}
parkinson_add<-c(as.numeric(ordered_exprSet[,1]),3)
parkinson_add<-as.factor(parkinson_add)
fit<-data.frame(parkinson_add,ordered_exprSet_add)
colnames(fit)<-colnames(ordered_exprSet)
#ordered_exprSet_add[,1]<-factor(ordered_exprSet_add[,1])
oadd_data1<-melt(fit[,1:(top_number+1)],id.vars = c("parkinson") )
#oadd_data2<-melt(fit[,1:(20+1)],id.vars = c("parkinson") )
#oadd_data3<-melt(fit[,1:(15+1)],id.vars = c("parkinson") )
#oadd_data4<-melt(fit[,1:(10+1)],id.vars = c("parkinson") )



######### boxplot pf normalized OTU ############
####################
pdf(file=paste0("./__final_plots/boxplot of top 20 normalized OTU vs parkinson.pdf"),width=6,height = 6)
boxplot(value~parkinson+variable, data=oadd_data1,col=rep(c(4,2,1),2),names=rep("",time=6),par(mar=c(5,5,1,1)),outcol=rep(c(4,2,1),2),xaxt="n")
labels <- top_genera[c(12,2)]
text(seq(2.5,8,by=3), par("usr")[3] - 0.2,  adj = 1, labels =labels, xpd = TRUE,cex=0.8)
mtext(side=2, line=3, "Standardized Log Relative Abundances", cex=1)
abline (v=1.5, col = "lightgrey", lty=2,lwd=2)
grid()
legend('topright', c('-positive', '-negative'),fill =c('red', 'blue'),bty = "o",cex=0.8,bg="white")
##
dev.off()
#####################

##################
pdf(file=paste0("./__final_plots/L1SVM+logistic/recall_precision_f.pdf"),width=6,height = 6)
par(mar=c(4,4,2,2))
par(mfrow = c(1, 3))
plot(s_l_recall,type="l",col=4,ylab="recall L1SVM+logistic")
plot(s_l_precision,type="l",col=5,ylab="precision L1SVM+logistic")

plot(s_l_F_score,type="l",col=6,ylab="F1 L1SVM+logistic")

dev.off()

pdf(file=paste0("./__final_plots/L1SVM+logistic/L1SVM+logistic_amlp.pdf"),width=6,height = 6)
par(mar=c(3,4,1,4))

plot(log(lambda), amlp_s_l,col="blue",type="l",lwd=2,xlab="", ylab="")
abline(h=amlp_non_otu,col="darkgreen",lty=2,lwd=2)
abline(h=amlp_base,col="red",lty=3,lwd=2)
#axis(side = 4, at=seq(0,130/327,by=130/327/10), labels=seq(1,0,by=-0.1))
#mtext(side=4, line=2, expression(R^2), cex=1.5,las=2)
mtext(side=2, line=2, "amlp", cex=1.5,las=3)
mtext(side=1, line=2, expression(log(lambda)), cex=1)
abline(h=seq(0.55,0.8,by=0.005),v=seq(-4.5,-1.5,by=1),col="lightgrey",lty=3)
legend('bottomleft', c('OTU+age+sex', 'age+sex','baseline'),col  =c('blue', 'darkgreen','red'),cex=0.8,lty=c(1,2,3),lwd=2,bg="white")
##
dev.off()




