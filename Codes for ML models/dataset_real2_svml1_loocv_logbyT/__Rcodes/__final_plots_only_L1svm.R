

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


SVM_test_error<-rep(0,times=nlambda)

for(i in 1:nlambda){
  SVM_test_error[i]<-readRDS(paste0("./_test_error/SVM_lambda_",i,".rds"))
  
}
min(SVM_test_error)
best_lam.n<-which.min(SVM_test_error)
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



n <- one_dataset$n



## find predictive metrics against lambda
y_true <- as.numeric(one_dataset$dataset$parkinson)
non_otu_er<-mean(y_true!=max.col(non_otu_predmatrix))

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
top_number<-25
rank<-ab_me_order.n
ordered_exprSet <- exprSet[sample_order, c("parkinson",ab_me_order)]
top_genera<-readRDS("./_top_genera_names.rds")
b_data<-data.frame(beta[,rank[1:top_number]],id=rownames(dataset))
beta_data<-melt(b_data,id.vars = c('id') )


pdf(file=paste0("./__final_plots/L1SVM/top 25 beta_absolute of L1SVM.pdf"),width=6,height = 6)
boxplot(abs(b_data[1:25]),col = -1*sign_of_top_beta+3,outcol=-1*sign_of_top_beta+3,names=rep("",time=25),par(mar=c(5,5,1,1)),xaxt='n')
labels <- top_genera[1:25]
text(seq(1,25,by=1), par("usr")[3] - 0.004, srt = 25, adj = 1, labels =labels, xpd = TRUE,cex=0.6)
mtext(side=2, line=3, "Absolute Value Coefficients of Genura Features of L1SVM", cex=1)
abline (v=seq(1.5,25,by = 1), col = "lightgrey", lty=2,lwd=2)
grid()
legend('topright', c('-positive', '-negative'),fill =c('red', 'blue'),bty = "o",cex=0.8,bg="white")
dev.off()


#top 20 beta
pdf(file=paste0("./__final_plots/L1SVM/top 20 beta_absolute of L1SVM.pdf"),width=6,height = 6)
boxplot(abs(b_data[1:20]),col = -1*sign_of_top_beta+3,outcol=-1*sign_of_top_beta+3,names=rep("",time=20),par(mar=c(5,5,1,1)),xaxt='n')
labels <- top_genera[1:20]
text(seq(1,20,by=1), par("usr")[3] - 0.004, srt = 23, adj = 1, labels =labels, xpd = TRUE,cex=0.6)
mtext(side=2, line=3, "Absolute Value Coefficients of Genura Features of L1SVM", cex=1)
abline (v=seq(1.5,20,by = 1), col = "lightgrey", lty=2,lwd=2)
grid()
legend('topright', c('-positive', '-negative'),fill =c('red', 'blue'),bty = "o",cex=0.8,bg="white")
dev.off()
#top 10 beta
pdf(file=paste0("./__final_plots/L1SVM/top 10 beta_absolute of L1SVM.pdf"),width=6,height = 6)
boxplot(abs(b_data[1:10]),col = -1*sign_of_top_beta+3,outcol=-1*sign_of_top_beta+3,names=rep("",time=10),par(mar=c(5,5,1,1)),xaxt='n')
labels <- top_genera[1:10]
text(seq(1,10,by=1), par("usr")[3] - 0.003, srt = 21, adj = 1, labels =labels, xpd = TRUE,cex=0.6)
mtext(side=2, line=3, "Absolute Value Coefficients of Genura Features of L1SVM", cex=1)
abline (v=seq(1.5,10,by = 1), col = "lightgrey", lty=2,lwd=2)
grid()
legend('topright', c('-positive', '-negative'),fill =c('red', 'blue'),bty = "o",cex=0.8,bg="white")
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
oadd_data2<-melt(fit[,1:(20+1)],id.vars = c("parkinson") )
oadd_data3<-melt(fit[,1:(15+1)],id.vars = c("parkinson") )
oadd_data4<-melt(fit[,1:(10+1)],id.vars = c("parkinson") )



######### boxplot pf normalized OTU ############
####################

pdf(file=paste0("./__final_plots/L1SVM/boxplot of top 25 normalized OTU vs parkinson L1SVM.pdf"),width=6,height = 6)
boxplot(value~parkinson+variable, data=oadd_data1,col=rep(c(4,2,1),25),names=rep("",time=75),par(mar=c(5,5,1,1)),outcol=rep(c(4,2,1),25),xaxt="n")
## Create some text labels
labels <- top_genera[1:25]
## Plot x axis labels at default tick marks
text(seq(1.5,75,by=3), par("usr")[3] - 0.25, srt = 25, adj = 1, labels =labels, xpd = TRUE,cex=0.6)
## Plot x axis label at line 6 (of 7)
mtext(side=2, line=3, "Standardlized Log Relative Abundance of L1SVM", cex=1)
abline (v=seq(3,72,by = 3), col = "grey", lty=2,lwd=2)
grid()
legend('topright', c('-parkinson', '-control'),fill =c('red', 'blue'),bty = "o",cex=0.7,bg='white')
dev.off()

pdf(file=paste0("./__final_plots/L1SVM/boxplot of top 20 normalized OTU vs parkinson L1SVM.pdf"),width=6,height = 6)
boxplot(value~parkinson+variable, data=oadd_data2,col=rep(c(4,2,1),20),names=rep("",time=60),par(mar=c(5,5,1,1)),outcol=rep(c(4,2,1),20),xaxt="n")
labels <- top_genera[1:20]
text(seq(1.5,60,by=3), par("usr")[3] - 0.25, srt = 24, adj = 1, labels =labels, xpd = TRUE,cex=0.6)
mtext(side=2, line=3, "Standardlized Log Relative Abundance of L1SVM", cex=1)
abline (v=seq(3,57,by = 3), col = "grey", lty=2,lwd=2)
grid()
legend('topright', c('-parkinson', '-control'),fill =c('red', 'blue'),bty = "o",cex=0.8,bg='white')
dev.off()

pdf(file=paste0("./__final_plots/L1SVM/boxplot of top 15 normalized OTU vs parkinson L1SVM.pdf"),width=6,height = 6)
boxplot(value~parkinson+variable, data=oadd_data3,col=rep(c(4,2,1),15),names=rep("",time=45),par(mar=c(5,5,1,1)),outcol=rep(c(4,2,1),15),xaxt="n")
labels <- top_genera[1:15]
text(seq(1.5,45,by=3), par("usr")[3] - 0.25, srt = 23, adj = 1, labels =labels, xpd = TRUE,cex=0.6)
mtext(side=2, line=3, "Standardlized Log Relative Abundance of L1SVM", cex=1)
abline (v=seq(3,42,by = 3), col = "grey", lty=2,lwd=2)
grid()
legend('topright', c('-parkinson', '-control'),fill =c('red', 'blue'),bty = "o",cex=0.8,bg="white")
dev.off()

pdf(file=paste0("./__final_plots/L1SVM/boxplot of top 10 normalized OTU vs parkinson L1SVM.pdf"),width=6,height = 6)
boxplot(value~parkinson+variable, data=oadd_data4,col=rep(c(4,2,1),10),names=rep("",time=30),par(mar=c(5,5,1,1)),outcol=rep(c(4,2,1),10),xaxt="n")
labels <- top_genera[1:10]
text(seq(1.5,30,by=3), par("usr")[3] - 0.15, srt = 21, adj = 1, labels =labels, xpd = TRUE,cex=0.6)
mtext(side=2, line=3, "Standardlized Log Relative Abundance of L1SVM", cex=1)
abline (v=seq(3,27,by = 3), col = "grey", lty=2,lwd=2)
grid()
legend('topright', c('-parkinson', '-control'),fill =c('red', 'blue'),bty = "o",cex=0.8,bg="white")
dev.off()




