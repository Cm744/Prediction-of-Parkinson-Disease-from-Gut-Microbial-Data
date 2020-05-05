

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
library(glmnet)
library('gplots')
library('limma')
library('WGCNA')
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

Y <- dataset$parkinson
if (!exists("irep")) irep <- 1
i <- irep #folds=327




X_tr=X[-i,]
#exprSet=t(X_tr[,-c(1:2)]) 
X_tr_c <- sweep(X_tr, 2, colMeans(X_tr), "-")
X_tr_n <- sweep(X_tr_c,2, apply(X_tr_c,2,sd), "/")
exprSet=data.frame(X_tr_n[,-c(1:2)])

TR_order<-order(TR)
selected_samples<-Y[-i][TR_order[c(1:20)]]
group <- sort(selected_samples)


test_error<-vector()

for(m in 1:nlambda){
  test_error[m]<-readRDS(paste0("./_test_error/lambda_",m,".rds"))
  
  
}


min_train_er_index<-which.min(test_error)#156
shownum<-readRDS(paste0("./fold_",i,"/lambda_",min_train_er_index,"num_otu_retained.rds"))
index_otu_retained<-readRDS(paste0("./fold_",i,"/lambda_",min_train_er_index,"index_otu_retained.rds"))
beta<-readRDS(paste0("./fold_",i,"/lambda_",min_train_er_index,"beta.rds"))

o.beta <- order(-abs(beta))

write.table(cbind(beta[o.beta],one_dataset$genus_names[o.beta]),
            file = sprintf("%sfold_%s/ordered.otu.txt",result_path, i), row.names = FALSE,col.names = FALSE)

plot(abs(beta[o.beta[1:200]]), type = "h")

heatmapData <- exprSet

#hmExp<-heatmapData[o.beta[1:20],]
hmExp0<-heatmapData[o.beta[1:10],TR_order[c(1:20)]]
#hmExp<-hmExp0[,order(TR)[1:20]]

mergedColors = as.character(labels2colors(group,colorSeq=c("gray","blue")))

hmMat=as.matrix(t(hmExp0))
#
#pdf(file=paste0('./fold_',i,'/_map_best_test_lambda.pdf'))
par(oma=c(2,2,2,3))

heatmap.2(hmMat,col='greenred',  trace="none", cexCol=1, dendrogram="col",Rowv=F,
          RowSideColors = mergedColors)
heatmap.2(hmMat,col='greenred',key=T,keysize = 0.7,
          trace="none", cexCol=0.6, cexRow=2.5,dendrogram="col",Rowv=F,
          RowSideColors = mergedColors,scale="col")

#par(mfrow = c(1, 1))
#par(las=2)
#barplot(beta,names.arg=otu_names[index_otu_retained],cex.names = 0.5,border=T)
#dev.off()

name<-paste0(otu_names[index_otu_retained],'=',genus_names[index_otu_retained])
write.table(name,file = paste0('./fold_',i,'/_remain_genus_names2.txt'))

#dev.off()
