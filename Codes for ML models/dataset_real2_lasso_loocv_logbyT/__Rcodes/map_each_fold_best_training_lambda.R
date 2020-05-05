

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
library('gplots')
library('limma')
library('WGCNA')
#dataset[,1:382]<-apply(dataset[,1:382],1,function(x)(log(x+1)-log(sum(x)+length(x))) )
lambda <- scan("/home/xiw378/canola/lassoglmm/realdata/ethnicity/lassoout/lambda.txt")
nlambda<-length(lambda)
genus_names<-one_dataset$genus_names
nsample <- n

Y_name <- "parkinson"
X_name <- otu_names
X_cov<-c("age","sex")
fit_formula <- as.formula(paste0(Y_name, "~", paste(c(X_cov, X_name), collapse = "+")))
X <- model.matrix(fit_formula,data=dataset)[,-1]

Y <- dataset$parkinson
if (!exists("irep")) irep <- 1
i <- irep #folds=327



X_tr <- X[-i, ]
Y_tr <- as.numeric(Y[-i])
X_ts <- X[i, ]
Y_ts <- as.numeric(Y[i])


exprSet=t(X_tr[,-c(1:2)]) 


group <- Y[-i][order(Y[-i],decreasing = T)]

train_er<-vector()
for(j in 1:nlambda){train_er[j]<-readRDS(paste0("./fold_",i,"/lambda_",j,"train_error_rate.rds"))}

min_train_er_index<-which.min(train_er)
shownum<-readRDS(paste0("./fold_",i,"/lambda_",min_train_er_index,"num_otu_retained.rds"))
index_otu_retained<-readRDS(paste0("./fold_",i,"/lambda_",min_train_er_index,"index_otu_retained.rds"))
beta<-readRDS(paste0("./fold_",i,"/lambda_",min_train_er_index,"beta.rds"))[index_otu_retained]

heatmapData <- exprSet

heatmapData2<-heatmapData[index_otu_retained,]
hmExp=log10(heatmapData2+1)

mergedColors = as.character(labels2colors(group,colorSeq=c("gray","blue")))

hmMat=as.matrix(t(hmExp))
#
pdf(file=paste0('./fold_',i,'/_map_best_training_lambda.pdf'))
par(oma=c(2,2,2,3))

heatmap.2(hmMat,col='greenred',  trace="none", cexCol=1, dendrogram="col",Rowv=F,
          RowSideColors = mergedColors)

par(mfrow = c(1, 1))
par(las=2)
barplot(beta,names.arg=otu_names[index_otu_retained],cex.names = 0.5)
dev.off()

name<-paste0(otu_names[index_otu_retained],'=',genus_names[index_otu_retained])
write.table(name,file = paste0('./fold_',i,'/_remain_genus.names.txt'))

#dev.off()
