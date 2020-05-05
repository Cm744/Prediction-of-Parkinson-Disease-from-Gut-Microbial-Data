library(caret)

result_path<-"~/ParkinsonBayes/dataset_real2_lasso_loocv_logbyT/"
#####################################################################
setwd(result_path)

ds=read.table(paste0(result_path,"genus_level.txt"),sep="\t",header=T, row.names=1,comment.char = "@",stringsAsFactors = FALSE,quote = "")
attach(ds)

parkinson<-as.factor(parkinson)
sex<-as.factor(sex)


otus.zinb<- data.frame(ds[,1:382])
genus_names<-colnames(ds[,1:382])
otu_names<-paste("genus",1:382,sep="")
otu.num<-ncol(otus.zinb)
colnames(otus.zinb)<-otu_names
Totalread<-apply(otus.zinb,1,sum)
n<-length(Totalread)
logTR<-log(Totalread)
dataset<-data.frame(otus.zinb, logTR,sex,age,parkinson)
folds.index<-createFolds(dataset[,'parkinson'],k=327)
one_dataset <- list(dataset = dataset, n = n , otu.num = otu.num, otu_names=otu_names,genus_names=genus_names,folds=327,folds.index=folds.index)
saveRDS(one_dataset, "./dataset_real2_loocv.rds")
#save(one_dataset,file="./dataset_real2_10folds.Rdata")

