




















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
library(reshape2)
library(ggplot2)

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
X_tr_c <- sweep(X_tr, 2, colMeans(X_tr), "-")
X_tr_n <- sweep(X_tr_c,2, apply(X_tr_c,2,sd), "/")
exprSet=data.frame(X_tr_n[,-c(1:2)],parkinson=Y[-i])

sample_order<-order(Y[-i])
#selected_samples<-Y[-i][sample_order[c(1:20)]]
#group <- sort(selected_samples)


test_error<-vector()

for(m in 1:nlambda){
  test_error[m]<-readRDS(paste0("./_test_error/lambda_",m,".rds"))
  
  
}


min_train_er_index<-which.min(test_error)#156
shownum<-readRDS(paste0("./fold_",i,"/lambda_",min_train_er_index,"num_otu_retained.rds"))
index_otu_retained<-readRDS(paste0("./fold_",i,"/lambda_",min_train_er_index,"index_otu_retained.rds"))
beta<-readRDS(paste0("./fold_",i,"/lambda_",min_train_er_index,"beta.rds"))

OTU_order.n <- order(-abs(beta))
OTU_order <- names(beta)[order(-abs(beta))]

write.table(cbind(beta[OTU_order.n],one_dataset$genus_names[OTU_order.n]),
            file = sprintf("%sfold_%s/ordered.otu.txt",result_path, i), row.names = FALSE,col.names = FALSE)

plot(abs(beta[OTU_order[1:200]]), type = "h")


#pid <- which (colnames(exprSet)=="parkinson")
ordered_exprSet <- exprSet[sample_order, c("parkinson",OTU_order)]

par(mfrow=c(1,2))
i_otu <- 4
plot(ordered_exprSet[,OTU_order[c(i_otu)]], col=as.numeric(ordered_exprSet$parkinson)+1)
boxplot(ordered_exprSet[,OTU_order[i_otu]]~as.numeric(ordered_exprSet$parkinson))



top_ordered_exprSet <- ordered_exprSet[c(1:10, 129+1:10), 2:20]


melted_top_ordered_exprSet <- melt(t(top_ordered_exprSet), id.vars=colnames(top_ordered_exprSet))


colnames(melted_top_ordered_exprSet)<-c("ID","Class","value")



p <- ggplot(melted_top_ordered_exprSet, aes(x="ID",y="Class")) 


p <- p + geom_tile(aes(fill=value))+scale_fill_gradient(name = "value",low = "#FFFFFF",high = "#012345")

p

pdf(file=paste0('./fold_',i,'/_map_ggplot.pdf'))

plot(abs(beta[OTU_order[1:200]]), type = "h")
p
dev.off()


#X_tr <- X[-i, ]
#df<-data.frame(X_tr[,-c(1:2)],parkinson=Y[-i])
#df1<-df[order(df$parkinson),]
#r<-paste0(df$parkinson[order(df$parkinson)],1:326)
#df2<-data.frame(df1[,-383],row.names = r)
#df3<-melt(t(df2), id.vars=colnames(dataset))
#colnames(df3)<-c("ID","class","value")