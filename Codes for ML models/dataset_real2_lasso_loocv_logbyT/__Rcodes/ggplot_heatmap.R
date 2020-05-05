

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
TR_order<-order(TR[-i])


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
genura_name_order<-genus_names[order(-abs(beta))]
OTU_order <- names(beta)[order(-abs(beta))]

write.table(cbind(beta[OTU_order.n],OTU_order,one_dataset$genus_names[OTU_order.n]),
            file = sprintf("%sfold_%s/_ordered.otu.txt",result_path, i), row.names = FALSE,col.names = FALSE)

#plot(abs(beta[OTU_order[1:10]]), type = "h")


#pid <- which (colnames(exprSet)=="parkinson")
ordered_exprSet <- exprSet[sample_order, c("parkinson",OTU_order)]
pdf(file=paste0('./fold_',i,'/_ggplot_20.pdf'))
#pdf(file=paste0('./fold_',i,'/_ggplot_10.pdf'))
top_number<-20 #10
top_genera<-vector()
top_taxa<-one_dataset$genus_names[OTU_order.n][1:top_number]
for(i_otu in 1:top_number){
  genera<-as.vector(strsplit(top_taxa[i_otu],".",fixed = T)[[1]])
  top_genera[i_otu]<-paste0(genera[length(genera)-1],',',genera[length(genera)])
}
data<-melt(ordered_exprSet[,1:(top_number+1)],id.vars = c("parkinson") )
ggplot(data, aes(variable, value,fill=parkinson)) + geom_boxplot()+labs(x="genura",y="Standardlized value of OTU reads",title = "OTU vs parkinson")+scale_x_discrete(labels=top_genera)+theme(axis.text.x = element_text(face="bold",size=7, angle=90))

beta1<-data.frame(ID=names(beta[OTU_order[1:top_number]]),value=unname(beta[OTU_order[1:top_number]]),n=OTU_order.n[1:top_number],num=as.factor(c(1:top_number)))
beta2<- transform(beta1,judge = ifelse(value > 0,"Positive","Negative"))
ggplot(data = beta2, mapping = aes(x = num,y = abs(value),fill=judge)) + geom_bar(stat= 'identity')+scale_x_discrete(labels=top_genera)+theme(axis.text.x = element_text(face="bold",size=7, angle=90))+labs(x="genura",y="absolute value of lasso coefficients",title = "rank of absolute value of coefficients")+scale_fill_manual(values = c("darkred","blue"))

d<-data.frame(X_tr_n [,-c(1:2)][,c(OTU_order[1:top_number])],parkinson=Y[-i])
df<-d[TR_order[c(1:20)],]
df1<-df[order(df$parkinson),]
r<-paste0(df$parkinson[order(df$parkinson)],1:20)

df2<-data.frame(df1[,-(top_number+1)])
df2<-data.frame(df1[,-(top_number+1)],row.names = r)


df3<-melt(t(df2), id.vars=c(OTU_order[1:top_number]))
colnames(df3)<-c("ID","class","value")

p <- ggplot(df3, aes(x=class,y=ID))


p <- p + geom_tile(aes(fill=value))+scale_fill_gradient2("legend name",low = "white",high = "red")


print(p)
dev.off()
pdf(file=paste0('./fold_',i,'/_box_and_scatter_plots_20.pdf'))
#pdf(file=paste0('./fold_',i,'/_box_and_scatter_plots_10.pdf'))
for(i_otu in c(1:top_number)){
par(mfrow=c(1,2))

plot(ordered_exprSet[,OTU_order[c(i_otu)]], col=as.numeric(ordered_exprSet$parkinson)+1,ylab=OTU_order[i_otu])

boxplot(ordered_exprSet[,OTU_order[i_otu]]~ordered_exprSet$parkinson,ylab=OTU_order[i_otu],xlab="parkinson")


}
dev.off()





