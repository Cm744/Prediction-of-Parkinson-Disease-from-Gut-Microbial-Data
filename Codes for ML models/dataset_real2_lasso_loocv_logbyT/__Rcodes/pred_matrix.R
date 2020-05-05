
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

if (!exists("irep")) irep <- 1
i <- irep #folds=327

pred_matrix <- array (0, dim = c(327,2,200))

  for (j in 1:200){
       pred_matrix[i,,j] <- readRDS(paste0(result_path,"fold_",i,"/lambda_",j,"test_prob.rds"))}
saveRDS(pred_matrix, file = paste0(result_path, "all_pred_matrix.rds"))

pred_matrix_cov <- array (0, dim = c(n,2))
pred_matrix_cov[i,] <- exp(readRDS(paste0(result_path,'fold_',i,'/log_prior_test','.rds')))
saveRDS(pred_matrix_cov, file = paste0(result_path, "pred_matrix-cov.rds"))