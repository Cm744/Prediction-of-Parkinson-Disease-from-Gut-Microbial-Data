#/bin/sh
#PBS -o /home/mac744/ParkinsonBayes/dataset_real2_svmela_loocv_logbyT/__Rcodes/outs/1_fit_logistic_regression_for_prior_logp.R.o$ifold
#PBS -e /home/mac744/ParkinsonBayes/dataset_real2_svmela_loocv_logbyT/__Rcodes/errors/1_fit_logistic_regression_for_prior_logp.R.e$ifold
#PBS -l nodes=1:ppn=1
#PBS -l mem=2GB
#PBS -l walltime=48:00:00
cd /home/mac744/ParkinsonBayes/dataset_real2_svmela_loocv_logbyT/__Rcodes
echo "Current working directory is `pwd`"
echo "Starting run at: `date`"
R --no-save --no-restore  <<EOI
irep <- $ifold
source("1_fit_logistic_regression_for_prior_logp.R", echo = TRUE)
EOI
echo "Program finished with exit code $? at: `date`"
exit 0
