SAMPLE_ITERS = 500
BURNOUT_ITERS = 500
THIN = 4
CHAIN = 1

print("Loading stan")

R_path="~/R/x86_64-redhat-linux-gnu-library/4.0"
.libPaths(R_path)
library(rstan)
rstan_options(auto_write = TRUE)

args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

if (length(args)==0) {
  TRAIN_START_YEAR = 1
  TRAIN_END_YEAR = 31
  TEST_YEAR = 41
  SEED = 1
  DROP_RATIO = 10
  TYPE = "Matern"
  TYPE = "NIRT"
}

if (length(args)==6){
  TRAIN_START_YEAR = as.integer(args[1])
  TRAIN_END_YEAR = as.integer(args[2])
  TEST_YEAR = as.integer(args[3])
  SEED = as.integer(args[4])
  DROP_RATIO = as.integer(args[5])
  TYPE = args[6]
}

library(gpirt)
library(dplyr)
library(ggplot2)
library(stats)
library(ltm)

source("load_TAPS.R")

gpirt_data_train = gpirt_data

# set random seed
set.seed(SEED)

# random drop one observation for every case from TRAIN_END_YEAR+1 to TEST_YEAR
for(h in (TRAIN_END_YEAR+1):(TEST_YEAR)){
  for (j in 1:m){
    mask = !is.na(gpirt_data[,j,h])
    if (sum(mask==1)>0){
      drop_unit = sample(which(mask==1), as.integer(DROP_RATIO/100*length(mask)))
    }
    gpirt_data_train[drop_unit, j, h] = NA
  }
}

theta_os = 1
theta_ls = 12 # length scale is set to a year

if(TYPE=="LIRT"){
  TYPE1 = TYPE
  TYPE = "RBF"
  theta_os = 1
  theta_ls = 10*horizon
}

if(TYPE=="NIRT"){
  TYPE1 = TYPE
  TYPE = "RBF"
  theta_os = 1
  theta_ls = 0.1
}

beta_prior_sds =  matrix(3.0, nrow = 3, ncol = ncol(gpirt_data))
theta_prior_sds =  matrix(1.0, nrow = 2, ncol = nrow(gpirt_data))
theta_prior_sds[2,] = 0
beta_prior_sds[1,] = 0.0
beta_prior_sds[3,] = 0.0
samples_all <- gpirtMCMC(gpirt_data_train, SAMPLE_ITERS,BURNOUT_ITERS,
                         THIN, CHAIN, theta_init = theta_init,
                         beta_prior_sds = beta_prior_sds,
                         theta_prior_sds = theta_prior_sds,
                         theta_os = theta_os, theta_ls = theta_ls,
                         vote_codes = NULL, thresholds=NULL,
                         SEED=SEED, constant_IRF = 1, KERNEL=TYPE)

TYPE = TYPE1

SAMPLE_ITERS = SAMPLE_ITERS/THIN
samples = samples_all[[1]]

xs = seq(-5,5,0.01)

for(it in 1:SAMPLE_ITERS){
  for(h in 1:horizon){
    for(j in 1:m){
      samples$f[[it]][,j,h] = samples$f[[it]][,j,h] + samples$beta[[it]][1,j,h] + samples$beta[[it]][2,j,h]*samples$theta[it,,h]
      samples$fstar[[it]][,j,h] = samples$fstar[[it]][,j,h] + samples$beta[[it]][1,j,h] + samples$beta[[it]][2,j,h]*xs
    }
  }
}

pred_theta = matrix(0, nrow=n, ncol=horizon)
pred_theta_sd = matrix(0, nrow=n, ncol=horizon)
for (h in 1:horizon) {
  mask = rep(0, SAMPLE_ITERS)
  for(iter in 1:SAMPLE_ITERS){
    if(cor(samples$theta[1,,h],  samples$theta[iter,,h])>0){
      mask[iter] = 1
    }
  }
  for(i in 1:n){
    if(sum(mask)>=sum(!mask)){
      tmp = samples$theta[mask==1,i,h]
    }
    else{
      tmp = samples$theta[mask==0,i,h]
    }
    pred_theta[i,h] = mean(tmp)
    pred_theta_sd[i,h] = sd(tmp)
  }
}

xs = seq(-5,5,0.01)
source("../getprob_gpirt.R")
gpirt_iccs = array(array(0, length(xs)*m*horizon),
                   c(length(xs),m, horizon))

C=5

# train/test statistics
# train
train_lls = c()
train_acc = c()
train_response = c()
train_prediction = c()

# test
test_lls = vector(mode='list', length=(TEST_YEAR-TRAIN_END_YEAR))
test_acc = vector(mode='list', length=(TEST_YEAR-TRAIN_END_YEAR))
test_response = vector(mode='list', length=(TEST_YEAR-TRAIN_END_YEAR))
test_prediction = vector(mode='list', length=(TEST_YEAR-TRAIN_END_YEAR))

for (h in 1:horizon) {
  wave = unique(data$wave)[h]
  h_ = h-TRAIN_END_YEAR
  for (j in 1:m) {
    IRFs = matrix(0, nrow=SAMPLE_ITERS, ncol=length(xs))
    mask = rep(0, SAMPLE_ITERS)
    for(iter in 1:SAMPLE_ITERS){
      IRFs[iter, ] = samples$fstar[[iter]][, j, h]
      if (cor(samples$fstar[[1]][, j, h],samples$fstar[[iter]][, j, h])>0){
        mask[iter] = 1
      }
    }
    thresholds = matrix(0,nrow=SAMPLE_ITERS,ncol=C+1)
    for(iter in 1:SAMPLE_ITERS){
      thresholds[iter, ] = samples$threshold[[iter]][j,,h]
    }
    if(sum(mask)>=sum(!mask)){
      probs = getprobs_gpirt(xs, t(IRFs[mask==1,]), thresholds[mask==1,])
    }
    else{
      probs = getprobs_gpirt(xs, t(IRFs[mask==0,]), thresholds[mask==0,])
    }
    tmp = probs %>%
      group_by(xs) %>%
      summarize(icc=sum(order*p))
    gpirt_iccs[,j,h] = tmp$icc
    
    # test/train statistic
    for (i in 1:n) {
      if(!is.na(gpirt_data_train[i,j,h]) & !is.na(gpirt_data[i,j,h])){
        # train
        pred_idx = 1+as.integer((pred_theta[i,h]+5)*100)
        ll = log(probs$p[probs$xs==xs[pred_idx]])
        y_pred = which.max(ll)
        train_acc = c(train_acc, y_pred==(gpirt_data[i,j, h]))
        train_lls = c(train_lls, ll[gpirt_data[i,j, h]])
        train_response = c(train_response, gpirt_data[i,j, h])
        train_prediction = c(train_prediction, y_pred)
      }
      if(is.na(gpirt_data_train[i,j,h]) & !is.na(gpirt_data[i,j,h])){
        # test
        pred_idx = 1+as.integer((pred_theta[i,h]+5)*100)
        ll = log(probs$p[probs$xs==xs[pred_idx]])
        y_pred = which.max(ll)
        test_acc[[h_]] = c(test_acc[[h_]], y_pred==(gpirt_data[i,j,h]))
        test_lls[[h_]] = c(test_lls[[h_]], ll[gpirt_data[i,j, h]])
        test_response[[h_]] = c(test_response[[h_]],gpirt_data[i,j, h])
        test_prediction[[h_]] = c(test_prediction[[h_]], y_pred)
      }
    }
  }
}


print("gpirt finished!")

if(TYPE=="RBF"){
  file_name = paste("./results/SEirt_TAPS_holdout_", "DR_", DROP_RATIO, "_SEED_", SEED, ".RData" , sep="")
}

if(TYPE=="NIRT"){
  file_name = paste("./results/Nirt_TAPS_holdout_", "DR_", DROP_RATIO, "_SEED_", SEED, ".RData" , sep="")
}

if(TYPE=="LIRT"){
  file_name = paste("./results/Lirt_TAPS_holdout_", "DR_", DROP_RATIO, "_SEED_", SEED, ".RData" , sep="")
}

if(TYPE=="Matern"){
  file_name = paste("./results/gpirt_TAPS_holdout_SEED_", "DR_", DROP_RATIO, "_", SEED, ".RData" , sep="")
}

save(gpirt_data_train, gpirt_data, pred_theta,pred_theta_sd,train_lls,
     train_acc, train_response, train_prediction,test_lls,
     test_acc, test_response, test_prediction,
     file=file_name)