#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

SAMPLE_ITERS = 3000
BURNOUT_ITERS = 3000
THIN = 4
CHAIN = 1
TYPE = "GP"

R_path="~/R/x86_64-redhat-linux-gnu-library/4.0"
.libPaths(R_path)
library(rstan)

if (length(args)==0) {
  TRAIN_START_YEAR = 2000
  TRAIN_END_YEAR = 2015
  TEST_YEAR = 2020
  SEED = 1
}

if (length(args)==4){
  TRAIN_START_YEAR = as.integer(args[1])
  TRAIN_END_YEAR = as.integer(args[2])
  TEST_YEAR = as.integer(args[3])
  SEED = as.integer(args[4])
}

library(gpirt)
library(dplyr)
library(ggplot2)
library(stats)
source("supremecourt_data.R")

if(TYPE=="GP"){
  theta_os = 1
  theta_ls = 7
}else if(TYPE=="CST"){
  theta_os = 1
  theta_ls = 10*horizon
}else{
  theta_os = 1
  theta_ls = 0.1
}

gpirt_data_train = gpirt_data

# set random seed
set.seed(SEED)

# random drop one observation for every case from TRAIN_END_YEAR+1 to TEST_YEAR
for(h in (TRAIN_END_YEAR-2000+1+1):(TEST_YEAR-2000+1)){
  year = unique(data$term)[h]
  caseIds = as.character(unique(data[data$term==year,"caseId"]))
  for (j in 1:length(caseIds)){
    mask = !is.na(gpirt_data[,j,h])
    drop_justice = sample(which(mask==1), 1)
    gpirt_data_train[drop_justice, j, h] = NA
  }
}

gpirt_data_train[is.na(gpirt_data_train)] = 0

stan_data <- list(n=n,
                  m=m,
                  horizon=horizon,
                  K=C,
                  sigma=0.1,
                  y=gpirt_data_train)

gpirt_data_train[gpirt_data_train==0] = NA

# train stan model
fit <- stan(file = "doirt-synthetic.stan",
            data = stan_data, 
            warmup = BURNOUT_ITERS, 
            iter = BURNOUT_ITERS + SAMPLE_ITERS, 
            chains = CHAIN, 
            cores = 1, 
            thin = THIN,
            control=list(adapt_delta=.98, max_treedepth = 15),
            seed = SEED,
            refresh=1
)

fit_params <- as.data.frame(fit)

SAMPLE_ITERS = SAMPLE_ITERS / THIN

samples = list()
samples[["theta"]] = array(array(0, SAMPLE_ITERS*n*horizon), 
                           c(SAMPLE_ITERS,n, horizon))
samples[["threshold"]] = array(array(0, SAMPLE_ITERS*(C+1)), 
                               c(SAMPLE_ITERS,(C+1)))

xs = seq(-5,5,0.01)
idx = 1:1001
for (it in 1:SAMPLE_ITERS) {
  samples[["fstar"]][[it]] =  array(array(0, length(xs[idx])*m*horizon), 
                                    c(length(xs[idx]),m, horizon))
  for(j in 1:m){
    for(h in 1:horizon){
      samples[["fstar"]][[it]][,j,h] = xs[idx]*fit_params[[paste("beta[",j,",",h,"]",sep="")]][it]
    }
  }
  for(i in 1:n){
    for(h in 1:horizon){
      samples[["theta"]][it,i,h] = fit_params[[paste("theta[",i,",",h,"]",sep="")]][it]
    }
  }
}

xs = seq(-5,5,0.01)
pred_theta = matrix(0, nrow=n, ncol=horizon)
pred_theta_sd = matrix(0, nrow=n, ncol=horizon)
for(i in 1:n){
  for (h in 1:horizon) {
    tmp = samples$theta[-1,i,h]
    pred_theta[i,h] = mean(tmp)
    pred_theta_sd[i,h] = sd(tmp)
  }
}

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

source("getprob_gpirt.R")

for(h in 1:horizon) {
  h_ = h+1999-TRAIN_END_YEAR
  year = unique(data$term)[h]
  caseIds = as.character(unique(data[data$term==year,"caseId"]))
  for (j in 1:length(caseIds)){
    IRFs = matrix(0, nrow=SAMPLE_ITERS, ncol=length(xs))
    for(iter in 1:SAMPLE_ITERS){
      IRFs[iter, ] = samples$fstar[[iter]][, j, h]
    }
    thresholds = array(array(0, SAMPLE_ITERS*(C+1)), 
                           c(SAMPLE_ITERS,(C+1)))
    thresholds[,1] = -Inf
    thresholds[,C+1] = Inf
    for(it in 1:SAMPLE_ITERS){
      for(c in 1:(C-1)){
        thresholds[it,1+c] = fit_params[[paste("alpha[",j,",",h,",",c,"]",sep="")]][it]
      }
    }
    probs = getprobs_gpirt(xs, t(IRFs), thresholds)
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

save(gpirt_data_train, gpirt_data, pred_theta,pred_theta_ll,pred_theta_sd,train_lls,
     train_acc, train_response, train_prediction,test_lls,
     test_acc, test_response, test_prediction,
     file=paste("./results/doirt_SupremeCourt_holdout_SEED_", SEED, ".RData" , sep=""))
