#!/usr/bin/env Rscript

# gpirt_path = "~/Documents/Github/OrdGPIRT"
# setwd(gpirt_path)
# setwd("../TAPS")
SAMPLE_ITERS = 1000
BURNOUT_ITERS = 1000
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
  TRAIN_END_YEAR = 42
  TEST_YEAR = 41
  DROP_RATIO = 0
  SEED = 12345
}

if (length(args)==5){
  TRAIN_START_YEAR = as.integer(args[1])
  TRAIN_END_YEAR = as.integer(args[2])
  TEST_YEAR = as.integer(args[3])
  SEED = as.integer(args[4])
  DROP_RATIO = as.integer(args[5])
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
if(TRAIN_END_YEAR<=41){
  for(h in (TRAIN_END_YEAR+1):(TEST_YEAR)){
    for (j in 1:m){
      mask = !is.na(gpirt_data[,j,h])
      if (sum(mask==1)>0){
        drop_unit = sample(which(mask==1), as.integer(DROP_RATIO/100*length(mask)))
      }
      gpirt_data_train[drop_unit, j, h] = NA
    }
  }
}


# code na as 0 for stan to ignore
gpirt_data_train[is.na(gpirt_data_train)] = 0

stan_data <- list(n=n,
                  m=m,
                  horizon=horizon,
                  K=C,
                  sigma=0.1,
                  y=gpirt_data_train)

# train stan model
fit <- stan(file = "../doirt-synthetic.stan",
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

if(TRAIN_END_YEAR==42){
  save.image(file='doirt_TAPS_2014.RData')
}


# gpirt_data[na_mask] = NA

fit_params <- as.data.frame(fit)

SAMPLE_ITERS = SAMPLE_ITERS / THIN * CHAIN

samples = list()
samples[["theta"]] = array(array(0, SAMPLE_ITERS*n*horizon), 
                           c(SAMPLE_ITERS,n, horizon))
xs = seq(-5,5,0.01)
idx = 1:length(xs)
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
  
  samples[["threshold"]][[it]] =  array(array(0, m*(C+1)*horizon), 
                                    c(m, C+1, horizon))
  for(j in 1:m){
    for(h in 1:horizon){
      samples[["threshold"]][[it]][j,1,h] = -Inf
      samples[["threshold"]][[it]][j,C+1,h] = Inf
      for(c in 1:(C-1)){
        samples[["threshold"]][[it]][j,1+c,h] = fit_params[[paste("alpha[",j,",",h,",",c, "]",sep="")]][it]
      }
    }
  }
}

# summarize pred theta over iterations
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

# get icc
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

print("doirt finished!")

if(TRAIN_END_YEAR==42){
  save.image(file='doirt_TAPS_2014.RData')
}else{
  save(gpirt_data_train, gpirt_data, pred_theta,pred_theta_sd,train_lls,
       train_acc, train_response, train_prediction,test_lls,
       test_acc, test_response, test_prediction,
       file=paste("./results/doirt_TAPS_holdout_", "DR_", DROP_RATIO, "_SEED_", SEED, ".RData" , sep=""))
  
}


# save.image(file='doirt_TAPS_2014.RData')
