#!/usr/bin/env Rscript

# gpirt_path = "~/Documents/Github/OrdGPIRT"
# setwd(gpirt_path)
print("Loading stan")

R_path="~/R/x86_64-redhat-linux-gnu-library/4.0"
.libPaths(R_path)
library(rstan)
rstan_options(auto_write = TRUE)

args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

library(gpirt)
library(dplyr)
library(ggplot2)
library(stats)
library(haven)
library(ltm)

setwd("../TAPS")
source("load_TAPS.R")

# code na as 0 for stan to ignore
na_mask = is.na(gpirt_data)
gpirt_data[na_mask] = 0

SAMPLE_ITERS = 100
BURNOUT_ITERS = 100
SEED = 12345

THIN = 1
CHAIN = 1
stan_data <- list(n=n,
                  m=m,
                  horizon=horizon,
                  K=C,
                  sigma=0.1,
                  y=gpirt_data)

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

save.image(file='doirt_TAPS_2014.RData')

gpirt_data[na_mask] = NA

fit_params <- as.data.frame(fit)

SAMPLE_ITERS = SAMPLE_ITERS / THIN

samples = list()
samples[["theta"]] = array(array(0, SAMPLE_ITERS*n*horizon), 
                           c(SAMPLE_ITERS,n, horizon))
samples[["threshold"]] = array(array(0, SAMPLE_ITERS*(C+1)), 
                               c(SAMPLE_ITERS,(C+1)))
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
source("../OrdGPIRT/getprob_gpirt.R")
gpirt_iccs = array(array(0, length(xs)*m*horizon),
                   c(length(xs),m, horizon))

C=5

for (h in 1:horizon) {
  wave = unique(data$wave)[h]
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
  }
}

save.image(file='doirt_TAPS_2014.RData')