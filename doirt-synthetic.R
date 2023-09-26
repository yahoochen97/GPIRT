#!/usr/bin/env Rscript

# gpirt_path = "~/Documents/Github/OrdGPIRT"
# setwd(gpirt_path)
print("Loading stan")

R_path="~/R/x86_64-redhat-linux-gnu-library/4.0"
.libPaths(R_path)
library(rstan)

args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

if (length(args)==0) {
  SEED = 1
  C = 2
  n = 100
  m = 10
  horizon = 10
  TYPE = "GP"
  CONSTANT_IRF = 0
}

if (length(args)==7){
  SEED = as.integer(args[1])
  C = as.integer(args[2])
  n = as.integer(args[3])
  m = as.integer(args[4])
  horizon = as.integer(args[5])
  TYPE = args[6]
  CONSTANT_IRF = as.integer(args[7])
}

library(dplyr)
library(stats)

source("getprob_gpirt.R")
HYP = paste("GP_C_", C, '_n_', n, '_m_', m, '_h_', horizon,'_CSTIRF_', CONSTANT_IRF , '_SEED_', SEED, sep="")
print(HYP)
load(file=paste("./data/", HYP, ".RData" , sep=""))
HYP = paste("BRW_C_", C, '_n_', n, '_m_', m, '_h_', horizon,'_CSTIRF_', CONSTANT_IRF , '_SEED_', SEED, sep="")

SAMPLE_ITERS = 500
BURNOUT_ITERS = 500

data_train[is.na(data_train)] = 0

THIN = 4
CHAIN = 1
stan_data <- list(n=n,
                  m=m,
                  horizon=horizon,
                  K=C,
                  sigma=0.1,
                  y=data_train)

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

data_train[data_train==0] = NA

fit_params <- as.data.frame(fit)

SAMPLE_ITERS = SAMPLE_ITERS / THIN

samples = list()
samples[["theta"]] = array(array(0, SAMPLE_ITERS*n*horizon), 
                           c(SAMPLE_ITERS,n, horizon))
# samples[["threshold"]] = array(array(0, SAMPLE_ITERS*(C+1)), 
#                                c(SAMPLE_ITERS,(C+1)))
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

pred_theta = matrix(0, nrow=nrow(data), ncol=dim(data)[3])
pred_theta_sd = matrix(0, nrow=nrow(data), ncol=dim(data)[3])
# 0 to keep, 1 to drop
# drop_wrong_signs = array(array(0, n*horizon*SAMPLE_ITERS), 
#                          c(n, horizon, SAMPLE_ITERS))
for(i in 1:n){
  for (h in 1:horizon) {
    tmp = samples$theta[-1,i,h]
    pred_theta[i,h] = mean(tmp)
    pred_theta_sd[i,h] = sd(tmp)
  }
}

pred_theta_ll = matrix(0, nrow=nrow(data), ncol=dim(data)[3])

for(i in 1:n){
  for (h in 1:horizon) {
    pred_theta_ll[i,h] = log(dnorm(theta[i,h]*sign(cor(theta[,h],pred_theta[,h])), 
                                   mean=pred_theta[i,h], 
                                   sd=pred_theta_sd[i,h]))
  }
}

ordinal_lls = function(f, thresholds){
  result = c()
  for (c in 1:(length(thresholds)-1)) {
    z1 = thresholds[c] - f;
    z2 = thresholds[c+1] - f;
    result = c(result,log(pnorm(z2)-pnorm(z1)));
  }
  return(result)
}

pred_lls = c()
pred_acc = c()
train_lls = c()
train_acc = c()

for (i in 1:n) {
  for (j in 1:m) {
    for(h in 1:horizon){
      if(!is.na(data[[i,j,h]])){
        lls = matrix(0,nrow=SAMPLE_ITERS, ncol = C)
        y_pred = rep(0, SAMPLE_ITERS)
        pred_idx = 1+as.integer((pred_theta[i,h]+5)*100)
        for (iter in 1:SAMPLE_ITERS) {
          f_pred = samples$fstar[[iter]][pred_idx, j, h]
          ll = ordinal_lls(f_pred, samples$threshold[[iter]][j,,h])
          lls[iter,] = ll
          y_pred[iter] =  which.max(ll)
        }
        ll = log(exp(lls))
        y_pred = round(mean(y_pred))
        if(train_idx[i,j, h]==0){
          pred_acc = c(pred_acc, y_pred==(data[[i,j, h]]))
          pred_lls = c(pred_lls, ll[data[[i,j, h]]])
        }else{
          train_acc = c(train_acc, y_pred==(data[[i,j, h]]))
          train_lls = c(train_lls, ll[data[[i,j, h]]])
        }
      }
    }
  }
}

# cor of icc
idx = 1:1001
# gpirt_iccs = colMeans(IRFs, dim=1)
gpirt_iccs = array(array(0, length(xs[idx])*m*horizon), 
                   c(length(xs[idx]),m, horizon))
true_iccs = array(array(0, length(xs[idx])*m*horizon), 
                  c(length(xs[idx]),m, horizon))
cor_icc = matrix(0, nrow=m, ncol=horizon)
rmse_icc = matrix(0, nrow=m, ncol=horizon)

for (h in 1:horizon) {
  for (j in 1:m) {
    source("true_irf.R")
    # sign(cor(theta[,h],pred_theta[,h]))
    probs = getprobs_gpirt(sign(cor(theta[,h],pred_theta[,h]))*xs[idx], 
                           irfs, matrix(thresholds,nrow=1))
    tmp = probs %>% 
      group_by(xs) %>%
      summarize(icc=sum(order*p))
    true_iccs[,j,h] = tmp$icc
    IRFs = matrix(0, nrow=SAMPLE_ITERS, ncol=length(idx))
    for(iter in 1:SAMPLE_ITERS){
      IRFs[iter, ] = samples$fstar[[iter]][, j, h]
    }
    thresholds_tmp = array(array(0, SAMPLE_ITERS*(C+1)), 
                           c(SAMPLE_ITERS,(C+1)))
    thresholds_tmp[,1] = -Inf
    thresholds_tmp[,C+1] = Inf
    for(it in 1:SAMPLE_ITERS){
      for(c in 1:(C-1)){
        thresholds_tmp[it,1+c] = fit_params[[paste("alpha[",j,",",h,",",c,"]",sep="")]][it]
      }
    }
    probs = getprobs_gpirt(xs[idx], t(IRFs), 
                           thresholds_tmp)
    tmp = probs %>% 
      group_by(xs) %>%
      summarize(icc=sum(order*p))
    gpirt_iccs[,j,h] = tmp$icc
    cor_icc[j,h] = cor(gpirt_iccs[,j,h], true_iccs[,j,h])
    rmse_icc[j,h] = sqrt(mean((gpirt_iccs[,j,h]-true_iccs[,j,h])^2))
    
    # test/train statistic
    # for (i in 1:n) {
    #   if(!is.na(data_train[i,j,h])){
    #     # train
    #     pred_idx = 1+as.integer((pred_theta[i,h]+5)*100)
    #     ll = log(probs$p[probs$xs==xs[pred_idx]])
    #     y_pred = which.max(ll)
    #     train_acc = c(train_acc, y_pred==data[i,j, h])
    #     train_lls = c(train_lls, ll[data[i,j, h]])
    #   }else{
    #     # test
    #     pred_idx = 1+as.integer((pred_theta[i,h]+5)*100)
    #     ll = log(probs$p[probs$xs==xs[pred_idx]])
    #     y_pred = which.max(ll)
    #     pred_acc = c(pred_acc, y_pred==data[i,j,h])
    #     pred_lls = c(pred_lls, ll[data[i,j, h]])
    #   }
    # }
  }
}

cor_theta = c()
for (i in 1:horizon) {
  cor_theta = c(cor_theta, cor(theta[,i], pred_theta[,i]))
}

plot(xs[idx],gpirt_iccs[,1,1], ylim=c(1,2))
mask = is.na(data_train[,1,1])
points(pred_theta[!mask,1], data_train[!mask,1,1])

theta_rhats = c(1,1,1)

print("Average Rhat:")
print(mean(theta_rhats))
print(max(theta_rhats))
print(mean(abs(cor_theta)))
print(mean(pred_theta_sd))
print(mean(pred_theta_ll))
print(mean(train_lls[!is.infinite(train_lls) & !is.na(train_lls)]))
print(mean(train_acc[!is.infinite(train_acc) & !is.na(train_acc)]))
print(mean(pred_lls[!is.infinite(pred_lls) & !is.na(pred_lls)]))
print(mean(pred_acc[!is.infinite(pred_acc) & !is.na(pred_acc)]))
print(mean(array(abs(cor_icc), n*horizon)))
print(mean(array(rmse_icc, n*horizon)))

save(gpirt_iccs, true_iccs, theta, pred_theta,pred_theta_ll,pred_theta_sd,train_lls,
     train_acc, pred_lls, pred_acc,cor_icc, rmse_icc, theta_rhats,
     file=paste("./results/gpirt_", HYP, ".RData" , sep=""))

quit()