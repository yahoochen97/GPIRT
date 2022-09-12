#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

# gpirt_path = "~/Documents/Github/OrdGPIRT"
# setwd(gpirt_path)
# TYPE = "RDM"

# gpirt_path = "~/Documents/Github/gpirt"
# setwd(gpirt_path)
# library(Rcpp)
# Rcpp::compileAttributes()
# install.packages(gpirt_path, type="source", repos = NULL)#,lib=R_path, INSTALL_opts = '--no-lock')
# setwd("../OrdGPIRT")


if (length(args)==0) {
  SEED = 1
  C = 2
  n = 50
  m = 10
  horizon = 10
  TYPE = "GP"
  CONSTANT_IRF = 1
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

# install gpirt package
R_path="~/R/x86_64-redhat-linux-gnu-library/4.0"
.libPaths(R_path)
# options("install.lock"=FALSE)
# gpirt_path = "~/Documents/Github/gpirt"
# gpirt_path = "../gpirt"
# setwd(gpirt_path)
# library(Rcpp)
# Rcpp::compileAttributes()
# install.packages(gpirt_path, type="source", repos = NULL,lib=R_path, INSTALL_opts = '--no-lock')
# setwd("../OrdGPIRT")
library(gpirt)
library(dplyr)
library(stats)

source("getprob_gpirt.R")
HYP = paste("GP_C_", C, '_n_', n, '_m_', m, '_h_', horizon,'_CSTIRF_', CONSTANT_IRF , '_SEED_', SEED, sep="")
print(HYP)
load(file=paste("./data/", HYP, ".RData" , sep=""))
HYP = paste(TYPE, "_C_", C, '_n_', n, '_m_', m, '_h_', horizon,'_CSTIRF_', CONSTANT_IRF , '_SEED_', SEED, sep="")

SAMPLE_ITERS = 500
BURNOUT_ITERS = 500
if(TYPE=="GP"){
    theta_os = 1
    theta_ls = as.integer(horizon/2)
}else if(TYPE=="CST"){
# constant model
    theta_os = 1
    theta_ls = 10*horizon
}else{
# random model
    theta_os = 1
    theta_ls = 0.1
}

THIN = 4
CHAIN = 3
beta_prior_means = matrix(0, nrow = 2, ncol = ncol(data_train))
beta_prior_sds =  matrix(1.0, nrow = 2, ncol = ncol(data_train))
beta_proposal_sds =  matrix(0.1, nrow = 2, ncol = ncol(data_train))
theta_init = matrix(0, nrow = n, ncol = horizon)
theta_init[,1] = rnorm(n)
theta_init[,1] = theta_init[,1]*sign(theta_init[,1]*theta[,1])
for (h in 2:horizon) {
  theta_init[,h] = theta_init[,1]
}

all_samples <- gpirtMCMC(data_train, SAMPLE_ITERS,BURNOUT_ITERS,
                     THIN=THIN, CHAIN=CHAIN, vote_codes = NULL,
                     beta_prior_means = beta_prior_means,
                     beta_prior_sds = beta_prior_sds, 
                     beta_proposal_sds = beta_proposal_sds,
                     theta_os = theta_os, theta_ls = theta_ls, 
                     theta_init = theta_init,
                     thresholds=NULL, SEED=SEED, constant_IRF = CONSTANT_IRF)

library(rstan)
sims <- matrix(rnorm((1+SAMPLE_ITERS)*CHAIN), nrow = 1+SAMPLE_ITERS, ncol = CHAIN)
theta_rhats = matrix(rnorm(n*horizon), nrow = n, ncol = horizon)
for(i in 1:n){
    for(h in 1:horizon){
        for(c in 1:CHAIN){
            pred_theta = colMeans(all_samples[[c]]$theta)
            sims[,c] = sign(cor(theta[,h],pred_theta[,h]))*all_samples$samples[[c]]$theta[,i,h]
        }
        theta_rhats[i, h] = Rhat(sims)
    }
}

samples = all_samples[[1]]
SAMPLE_ITERS = SAMPLE_ITERS/THIN

xs = seq(-5,5,0.01)
# pred_theta = colMeans(samples$theta)

pred_theta = matrix(0, nrow=nrow(data), ncol=dim(data)[3])
pred_theta_sd = matrix(0, nrow=nrow(data), ncol=dim(data)[3])
# 0 to keep, 1 to drop
# drop_wrong_signs = array(array(0, n*horizon*SAMPLE_ITERS), 
#                          c(n, horizon, SAMPLE_ITERS))
for(i in 1:n){
    for (h in 1:horizon) {
      tmp = samples$theta[-1,i,h]
      # drop wrong sign
      # drop_wrong_sign = (sign(cor(theta[,h],colMeans(samples$theta)[,h]))*tmp*theta[i,h]<0)
      # drop_wrong_signs[i,h,] = drop_wrong_sign
      # tmp = tmp[drop_wrong_sign==0]
      # if(is.na(mean(tmp))){
      #   tmp = samples$theta[-1,i,h]
      # }
      pred_theta[i,h] = mean(tmp)
      pred_theta_sd[i,h] = sd(tmp)
    }
}

pred_theta_ll = matrix(0, nrow=nrow(data), ncol=dim(data)[3])

for(i in 1:n){
    for (h in 1:horizon) {
        pred_theta_ll[i,h] = log(dnorm(theta[i,h], mean=pred_theta[i,h], sd=pred_theta_sd[i,h]))
    }
}

for(it in 1:SAMPLE_ITERS){
  for(h in 1:horizon){
    for(j in 1:m){
      samples$f[[it]][,j,h] = samples$f[[it]][,j,h] + samples$beta[[it]][1,j,h] + samples$beta[[it]][2,j,h]*samples$theta[it,,h]
      samples$fstar[[it]][,j,h] = samples$fstar[[it]][,j,h] + samples$beta[[it]][1,j,h] + samples$beta[[it]][2,j,h]*xs
    }
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

sample_IRFs = vector(mode = "list", length = SAMPLE_ITERS)
for (iter in 1:SAMPLE_ITERS){
    # recover fstar
    # sample_IRFs[[iter]] = recover_fstar(BURNOUT_ITERS+iter-1,samples$f[[iter]],data_train, 
    #               as.matrix(samples$theta[iter,,]), samples$threshold[iter,],
    #             beta_prior_means, beta_prior_sds,
    #               constant_IRF=CONSTANT_IRF)$fstar
    sample_IRFs[[iter]] = samples$fstar[[iter]]
}

for (i in 1:n) {
    for (j in 1:m) {
        for(h in 1:horizon){
            if(!is.na(data[[i,j,h]])){
                lls = matrix(0,nrow=SAMPLE_ITERS, ncol = C)
                y_pred = rep(0, SAMPLE_ITERS)
                pred_idx = 1+as.integer((pred_theta[i,h]+5)*100)
                for (iter in 1:SAMPLE_ITERS) {
                    f_pred = sample_IRFs[[iter]][pred_idx, j, h]
                    ll = ordinal_lls(f_pred, samples$threshold[iter,])
                    lls[iter,] = ll
                    y_pred[iter] =  which.max(ll)
                }
                # if(sum(drop_wrong_signs[i,h,]==0)>1){
                #   ll = log(colMeans(exp(lls[drop_wrong_signs[i,h,]==0,])))
                # }else{
                #   ll = log(exp(lls[drop_wrong_signs[i,h,]==0,]))
                # }
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
idx = (as.integer(min(theta)*100+500)):(as.integer(max(theta)*100+500))
idx = 401:601
# IRFs = Reduce("+",samples$IRFs)/length(samples$IRFs)
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
            IRFs[iter, ] = sample_IRFs[[iter]][idx, j, h]
        }
        probs = getprobs_gpirt(xs[idx], t(IRFs), 
                samples$threshold[1:SAMPLE_ITERS,])
        tmp = probs %>% 
            group_by(xs) %>%
            summarize(icc=sum(order*p))
        gpirt_iccs[,j,h] = tmp$icc
        cor_icc[j,h] = cor(gpirt_iccs[,j,h], true_iccs[,j,h])
        rmse_icc[j,h] = sqrt(mean((gpirt_iccs[,j,h]-true_iccs[,j,h])^2))
    }
}

cor_theta = c()
for (i in 1:horizon) {
    cor_theta = c(cor_theta, cor(theta[,i], pred_theta[,i]))
}

plot(samples$theta[10,,10],samples$f[[10]][,1,10])

plot(xs[idx],gpirt_iccs[,1,1], ylim=c(1,2))
mask = is.na(data_train[,1,1])
points(pred_theta[!mask,1], data_train[!mask,1,1])

tmp=rep(0,SAMPLE_ITERS)
for(it in 1:SAMPLE_ITERS){
  tmp[it] = samples$beta[[it]][2,1,1]
}
plot(1:SAMPLE_ITERS,tmp)

print("Average Rhat:")
print(mean(theta_rhats))
print(max(theta_rhats))
print(mean(abs(cor_theta)))
print(mean(pred_theta_sd))
print(mean(pred_theta_ll))
print(mean(train_lls[!is.infinite(train_lls)]))
print(mean(train_acc[!is.infinite(train_lls)]))
print(mean(pred_lls[!is.infinite(pred_lls)]))
print(mean(pred_acc[!is.infinite(pred_lls)]))
print(mean(array(abs(cor_icc), n*horizon)))
print(mean(array(rmse_icc, n*horizon)))

save(gpirt_iccs, true_iccs, theta, pred_theta,pred_theta_ll,pred_theta_sd,train_lls,
      train_acc, pred_lls, pred_acc,cor_icc, rmse_icc, theta_rhats,
      file=paste("./results/gpirt_", HYP, ".RData" , sep=""))
