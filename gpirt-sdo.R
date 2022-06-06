#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

if (length(args)==0) {
    SEED = 1
    C = 2
    n = 20
    m = 5
    horizon = 10
    TYPE = "GP"
}
if (length(args)==6){
    SEED = as.integer(args[1])
    C = as.integer(args[2])
    n = as.integer(args[3])
    m = as.integer(args[4])
    horizon = as.integer(args[5])
    TYPE = args[6]
}

# install gpirt package
R_path="~/R/x86_64-redhat-linux-gnu-library/4.0"
.libPaths(R_path)
# options("install.lock"=FALSE)
# gpirt_path = "../gpirt"
# gpirt_path = "~/Documents/Github/gpirt"
# setwd(gpirt_path)
# library(Rcpp)
# Rcpp::compileAttributes()
# install.packages(gpirt_path, type="source", repos = NULL,lib=R_path, INSTALL_opts = '--no-lock')
# setwd("../OrdGPIRT")
library(gpirt)
library(dplyr)

source("getprob_gpirt.R")
HYP = paste("GP_C_", C, '_n_', n, '_m_', m, '_h_', horizon, '_SEED_', SEED, sep="")
load(file=paste("./data/", HYP, ".RData" , sep=""))
HYP = paste(TYPE, "_C_", C, '_n_', n, '_m_', m, '_h_', horizon, '_SEED_', SEED, sep="")

SAMPLE_ITERS = 500
BURNOUT_ITERS = 500
if(TYPE=="GP"){
    theta_os = 1
    theta_ls = horizon
}else if(TYPE=="CST"){
    theta_os = 0
    theta_ls = -1
}else{
    theta_os = 0
    theta_ls = 1 + horizon / 3
}

THIN = 1
beta_prior_sds =  matrix(0.5, nrow = 2, ncol = ncol(data_train))
samples <- gpirtMCMC(data_train, SAMPLE_ITERS,BURNOUT_ITERS, THIN,
                     beta_prior_sds = beta_prior_sds, theta_os = theta_os,
                     theta_ls = theta_ls, vote_codes = NULL, thresholds=NULL)

SAMPLE_ITERS = SAMPLE_ITERS/THIN

xs = seq(-5,5,0.01)
pred_theta = colMeans(samples$theta)

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
    sample_IRFs[[iter]] = recover_fstar(BURNOUT_ITERS+iter-1,samples$f[[iter]],data_train, 
                  as.matrix(samples$theta[iter,,]), samples$threshold[iter,])$fstar
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
                ll = log(colMeans(exp(lls)))
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
idx = 301:701
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
            # IRFs[iter,] = samples$IRFs[[iter]][idx,j,h]
            # IRFs[iter,] = recover_fstar(BURNOUT_ITERS+iter-1,samples$f[[iter]],data_train, 
            #                        as.matrix(samples$theta[iter,,]), samples$threshold[iter,])$fstar[idx,j,h]
            IRFs[iter, ] = sample_IRFs[[iter]][idx, j, h]
        }
        probs = getprobs_gpirt(xs[idx], t(IRFs), 
                               samples$threshold)
        tmp = probs %>% 
            group_by(xs) %>%
            summarize(icc=sum(order*p))
        gpirt_iccs[,j,h] = tmp$icc
        cor_icc[j,h] = cor(gpirt_iccs[,j,h], true_iccs[,j,h])
        rmse_icc[j,h] = sqrt(mean((gpirt_iccs[,j,h]-true_iccs[,j,h])^2))
    }
}

cor_theta = c()
for (i in 1:h) {
    cor_theta = c(cor_theta, cor(theta[,i], pred_theta[,i]))
}

print(mean(abs(cor_theta)))
print(mean(train_lls[!is.infinite(train_lls)]))
print(mean(train_acc[!is.infinite(train_lls)]))
print(mean(pred_lls[!is.infinite(pred_lls)]))
print(mean(pred_acc[!is.infinite(pred_lls)]))
print(mean(array(abs(cor_icc), n*horizon)))
print(mean(array(rmse_icc, n*horizon)))

save(gpirt_iccs, true_iccs, pred_theta,train_lls, train_acc, pred_lls, pred_acc,cor_icc, rmse_icc,
     file=paste("./results/gpirt_", HYP, ".RData" , sep=""))
