#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

if (length(args)==0) {
    SEED = 95
    C = 5
    n = 1000
    m = 50
    TYPE = "GP"
}
if (length(args)==5){
    SEED = as.integer(args[1])
    C = as.integer(args[2])
    n = as.integer(args[3])
    m = as.integer(args[4])
    TYPE = args[5]
}

# install gpirt package
R_path="~/R/x86_64-redhat-linux-gnu-library/4.0"
.libPaths(R_path)
# options("install.lock"=FALSE)
# gpirt_path = "../gpirt"
gpirt_path = "~/Documents/Github/gpirt"
setwd(gpirt_path)
library(Rcpp)
Rcpp::compileAttributes()
install.packages(gpirt_path, type="source", repos = NULL)#, lib=R_path)
setwd("../OrdGPIRT")
library(gpirt)
library(dplyr)

source("getprob_gpirt.R")
HYP = paste(TYPE, "_C_", C, '_n_', n, '_m_', m, '_SEED_', SEED, sep="")
load(file=paste("./data/", HYP, ".RData" , sep=""))
set.seed(SEED)

SAMPLE_ITERS = 100
BURNOUT_ITERS = 100
THIN = 1
beta_prior_sds =  matrix(0.5, nrow = 2, ncol = ncol(data_train))
samples <- gpirtMCMC(data_train, SAMPLE_ITERS,BURNOUT_ITERS, THIN,
                     beta_prior_sds = beta_prior_sds,
                     vote_codes = NULL, thresholds=NULL)
# save(samples, file = "vignettes/sdo.RData")
# load(file = "vignettes/sdo.RData")

THIN = 1
samples$theta = samples$theta[seq(1,SAMPLE_ITERS, THIN),]
samples$f = samples$f[,,seq(1,SAMPLE_ITERS, THIN)]
samples$threshold = samples$threshold[seq(1,SAMPLE_ITERS, THIN),]
samples$IRFs = samples$IRFs[,,seq(1,SAMPLE_ITERS, THIN)]

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
for (i in 1:nrow(data)) {
    for (j in 1:ncol(data)) {
        if(!is.na(data[[i,j]])){
          lls = matrix(0,nrow=SAMPLE_ITERS, ncol = C)
          y_pred = rep(0, SAMPLE_ITERS)
          for (iter in 1:SAMPLE_ITERS) {
            f_pred = samples$IRFs[as.integer((pred_theta[i]+5)*100), j, iter]
            ll = ordinal_lls(f_pred, samples$threshold[iter,])
            lls[iter,] = ll
            y_pred[iter] =  which.max(ll)
          }
          ll = colMeans(lls)
          y_pred = round(mean(y_pred))
            if(train_idx[i,j]==0){
                pred_acc = c(pred_acc, y_pred==(data[[i,j]]))
                pred_lls = c(pred_lls, ll[data[[i,j]]])
            }else{
                train_acc = c(train_acc, y_pred==(data[[i,j]]))
                train_lls = c(train_lls, ll[data[[i,j]]])
            }
        }
    }
}

# TODO: cor of icc
idx = (as.integer(min(theta)*100+500)):(as.integer(max(theta)*100+500))
gpirt_iccs = rowMeans(samples$IRFs[idx,,], dim=2)
true_iccs = matrix(0, nrow=nrow(gpirt_iccs), ncol=m)
cor_icc = rep(0, m)
rmse_icc = rep(0, m)
for (j in 1:m) {
    source("true_irf.R")
    probs = getprobs_gpirt(sign(cor(theta,pred_theta))*xs[idx], irfs, matrix(thresholds,nrow=1))
    tmp = probs %>% 
        group_by(xs) %>%
        summarize(icc=sum(order*p))
    true_iccs[,j] = tmp$icc
    probs = getprobs_gpirt(xs[idx], samples$IRFs[idx,j,], samples$threshold)
    tmp = probs %>% 
        group_by(xs) %>%
        summarize(icc=sum(order*p))
    gpirt_iccs[,j] = tmp$icc
    cor_icc[j] = cor(gpirt_iccs[,j], true_iccs[,j])
    rmse_icc[j] = sqrt(mean((gpirt_iccs[,j]-true_iccs[,j])^2))
}

print(cor(theta,pred_theta))
print(mean(train_lls[!is.infinite(train_lls)]))
print(mean(train_acc[!is.infinite(train_lls)]))
print(mean(pred_lls[!is.infinite(pred_lls)]))
print(mean(pred_acc[!is.infinite(pred_lls)]))
print(mean(cor_icc))
print(mean(rmse_icc))

save(samples, pred_theta,train_lls, train_acc, pred_lls, pred_acc,cor_icc, rmse_icc,
     file=paste("./results/gpirt_", HYP, ".RData" , sep=""))
