#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

if (length(args)==0) {
    SEED = 1
    C = 2
    n = 100
    m = 20
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
gpirt_path = "~/Documents/Github/gpirt"
setwd(gpirt_path)
library(Rcpp)
Rcpp::compileAttributes()
install.packages(gpirt_path, type="source", repos = NULL)#,lib=R_path, INSTALL_opts = '--no-lock')
setwd("../OrdGPIRT")
library(gpirt)
library(dplyr)

source("getprob_gpirt.R")
HYP = paste(TYPE, "_C_", C, '_n_', n, '_m_', m, '_h_', horizon, '_SEED_', SEED, sep="")
load(file=paste("./data/", HYP, ".RData" , sep=""))
set.seed(SEED)

SAMPLE_ITERS = 100
BURNOUT_ITERS = 100
THIN = 1
beta_prior_sds =  matrix(0.5, nrow = 2, ncol = ncol(data_train))
samples <- gpirtMCMC(data_train, SAMPLE_ITERS,BURNOUT_ITERS, THIN,
                     beta_prior_sds = beta_prior_sds,
                     vote_codes = NULL, thresholds=NULL, SEED=SEED)

# THIN = 1
# samples$theta = samples$theta[seq(1,SAMPLE_ITERS, THIN),]
# samples$f = samples$f[,,seq(1,SAMPLE_ITERS, THIN)]
# samples$threshold = samples$threshold[seq(1,SAMPLE_ITERS, THIN),]
# samples$IRFs = samples$IRFs[,,seq(1,SAMPLE_ITERS, THIN)]

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
for (i in 1:n) {
    for (j in 1:m) {
        for(h in 1:horizon){
            if(!is.na(data[[i,j,h]])){
              lls = matrix(0,nrow=SAMPLE_ITERS, ncol = C)
              y_pred = rep(0, SAMPLE_ITERS)
              for (iter in 1:SAMPLE_ITERS) {
                f_pred = samples$IRFs[[iter]][as.integer((pred_theta[i,h]+5)*100), j, h]
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

# TODO: cor of icc
idx = (as.integer(min(theta)*100+500)):(as.integer(max(theta)*100+500))
idx = 301:701
# IRFs = Reduce("+",samples$IRFs)/length(samples$IRFs)
# gpirt_iccs = colMeans(IRFs, dim=1)
gpirt_iccs = array(array(0, length(xs[idx])*m*horizon), 
                   c(length(xs[idx]),m, horizon))
cor_icc = matrix(0, nrow=m, ncol=horizon)
rmse_icc = matrix(0, nrow=m, ncol=horizon)

for (h in 1:horizon) {
    
    true_iccs = matrix(0, nrow=nrow(gpirt_iccs), ncol=m)
    for (j in 1:m) {
        source("true_irf.R")
        # sign(cor(theta[,h],pred_theta[,h]))
        probs = getprobs_gpirt(sign(cor(theta[,h],pred_theta[,h]))*xs[idx], 
                               irfs, matrix(thresholds,nrow=1))
        tmp = probs %>% 
            group_by(xs) %>%
            summarize(icc=sum(order*p))
        true_iccs[,j] = tmp$icc
        IRFs = matrix(0, nrow=SAMPLE_ITERS, ncol=length(idx))
        for(iter in 1:SAMPLE_ITERS){
            IRFs[iter,] = samples$IRFs[[iter]][idx,j,h]
        }
        probs = getprobs_gpirt(xs[idx], t(IRFs), 
                               samples$threshold)
        tmp = probs %>% 
            group_by(xs) %>%
            summarize(icc=sum(order*p))
        gpirt_iccs[,j,h] = tmp$icc
        cor_icc[j,h] = cor(gpirt_iccs[,j,h], true_iccs[,j])
        rmse_icc[j,h] = sqrt(mean((gpirt_iccs[,j,h]-true_iccs[,j])^2))
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

save(gpirt_iccs, pred_theta,train_lls, train_acc, pred_lls, pred_acc,cor_icc, rmse_icc,
     file=paste("./results/gpirt_", HYP, ".RData" , sep=""))

# plot(1:10,samples$theta[1,1,],lty=1,ylim=c(-6,6))
# for(i in 2:100){
#     lines(1:10,samples$theta[i,1,])
# }

