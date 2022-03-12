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
library(ltm)
source("getprob_gpirt.R")
library(dplyr)
HYP = paste(TYPE, "_C_", C, '_n_', n, '_m_', m, '_SEED_', SEED, sep="")
load(file=paste("./data/", HYP, ".RData" , sep=""))
set.seed(SEED)

results = grm(data = data_train, na.action = NULL)
betas = results$coefficients
pred_theta = factor.scores(results, resp.patterns = data_train)
pred_theta = pred_theta$score.dat[,m+3]


pred_lls = c()
pred_acc = c()
train_lls = c()
train_acc = c()
xs = seq(-5,5,0.01)
idx = (as.integer(min(theta)*100+500)):(as.integer(max(theta)*100+500))
grm_iccs = matrix(0, nrow=length(idx), ncol=m)
true_iccs = matrix(0, nrow=nrow(grm_iccs), ncol=m)
cor_icc = rep(0, m)
rmse_icc = rep(0, m)

for (i in 1:nrow(data)) {
  for (j in 1:ncol(data)) {
    if(!is.na(data[[i,j]])){
      ps = rep(0, C+1)
      ps[C+1] = 1
      beta_coef = betas[[paste('Item ', j, sep='')]]
      for (c in 1:(C-1)){
        lp = beta_coef[[c]] - pred_theta[i]*beta_coef[[C]]
        ps[c+1] = 1 / (1+ exp(-lp))
      }
      
      ll = rep(0, C)
      for (c in 1:C) {
        ll[c] = log(ps[c+1] - ps[c])
      }
      
      y_pred = which.max(ll)
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

for (j in 1:m) {
  source("true_irf.R")
  probs = getprobs_gpirt(sign(cor(theta,pred_theta))*xs[idx], irfs, matrix(thresholds,nrow=1))
  tmp = probs %>% 
    group_by(xs) %>%
    summarize(icc=sum(order*p))
  true_iccs[,j] = tmp$icc
  tmp = c()
  for (i in 1:length(xs[idx])) {
    ps = rep(0, C+1)
    ps[C+1] = 1
    beta_coef = betas[[paste('Item ', j, sep='')]]
    for (c in 1:(C-1)){
      lp = beta_coef[[c]] - xs[idx][i]*beta_coef[[C]]
      ps[c+1] = 1 / (1+ exp(-lp))
    }
    
    p = rep(0, C)
    for (c in 1:C) {
      p[c] = ps[c+1] - ps[c]
    }
    tmp = c(tmp, sum(p*(1:C)))
  }
  grm_iccs[,j] = tmp
  cor_icc[j] = cor(grm_iccs[,j], true_iccs[,j])
  rmse_icc[j] = sqrt(mean((grm_iccs[,j]-true_iccs[,j])^2))
}

print(cor(theta,pred_theta))
print(mean(train_lls))
print(mean(train_acc))
print(mean(pred_lls))
print(mean(pred_acc))
print(mean(cor_icc))
print(mean(rmse_icc))

save(pred_theta,train_lls, train_acc, pred_lls, pred_acc, cor_icc, rmse_icc,
     file=paste("./results/grm_", HYP, ".RData" , sep=""))
