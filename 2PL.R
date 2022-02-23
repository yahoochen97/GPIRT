#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

if (length(args)==0) {
  SEED = 1
  C = 5
  n = 1000
  m = 10
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
HYP = paste(TYPE, "_C_", C, '_n_', n, '_m_', m, '_SEED_', SEED, sep="")
load(file=paste("./data/", HYP, ".RData" , sep=""))
set.seed(SEED)

# data = data.matrix(SDO)[1:500,]
# unique_ys = unique(as.vector(data))
# C = length(unique(unique_ys[!is.na(unique_ys)]))

results = grm(data = data_train, na.action = NULL)
betas = results$coefficients
pred_theta = factor.scores(results, resp.patterns = data_train)
pred_theta = pred_theta$score.dat[,m+3]


pred_lls = c()
pred_acc = c()
train_lls = c()
train_acc = c()
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

print(cor(theta,pred_theta))
print(mean(train_lls))
print(mean(train_acc))
print(mean(pred_lls))
print(mean(pred_acc))
