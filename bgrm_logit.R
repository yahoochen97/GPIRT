#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

if (length(args)==0) {
  SEED = 1
  C = 5
  n = 100
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
HYP = paste(TYPE, "_C_", C, '_n_', n, '_m_', m, '_SEED_', SEED, sep="")
load(file=paste("./data/", HYP, ".RData" , sep=""))
set.seed(SEED)

library(dplyr)
source("getprob_gpirt.R")
library(rstan)
rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores())
xs = seq(-3,3,0.01)
data_train[is.na(data_train)] = 0
stan_data <- list(N=n,
                  M=m,
                  C=C,
                  Nstar=length(xs),
                  xstar=xs,
                  y=data_train)

# train stan model
fit <- stan(file = "bgrm_logit.stan",
            data = stan_data, 
            warmup = 500, 
            iter = 1000, 
            chains = 1, 
            cores = 1, 
            thin = 4,
            control=list(adapt_delta=.98, max_treedepth = 15),
            seed = SEED,
            refresh= 1
)

# saveRDS(fit, paste("./results/bgrm_", HYP, ".rds" , sep=""))

samples <- as.data.frame(fit)
SAMPLE_ITERS = length(samples[["theta[1]"]])
pred_theta = rep(0, n)
pred_lls = c()
pred_acc = c()
train_lls = c()
train_acc = c()
for (i in 1:nrow(data)) {
  pred_theta[i] = mean(samples[[paste("theta[",as.character(i), "]", sep="")]])
  i_star = as.integer(pred_theta[i]*100+300)
  for (j in 1:ncol(data)) {
    if(!is.na(data[[i,j]])){
      lls = matrix(0,nrow=SAMPLE_ITERS, ncol = C)
      ps = matrix(0, nrow=C, ncol=SAMPLE_ITERS)
      for(c in 1:C){
        ps[c,] = samples[[paste("p[",as.character(i_star), ",", as.character(j), ",", as.character(c) ,"]", sep="")]]
      }
      ll = log(rowMeans(ps))
      y_pred =  which.max(ll)
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

# TODO: cor of icc
idx = (as.integer(min(theta)*100+300)):(as.integer(max(theta)*100+300))
bgrm_iccs = matrix(0, nrow=length(idx), ncol=m)
true_iccs = matrix(0, nrow=nrow(bgrm_iccs), ncol=m)
cor_icc = rep(0, m)
rmse_icc = rep(0, m)
for (j in 1:m) {
  source("true_irf.R")
  probs = getprobs_gpirt(sign(cor(theta,pred_theta))*xs[idx], irfs, matrix(thresholds,nrow=1))
  tmp = probs %>% 
    group_by(xs) %>%
    summarize(icc=sum(order*p))
  true_iccs[,j] = tmp$icc
  tmp = c()
  for (i in 1:length(xs[idx])) {
    tmp = c(tmp, mean(samples[[paste("icc[", as.integer(xs[idx][i]*100+300),",",j,"]", sep="")]]))
  }
  bgrm_iccs[,j] = tmp
  cor_icc[j] = cor(bgrm_iccs[,j], true_iccs[,j])
  rmse_icc[j] = sqrt(mean((bgrm_iccs[,j]-true_iccs[,j])^2))
}

print(mean(cor_icc))
print(mean(rmse_icc))

save(pred_theta,train_lls, train_acc, pred_lls, pred_acc,cor_icc, rmse_icc,
     file=paste("./results/bgrm_", HYP, ".RData" , sep=""))
