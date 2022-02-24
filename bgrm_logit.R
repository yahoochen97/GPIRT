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
HYP = paste(TYPE, "_C_", C, '_n_', n, '_m_', m, '_SEED_', SEED, sep="")
load(file=paste("./data/", HYP, ".RData" , sep=""))
set.seed(SEED)

library(rstan)

data_train[is.na(data_train)] = 0
stan_data <- list(N=n,
                  M=m,
                  C=C,
                  y=data_train)

# train stan model
fit <- stan(file = "bgrm_logit.stan",
            data = stan_data, 
            warmup = 1000, 
            iter = 3000, 
            chains = 1, 
            cores = 1, 
            thin = 4,
            control=list(adapt_delta=.98, max_treedepth = 15),
            seed = SEED,
            refresh=0
)

samples <- as.data.frame(fit)
SAMPLE_ITERS = length(samples[["theta[1]"]])
pred_theta = rep(0, n)
pred_lls = c()
pred_acc = c()
train_lls = c()
train_acc = c()
for (i in 1:nrow(data)) {
  pred_theta[i] = mean(samples[[paste("theta[",as.character(i), "]", sep="")]])
  for (j in 1:ncol(data)) {
    if(!is.na(data[[i,j]])){
      lls = matrix(0,nrow=SAMPLE_ITERS, ncol = C)
      ps = matrix(0, nrow=C, ncol=SAMPLE_ITERS)
      for(c in 1:C){
        ps[c,] = samples[[paste("p[",as.character(i), ",", as.character(j), ",", as.character(c) ,"]", sep="")]]
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

save(pred_theta,train_lls, train_acc, pred_lls, pred_acc,
     file=paste("./results/bgrm_", HYP, ".RData" , sep=""))
