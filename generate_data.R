#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

# gpirt_path = "~/Documents/Github/gpirt"
# setwd(gpirt_path)
# setwd("../OrdGPIRT")

if (length(args)==0) {
  SEED = 15
  C = 2
  n = 100
  m = 50
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

R_path="~/R/x86_64-redhat-linux-gnu-library/4.0"
.libPaths(R_path)
SIGMA = 1
source("getprob_gpirt.R")
HYP = paste(TYPE, "_C_", C, '_n_', n, '_m_', m, '_h_', horizon, '_SEED_', SEED, sep="")

set.seed(SEED)
thresholds = rep(0, C+1)
thresholds[1] = -Inf
thresholds[C+1] = Inf
if(C>2){
  thresholds[2:C] = seq(-2,2,length.out=C-1)
  thresholds[2:C] = thresholds[2:C] + rnorm(C-1,0,0.25)
  thresholds[2:C] = thresholds[2:C] - mean(thresholds[2:C])
  thresholds[2:C] = (thresholds[2:C])/sd(thresholds[2:C])
  thresholds[2:C] = sort(thresholds[2:C])
}else{
  thresholds[2] = rnorm(1,0,0.25)
}

if(TYPE=="2PL"){
  gen_responses <- function(theta, alpha, beta, thresholds) {
    # ordinal regression
    C <- length(thresholds) - 1
    n <- length(theta)
    m <- length(alpha)
    responses <- matrix(0, n, m)
    for ( j in 1:m ) {
      for ( i in 1:n ) {
        f = alpha[j] + beta[j] * theta[i]
        ps = rep(0, C)
        for (c in 1:C) {
          z1 = thresholds[c] - f
          z2 = thresholds[c+1] -f
          ps[c] = pnorm(z2) - pnorm(z1)
        }
        responses[i, j] <- sample(1:C, 1, prob = ps)
      }
    }
    return(responses)
  }
  
  theta <- runif(n, -3, 3) # Respondent ability parameters
  alpha <- runif(m, -2, 2) # Item difficulty parameters
  beta  <- runif(m, -2, 2) # Item discrimination parameters
  data <- gen_responses(theta, alpha, beta, thresholds)
}
if(TYPE=="GP"){
  library(MASS)
  gen_responses <- function(theta, anchor_xs, anchor_ys, thresholds) {
    # ordinal regression
    C <- length(thresholds) - 1
    n <- nrow(theta)
    horizon <- ncol(theta)
    m <- nrow(anchor_xs)
    NUM_ANCHOR = ncol(anchor_xs)
    responses <- array(rep(0, n*m*horizon), c(n, m, horizon))
    for (h in 1:horizon) {
      for ( j in 1:m ) {
        K = SEKernel(anchor_xs[j,,h], sigma=SIGMA)
        K = K + diag(1e-6, NUM_ANCHOR,NUM_ANCHOR)
        inv_K = ginv(K)
        for ( i in 1:n ) {
          K1 = dnorm(theta[i,h]-anchor_xs[j,,h], sd=SIGMA)/dnorm(0, sd=SIGMA)
          f = K1 %*% inv_K %*% (anchor_ys[j,,h])
          ps = rep(0, C)
          for (c in 1:C) {
            z1 = thresholds[c] - f
            z2 = thresholds[c+1] -f
            ps[c] = pnorm(z2) - pnorm(z1)
          }
          responses[i,j,h] <- sample(1:C, 1, prob = ps)
        }
      }
    }
    return(responses)
  }
  theta = matrix(0, nrow = n, ncol = horizon)
  if(horizon>1){
    for (i in 1:n) {
      K = SEKernel(1:horizon, sigma=as.integer(1+horizon/2))
      K = K + diag(1e-6, horizon, horizon)
      theta[i,] <- t(chol(K))%*%rnorm(horizon)  # Respondent ability parameters
    }
  }else{
    theta[,1] = rnorm(n)
  }
  
  xs = seq(-5,5,0.01)
  NUM_ANCHOR = 50
  anchor_xs <- array(rep(0,m*NUM_ANCHOR,horizon), c(m, NUM_ANCHOR, horizon))
  anchor_ys <- array(rep(0,m*NUM_ANCHOR,horizon), c(m, NUM_ANCHOR, horizon))
  for (h in 1:horizon){
    idx = (as.integer(min(theta)*100+500)):(as.integer(max(theta)*100+500))
    idx = 301:701
    for (j in 1:m) {
      anchor_xs[j,,h] = seq(-2,2, length.out = NUM_ANCHOR) # anchor points
      K = SEKernel(anchor_xs[j,,h], sigma=SIGMA)
      K = K + diag(1e-6, NUM_ANCHOR,NUM_ANCHOR)
      anchor_ys[j,,h]  <- t(chol(K))%*%rnorm(NUM_ANCHOR) 
      anchor_ys[j,,h] = anchor_ys[j,,h] - mean(anchor_ys[j,,h])
      anchor_ys[j,,h] = anchor_ys[j,,h] / sd(anchor_ys[j,,h])
    }
  }
  data <- gen_responses(theta, anchor_xs, anchor_ys, thresholds)
}

# split into train and test data
N = n*m*horizon
na_mask = matrix(is.na(data), nrow = N)
train_idx = rep(0,N)
train_idx[sample((1:N)[!na_mask],as.integer(0.8*N),replace=FALSE)] = 1
train_idx = array(train_idx, c(n,m,horizon))
data_train = data

for (i in 1:nrow(data)) {
  for (j in 1:ncol(data)) {
    for (h in 1:horizon) {
      if(train_idx[i,j,h]==0){
        data_train[i,j,h] = NA
      }
    }
  }
}

if(TYPE=="2PL"){
  save(data,data_train, train_idx,thresholds,theta,alpha,beta, 
       file=paste("./data/", HYP, ".RData" , sep=""))
}
if(TYPE=="GP"){
  save(data,data_train, train_idx,thresholds,theta,NUM_ANCHOR, anchor_xs,anchor_ys,SIGMA,
       file=paste("./data/", HYP, ".RData" , sep=""))
}
