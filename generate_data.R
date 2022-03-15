#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

if (length(args)==0) {
  SEED = 1
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

R_path="~/R/x86_64-redhat-linux-gnu-library/4.0"
.libPaths(R_path)
SIGMA = 1
source("getprob_gpirt.R")
HYP = paste(TYPE, "_C_", C, '_n_', n, '_m_', m, '_SEED_', SEED, sep="")

set.seed(SEED)
thresholds = rep(0, C+1)
thresholds[1] = -Inf
thresholds[C+1] = Inf
thresholds[2:C] = seq(-2,2,length.out=C-1)
thresholds[2:C] = thresholds[2:C] + rnorm(C-1,0,0.25)
thresholds[2:C] = thresholds[2:C] - mean(thresholds[2:C])
thresholds[2:C] = (thresholds[2:C])/sd(thresholds[2:C])
thresholds[2:C] = sort(thresholds[2:C])

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
  library(KRLS)
  library(MASS)
  gen_responses <- function(theta, anchor_xs, anchor_ys, thresholds) {
    # ordinal regression
    C <- length(thresholds) - 1
    n <- length(theta)
    m <- nrow(anchor_xs)
    NUM_ANCHOR = ncol(anchor_xs)
    responses <- matrix(0, n, m)
    for ( j in 1:m ) {
      K = gausskernel(anchor_xs[j,], sigma=SIGMA)
      K = K + diag(1e-6, NUM_ANCHOR,NUM_ANCHOR)
      inv_K = ginv(K)
      for ( i in 1:n ) {
        K1 = dnorm(theta[i]-anchor_xs[j,], sd=SIGMA)/dnorm(0, sd=SIGMA)
        f = K1 %*% inv_K %*% (anchor_ys[j,])
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
  theta <- seq(-1,1, length.out = n) # Respondent ability parameters
  xs = seq(-5,5,0.01)
  idx = (as.integer(min(theta)*100+500)):(as.integer(max(theta)*100+500))
  NUM_ANCHOR = 50
  anchor_xs <- matrix(0, nrow=m,ncol=NUM_ANCHOR)
  anchor_ys <- matrix(0, nrow=m,ncol=NUM_ANCHOR)
  for (j in 1:m) {
    anchor_xs[j,] = seq(-2,2, length.out = NUM_ANCHOR) # anchor points
    K = gausskernel(anchor_xs[j,], sigma=SIGMA)
    K = K + diag(1e-6, NUM_ANCHOR,NUM_ANCHOR)
    anchor_ys[j,]  <- t(chol(K))%*%rnorm(NUM_ANCHOR) 
    anchor_ys[j,] = anchor_ys[j,] - mean(anchor_ys[j,])
    anchor_ys[j,] = anchor_ys[j,] / sd(anchor_ys[j,])
  }
  data <- gen_responses(theta, anchor_xs, anchor_ys, thresholds)
}

# split into train and test data
N = nrow(data)*ncol(data)
na_mask = matrix(is.na(data), nrow = N)
train_idx = rep(0,N)
train_idx[sample((1:N)[!na_mask],as.integer(0.8*N),replace=FALSE)] = 1
train_idx = matrix(train_idx, nrow = nrow(data))
data_train = data

for (i in 1:nrow(data)) {
  for (j in 1:ncol(data)) {
    if(train_idx[i,j]==0){
      data_train[i,j] = NA
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
