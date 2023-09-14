R_path="~/R/x86_64-redhat-linux-gnu-library/4.0"
.libPaths(R_path)
rstan_options(auto_write = TRUE)

args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

library(rstan)
library(gpirt)
library(dplyr)
library(ggplot2)
library(stats)

# gpirt_path = "~/Documents/Github/OrdGPIRT"
# setwd(gpirt_path)
load(file="./data/senate_data_90.RData")
gpirt_data = data

SAMPLE_ITERS = 500
BURNOUT_ITERS = 500
TYPE = "GP"
n = nrow(data)
m = ncol(data)
horizon = length(data[1,1,])
C = 2

if(TYPE=="GP"){
  theta_os = 1
  theta_ls = 7
}else if(TYPE=="CST"){
  theta_os = 0
  theta_ls = -1
}else{
  theta_os = 0
  theta_ls = 7
}


set.seed(12345)
theta_init = matrix(0, nrow = n, ncol = horizon)
NOMINATE_SCORE_DIM1 = matrix(NA, nrow = n, ncol = horizon)
NOMINATE_SCORE_DIM2 = matrix(NA, nrow = n, ncol = horizon)
for(h in 1:length(session_ids)){
  session_id = session_ids[h]
  members = read.csv(paste("./data/S", session_id, "_members.csv", sep=""))
  members = members[members$chamber=="Senate", c("icpsr", 
                                                 "nokken_poole_dim1", "nokken_poole_dim2")]
  members = members[members$icpsr!=40106, ]
  # excluse RUSSELL, Richard Brevard, Jr. of GA
  members = members[members$icpsr!=8138, ]
  # exclude James Danforth
  members = members[members$icpsr!=14447, ]
  current_unique_icpsrs = unique(members$icpsr)
  # nominate scores 
  nominate_scores = matrix(0, nrow=length(current_unique_icpsrs), ncol=2)
  idx = c()
  for(j in 1:length(current_unique_icpsrs)){
    icpsr = current_unique_icpsrs[j]
    nominate_scores[j,1] = members[members$icpsr==icpsr, "nokken_poole_dim1"]
    nominate_scores[j,2] = members[members$icpsr==icpsr, "nokken_poole_dim2"]
    idx = c(idx, which(icpsr==unique_icpsr))
  }
  theta_init[idx,h] = nominate_scores[,1] + 0.1*rnorm(n=length(idx))
  NOMINATE_SCORE_DIM1[idx,h] = nominate_scores[,1]
  NOMINATE_SCORE_DIM2[idx,h] = nominate_scores[,2]
}

SEED = 1
THIN = 4
CHAIN = 1

# code na as 0 for stan to ignore
gpirt_data_train = gpirt_data
gpirt_data_train[is.na(gpirt_data_train)] = 0

stan_data <- list(n=n,
                  m=m,
                  horizon=horizon,
                  K=C,
                  sigma=0.1,
                  y=gpirt_data_train)

# train stan model
fit <- stan(file = "doirt-synthetic.stan",
            data = stan_data, 
            warmup = BURNOUT_ITERS, 
            iter = BURNOUT_ITERS + SAMPLE_ITERS, 
            chains = CHAIN, 
            cores = 1, 
            thin = THIN,
            control=list(adapt_delta=.98, max_treedepth = 15),
            seed = SEED,
            refresh=1
)

fit_params <- as.data.frame(fit)

SAMPLE_ITERS = SAMPLE_ITERS / THIN

samples = list()
samples[["theta"]] = array(array(0, SAMPLE_ITERS*n*horizon), 
                           c(SAMPLE_ITERS,n, horizon))
xs = seq(-5,5,0.01)
idx = 1:length(xs)
for (it in 1:SAMPLE_ITERS) {
  samples[["fstar"]][[it]] =  array(array(0, length(xs[idx])*m*horizon), 
                                    c(length(xs[idx]),m, horizon))
  for(j in 1:m){
    for(h in 1:horizon){
      samples[["fstar"]][[it]][,j,h] = xs[idx]*fit_params[[paste("beta[",j,",",h,"]",sep="")]][it]
    }
  }
  for(i in 1:n){
    for(h in 1:horizon){
      samples[["theta"]][it,i,h] = fit_params[[paste("theta[",i,",",h,"]",sep="")]][it]
    }
  }
  
  samples[["threshold"]][[it]] =  array(array(0, m*(C+1)*horizon), 
                                        c(m, C+1, horizon))
  for(j in 1:m){
    for(h in 1:horizon){
      samples[["threshold"]][[it]][j,1,h] = -Inf
      samples[["threshold"]][[it]][j,C+1,h] = Inf
      for(c in 1:(C-1)){
        samples[["threshold"]][[it]][j,1+c,h] = fit_params[[paste("alpha[",j,",",h,",",c, "]",sep="")]][it]
      }
    }
  }
}

save.image(file='./results/doirt_senate_90.RData')

# summarize pred theta over iterations
pred_theta = matrix(0, nrow=n, ncol=horizon)
pred_theta_sd = matrix(0, nrow=n, ncol=horizon)

# predicted ideology
pred_theta = matrix(0, nrow=nrow(data), ncol=dim(data)[3])
pred_theta_sd = matrix(0, nrow=nrow(data), ncol=dim(data)[3])
for(h in 1:horizon){
  session_id = session_ids[h]
  members = read.csv(paste("./data/S", session_id, "_members.csv", sep=""))
  members = members[members$chamber=="Senate", c("icpsr", 
                                                 "nokken_poole_dim1", "nokken_poole_dim2")]
  # exlude BARKLEY, Dean
  members = members[members$icpsr!=40106, ]
  # excluse RUSSELL, Richard Brevard, Jr. of GA
  members = members[members$icpsr!=8138, ]
  # exclude James Danforth
  members = members[members$icpsr!=14447, ]
  current_unique_icpsrs = unique(members$icpsr)
  # nominate scores 
  nominate_scores = matrix(0, nrow=length(current_unique_icpsrs), ncol=2)
  idx = c()
  for(j in 1:length(current_unique_icpsrs)){
    icpsr = current_unique_icpsrs[j]
    # idx = which(icpsr == current_unique_icpsrs)
    nominate_scores[j,1] = members[members$icpsr==icpsr, "nokken_poole_dim1"]
    nominate_scores[j,2] = members[members$icpsr==icpsr, "nokken_poole_dim2"]
    idx = c(idx, which(icpsr==unique_icpsr))
  }
  mask = rep(0, SAMPLE_ITERS)
  for(iter in 1:SAMPLE_ITERS){
    if(cor(samples$theta[1,,h],  samples$theta[iter,,h])>0){
      mask[iter] = 1
    }
  }
  for (j in 1:length(current_unique_icpsrs)) {
    i = idx[j]
    if(sum(mask)>=sum(!mask)){
      tmp = samples$theta[mask==1,i,h]
    }
    else{
      tmp = samples$theta[mask==0,i,h]
    }
    pred_theta[i,h] = mean(tmp)
    pred_theta_sd[i,h] = sd(tmp)
  }
}

###################################
# effective sample size diagonis
library(coda)

ESS_theta = matrix(0, nrow = n, ncol=h)
for (h in 1:horizon) {
  for(i in 1:n){
    ESS_theta[i,h] = coda::effectiveSize(samples$theta[,i,h])
  }
}


###################################
# effective sample size diagonis
library(coda)

ESS_f = array(rep(0,n*m*horizon), c(n,m,horizon))
for(h in 1:horizon){
  for(j in 1:m){
    for(i in 1:n){
      tmp = rep(0, SAMPLE_ITERS)
      for(it in 1:SAMPLE_ITERS){
        tmp[it]= samples$f[[it]][i,j,h] 
      }
      ESS_f[i,j,h] = coda::effectiveSize(tmp)
    }
  }
}

# train/test statistics
# train
train_lls = c()
train_acc = c()
train_response = c()
train_prediction = c()

# gp IRFs
xs = seq(-5,5,0.01)
idx = 201:801
gpirt_iccs = array(array(0, length(xs[idx])*m*horizon),
                   c(length(xs[idx]),m, horizon))

source("getprob_gpirt.R")
for (h in 1:horizon){
  session_id = session_ids[h]
  
  rollcalls = read.csv(paste("./data/S", session_id, "_rollcalls.csv", sep=""))
  rollcalls = rollcalls[,c("congress", "rollnumber", "yea_count","nay_count","date" )]
  rollcalls = rollcalls[(rollcalls$yea_count!=0)&(rollcalls$nay_count!=0),]
  
  rollcall_ids = unique(rollcalls$rollnumber)
  for (j in 1:length(rollcall_ids)) {
    print(j)
    IRFs = matrix(0, nrow=SAMPLE_ITERS, ncol=length(xs))
    mask = rep(0, SAMPLE_ITERS)
    for(iter in 1:SAMPLE_ITERS){
      IRFs[iter, ] = samples$fstar[[iter]][, j, h]
      if (cor(samples$fstar[[1]][, j, h],samples$fstar[[iter]][, j, h])>0){
        mask[iter] = 1
      }
    }
    thresholds = matrix(0,nrow=SAMPLE_ITERS,ncol=C+1)
    for(iter in 1:SAMPLE_ITERS){
      thresholds[iter, ] = samples$threshold[[iter]][j,,h]
    }
    if(sum(mask)>=sum(!mask)){
      probs = getprobs_gpirt(xs, t(IRFs[mask==1,]), thresholds[mask==1,])
    }
    else{
      probs = getprobs_gpirt(xs, t(IRFs[mask==0,]), thresholds[mask==0,])
    }
    tmp = probs %>%
      group_by(xs) %>%
      summarize(icc=sum(order*p))
    gpirt_iccs[,j,h] = probs$p[probs$order==2]
    
    for (i in 1:n) {
      if(!is.na(gpirt_data[i,j,h])){
        # train
        pred_idx = 1+as.integer((pred_theta[i,h]+5)*100)
        ll = log(probs$p[probs$xs==xs[pred_idx]])
        y_pred = which.max(ll)
        train_acc = c(train_acc, y_pred==(gpirt_data[i,j, h]))
        train_lls = c(train_lls, ll[gpirt_data[i,j, h]])
        train_response = c(train_response, gpirt_data[i,j, h])
        train_prediction = c(train_prediction, y_pred)
      }
    }
  }
}

save.image(file='./results/doirt_senate_90.RData')