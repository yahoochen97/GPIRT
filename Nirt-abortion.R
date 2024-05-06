library(gpirt)
library(dplyr)
library(ggplot2)
library(stats)
library(haven)

gpirt_path = "~/Documents/Github/OrdGPIRT"
setwd(gpirt_path)
TYPE = "NIRT"

# read data
data = read_dta("./data/SenatePeriods.dta")
congresses = sort(unique(data$congress))
horizon = length(congresses)

# remove president
data = data[data$name!="REAGAN",]
data = data[data$name!="BUSH",]
data = data[data$name!="CLINTON",]

m = 1 # number of vote calls
all_senator_ids = unique(data$id) # id of all senators in this dataset
n = length(all_senator_ids) # number of senators
C = 2 # binary questions

# iterate each congress
for (congress in congresses) {
  num_rollcalls = nrow(unique(data[data$congress==congress, c("rollcall")]))
  print(num_rollcalls)
  if(m < num_rollcalls){
    m = num_rollcalls
  }
}

# n by m by horizon
rollcall_data =  array(array(NA, n*m*horizon), c(n, m, horizon))
for (h in 1:length(congresses)) {
  congress = congresses[h]
  rollcall_ids = unique(data[data$congress==congress, "rollcall"]$rollcall)
  senator_ids = unique(data[data$congress==congress,"id"]$id)
  for(j in 1:length(rollcall_ids)){
    rollcall_id = rollcall_ids[j]
    for(k in 1:length(senator_ids)){
      i = which(senator_ids[k]==all_senator_ids)
      tmp = data[data$congress==congress & data$id==senator_ids[k] & data$rollcall==rollcall_id, ]
      if(nrow(tmp)){
        rollcall_data[i,j,h] = tmp$dv + 1
      }
    }
  }
}

# n by horizon
theta_init = matrix(0, nrow = n, ncol = horizon)
nominate_theta =  array(array(NA, n*horizon), c(n, horizon))
nominate_theta2 =  array(array(NA, n*horizon), c(n, horizon))
for (h in 1:length(congresses)) {
  congress = congresses[h]
  members = read.csv(paste("./data/S", congress, "_members.csv", sep=""))
  members = members[members$chamber=="Senate", c("icpsr", "nokken_poole_dim1", "nokken_poole_dim2")]
  senator_ids = unique(data[data$congress==congress,"id"]$id)
  # nominate scores 
  nominate_scores = matrix(0, nrow=length(senator_ids), ncol=2)
  idx = c()
  for(j in 1:length(senator_ids)){
    icpsr = senator_ids[j]
    nominate_scores[j,1] = members[members$icpsr==icpsr, "nokken_poole_dim1"]
    nominate_scores[j,2] = members[members$icpsr==icpsr, "nokken_poole_dim2"]
    idx = c(idx, which(icpsr==all_senator_ids))
  }
  nominate_theta[idx,h] = nominate_scores[,1]
  nominate_theta2[idx,h] = nominate_scores[,2]
  theta_init[idx,h] = nominate_scores[,1] + 0.1*rnorm(length(idx))
}


theta_os = 0
theta_ls = 7

SAMPLE_ITERS = 100
BURNOUT_ITERS = 100
SEED = 1
THIN = 1
CHAIN = 1
beta_prior_sds =  matrix(1.0, nrow = 3, ncol = ncol(rollcall_data))
samples <- gpirtMCMC(rollcall_data, SAMPLE_ITERS,BURNOUT_ITERS,
                     THIN, CHAIN, beta_prior_sds = beta_prior_sds,
                     theta_init = theta_init, theta_os = theta_os,
                     theta_ls = theta_ls, vote_codes = NULL, thresholds=NULL,
                     SEED=SEED, constant_IRF = 0)

samples = samples[[1]]
SAMPLE_ITERS = SAMPLE_ITERS/THIN

# predicted ideology
pred_theta = matrix(0, nrow=nrow(rollcall_data), ncol=dim(rollcall_data)[3])
pred_theta_sd = matrix(0, nrow=nrow(rollcall_data), ncol=dim(rollcall_data)[3])
for(h in 1:horizon){
  congress = congresses[h]
  senator_ids = unique(data[data$congress==congress,"id"]$id)
  idx = c()
  for(j in 1:length(senator_ids)){
    icpsr = senator_ids[j]
    idx = c(idx, which(icpsr==all_senator_ids))
  }
  for (j in 1:length(senator_ids)) {
    i = idx[j]
    tmp = samples$theta[-1,i,h]
    pred_theta[i,h] = mean(tmp)
    pred_theta_sd[i,h] = sd(tmp)
  }
}

xs = seq(-5,5,0.01)
for(it in 1:SAMPLE_ITERS){
  for(h in 1:horizon){
    for(j in 1:m){
      samples$f[[it]][,j,h] = samples$f[[it]][,j,h] + samples$beta[[it]][1,j,h] + samples$beta[[it]][2,j,h]*samples$theta[it,,h]
      samples$fstar[[it]][,j,h] = samples$fstar[[it]][,j,h] + samples$beta[[it]][1,j,h] + samples$beta[[it]][2,j,h]*xs
    }
  }
}

save.image(file='./results/Nirt_abortion.RData')

source("getprob_gpirt.R")

# NIRT
train_lls = c()
train_acc = c()
response1 = c()
prediction1 = c()

# DO-IRT
train_lls2 = c()
train_acc2 = c()
response2 = c()
prediction2 = c()

for(h in 1:horizon) {
  congress = congresses[h]
  # votes = read.csv(paste("./data/S", toString(congress), "_votes.csv", sep=""))
  rollcall_ids = unique(data[data$congress==congress, "rollcall"]$rollcall)
  for (j in 1:length(rollcall_ids)){
    IRFs = matrix(0, nrow=SAMPLE_ITERS, ncol=length(xs))
    for(iter in 1:SAMPLE_ITERS){
       IRFs[iter, ] = samples$fstar[[iter]][, j, h]# *tmp
    }
    thresholds = matrix(0,nrow=SAMPLE_ITERS,ncol=C+1)
    for(iter in 1:SAMPLE_ITERS){
      thresholds[iter, ] = samples$threshold[[iter]][j,,h]
    }
    probs = getprobs_gpirt(xs, t(IRFs), thresholds)
    for (i in 1:n) {
      if(!is.na(rollcall_data[[i,j,h]])  & !is.na(nominate_theta[i,h]) ){
        # GP-IRT
        pred_idx = 1+as.integer((pred_theta[i,h]+5)*100)
        ll = log(probs$p[probs$xs==xs[pred_idx]])
        y_pred = which.max(ll)
        train_acc = c(train_acc, y_pred==(rollcall_data[[i,j, h]]))
        train_lls = c(train_lls, ll[rollcall_data[[i,j, h]]])
        response1 = c(response1, rollcall_data[[i,j, h]])
        prediction1 = c(prediction1, y_pred)
        
        # NOMINATE
        pred_idx = 1+as.integer((nominate_theta[i,h]+5)*100)
        ll = log(probs$p[probs$xs==xs[pred_idx]])
        y_pred = which.max(ll)
        train_acc2 = c(train_acc2, y_pred==(rollcall_data[[i,j, h]]))
        train_lls2 = c(train_lls2, ll[rollcall_data[[i,j, h]]])
        response2 = c(response2, rollcall_data[[i,j, h]])
        prediction2 = c(prediction2, y_pred)
      }
    }
  }
}

print(mean(train_acc))
print(mean(train_acc2))

print(t.test(train_acc,train_acc2,paired=TRUE))

print(mean(train_lls))
print(mean(train_lls2))
print(t.test(train_lls,train_lls2,paired=TRUE))

library(pROC)

print(auc(response1, prediction1))
print(auc(response2, prediction2))
