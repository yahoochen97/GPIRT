gpirt_path = "~/Documents/Github/gpirt"
setwd(gpirt_path)
library(Rcpp)
Rcpp::compileAttributes()
install.packages(gpirt_path, type="source", repos = NULL)#,lib=R_path, INSTALL_opts = '--no-lock')

library(gpirt)
library(dplyr)
library(ggplot2)
library(stats)
library(haven)

gpirt_path = "~/Documents/Github/OrdGPIRT"
setwd(gpirt_path)
SAMPLE_ITERS = 500
BURNOUT_ITERS = 500
TYPE = "GP"

# read data
data = read_dta("./data/SenatePeriods.dta")
congresses = sort(unique(data$congress))
horizon = length(congresses)

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
  # votes = read.csv(paste("./data/S", toString(congress), "_votes.csv", sep=""))
  rollcall_ids = unique(data$rollcall[data$congress==congress])
  senator_ids = unique(data$id[data$congress==congress])
  for(j in 1:length(rollcall_ids)){
    rollcall_id = rollcall_ids[j]
    for(k in 1:length(senator_ids)){
      i = which(senator_ids[k]==all_senator_ids)
      # cast_code = votes$cast_code[votes$rollnumber==rollcall_id & votes$icpsr==senator_ids[k]]
      tmp = data[data$congress==congress & data$id==senator_ids[k] & data$rollcall==rollcall_id, ]
      if(nrow(tmp)){
        rollcall_data[i,j,h] = tmp$dv + 1
      }
    }
  }
}


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

SAMPLE_ITERS = 500
BURNOUT_ITERS = 500
SEED = 1
THIN = 1
CHAIN = 1
beta_prior_sds =  matrix(0.5, nrow = 2, ncol = ncol(rollcall_data))
beta_proposal_sds =  matrix(0.1, nrow = 2, ncol = ncol(rollcall_data))
samples <- gpirtMCMC(rollcall_data, SAMPLE_ITERS,BURNOUT_ITERS,
                     THIN, CHAIN, beta_prior_sds = beta_prior_sds, 
                     beta_proposal_sds = beta_proposal_sds, theta_os = theta_os,
                     theta_ls = theta_ls, vote_codes = NULL, thresholds=NULL,
                     SEED=SEED, constant_IRF = 0)

samples = samples[[1]]
SAMPLE_ITERS = SAMPLE_ITERS/THIN

save.image(file='./results/gpirt_abortion.RData')

# predicted ideology
pred_theta = matrix(0, nrow=nrow(rollcall_data), ncol=dim(rollcall_data)[3])
pred_theta_sd = matrix(0, nrow=nrow(rollcall_data), ncol=dim(rollcall_data)[3])
for(h in 1:horizon){
  congress = congresses[h]
  current_unique_icpsrs = unique(data$id[data$congress==congress])
  idx = c()
  for(j in 1:length(current_unique_icpsrs)){
    icpsr = current_unique_icpsrs[j]
    idx = c(idx, which(icpsr==all_senator_ids))
  }
  for (j in 1:length(current_unique_icpsrs)) {
    i = idx[j]
    tmp = samples$theta[-1,i,h]
    # if(sd(tmp)<1){
    #   # drop wrong sign
    #   # drop_wrong_sign = (sign(cor(nominate_scores[,1],colMeans(samples$theta)[idx,h]))*tmp*nominate_scores[j,1]<0)
    #   # tmp = tmp[drop_wrong_sign==0]
    #   if(is.na(mean(tmp))){
    #     tmp = samples$theta[-1,i,h]
    #   }
    # }
    
    pred_theta[i,h] = mean(tmp)
    pred_theta_sd[i,h] = sd(tmp)
  }
}

xs = seq(-5,5,0.01)
idx = 301:701
gpirt_iccs = array(array(0, length(xs[idx])*m*horizon),
                   c(length(xs[idx]),m, horizon))

source("getprob_gpirt.R")
for (h in 1:horizon) {
  for (j in 1:m) {
    IRFs = matrix(0, nrow=SAMPLE_ITERS, ncol=length(idx))
    for(iter in 1:SAMPLE_ITERS){
      IRFs[iter, ] = samples$fstar[[iter]][idx, j, h]
      # IRFs[iter, ] = IRFs[iter, ] * sign(cor(IRFs[iter, ], IRFs[1, ]))
    }
    probs = getprobs_gpirt(xs[idx], t(IRFs),
                           samples$threshold)
    tmp = probs %>%
      group_by(xs) %>%
      summarize(icc=sum(order*p))
    gpirt_iccs[,j,h] = tmp$icc
  }
}

save.image(file='./results/gpirt_abortion.RData')