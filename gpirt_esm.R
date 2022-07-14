library(gpirt)
library(dplyr)
library(ggplot2)
library(stats)

gpirt_path = "~/Documents/Github/OrdGPIRT"
setwd(gpirt_path)
SAMPLE_ITERS = 100
BURNOUT_ITERS = 100
TYPE = "GP"

# esm.BFI21 and esm.BFI36
data = read.csv("./data/esm_w1_redacted.csv")
data = data[,c("esm.IDnum.w1","esm.PRO01.w1","esm.PRO03.w1","esm.PRO04.w1",
               "esm.BFI21.w1","esm.BFI36.w1")]
colnames(data) = c("SID", "freq", "hourblock","day","E_quiet", "E_outgoing")
data = data[!is.na(data$E_quiet),]
data = data[!is.na(data$E_outgoing),]

# data = data[(data$day<50) & (data$SID<5000), ]

num_days = length(unique(data$day))
num_hours = length(unique(data$hourblock))
unique_ids = unique(data$SID)
n = length(unique_ids)
m = 2 # two questions across all time periods
horizon = num_days*num_hours
C = 5

ems_data = array(array(NA, n*m*horizon), c(n,m, horizon))
for(iter in 1:nrow(data)){
  id = data$SID[iter]
  hourblock = data$hourblock[iter]
  day = data$day[iter]
  E_quiet = data$E_quiet[iter]
  E_outgoing = data$E_quiet[iter]
  i = which(id==unique_ids)
  h = (day-1)*num_hours + hourblock
  ems_data[i,1,h] = E_quiet
  ems_data[i,2,h] = E_outgoing
}

data = ems_data

if(TYPE=="GP"){
  theta_os = 1
  theta_ls = as.integer(horizon/2)
}else if(TYPE=="CST"){
  theta_os = 0
  theta_ls = -1
}else{
  theta_os = 0
  theta_ls = as.integer(horizon/2)
}

set.seed(12345)

SEED = 1
THIN = 1
CHAIN = 1
beta_prior_sds =  matrix(0.5, nrow = 2, ncol = ncol(data))
samples <- gpirtMCMC(data, SAMPLE_ITERS,BURNOUT_ITERS,
                     THIN, CHAIN,
                     beta_prior_sds = beta_prior_sds, theta_os = theta_os,
                     theta_ls = theta_ls, vote_codes = NULL, thresholds=NULL,
                     SEED=SEED)

samples = samples[[1]]
SAMPLE_ITERS = SAMPLE_ITERS/THIN

save.image(file='./results/gpirt_esm.RData')

# predicted ideology
pred_theta = matrix(0, nrow=nrow(data), ncol=dim(data)[3])
pred_theta_sd = matrix(0, nrow=nrow(data), ncol=dim(data)[3])
for(h in 1:horizon){
  for (i in 1:n) {
    tmp = samples$theta[-1,i,h]
    pred_theta[i,h] = mean(tmp)
    pred_theta_sd[i,h] = sd(tmp)
  }
}


sample_IRFs = vector(mode = "list", length = SAMPLE_ITERS)
for (iter in 1:SAMPLE_ITERS){
  # recover fstar
  sample_IRFs[[iter]] = recover_fstar(BURNOUT_ITERS+iter-1,samples$f[[iter]],data,
                                      as.matrix(samples$theta[iter,,]), samples$threshold[iter,])$fstar
}

# icc
xs = seq(-5,5,0.01)
source("getprob_gpirt.R")
idx = (as.integer(min(theta)*100+500)):(as.integer(max(theta)*100+500))
idx = 301:701
gpirt_iccs = array(array(0, length(xs[idx])*m*horizon),
                   c(length(xs[idx]),m, horizon))

for (h in 1:horizon) {
  for (j in 1:m) {
    IRFs = matrix(0, nrow=SAMPLE_ITERS, ncol=length(idx))
    for(iter in 1:SAMPLE_ITERS){
      IRFs[iter, ] = sample_IRFs[[iter]][idx, j, h]
    }
    probs = getprobs_gpirt(xs[idx], t(IRFs),
                           samples$threshold[1:SAMPLE_ITERS,])
    tmp = probs %>%
      group_by(xs) %>%
      summarize(icc=sum(order*p))
    gpirt_iccs[,j,h] = tmp$icc
  }
}
