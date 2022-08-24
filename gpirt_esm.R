library(gpirt)
library(dplyr)
library(ggplot2)
library(stats)

gpirt_path = "~/Documents/Github/OrdGPIRT"
setwd(gpirt_path)
SAMPLE_ITERS = 100
BURNOUT_ITERS = 100
TYPE = "CST"

source("2PL_esm.R")

# esm.BFI09, esm.BFI19 and esm.BFI04 
# data = read.csv("./data/esm_w1_redacted.csv")
# data = data[,c("esm.IDnum.w1","esm.PRO01.w1","esm.PRO03.w1","esm.PRO04.w1",
#                "esm.BFI09.w1","esm.BFI04.w1","esm.BFI19.w1")]
# colnames(data) = c("SID", "freq", "hourblock","day","N_relax", "N_depress", "N_worried")
# data = data[!is.na(data$N_relax),]
# data = data[!is.na(data$N_depress),]
# data = data[!is.na(data$N_worried),]
# irf_names = c("relax", "depress", "worried")
# 
# num_days = length(unique(data$day))
# num_hours = length(unique(data$hourblock))
# unique_ids = unique(data$SID)
# n = length(unique_ids)
# m = 3
# horizon = num_days*num_hours
# C = 5

# ems_data = array(array(NA, n*m*horizon), c(n, m, horizon))
# for(iter in 1:nrow(data)){
#   id = data$SID[iter]
#   hourblock = data$hourblock[iter]
#   day = data$day[iter]
#   N_relax = data$N_relax[iter]
#   N_depress = data$N_depress[iter]
#   N_worried = data$N_worried[iter]
#   i = which(id==unique_ids)
#   h = (day-1)*num_hours + hourblock
#   ems_data[i,1,h] = N_relax
#   ems_data[i,2,h] = N_depress
#   ems_data[i,3,h] = N_worried
# }

# data = ems_data

if(TYPE=="GP"){
  theta_os = 1
  theta_ls = as.integer(horizon)
}else if(TYPE=="CST"){
  theta_os = 0
  theta_ls = -1
}else{
  theta_os = 0
  theta_ls = as.integer(horizon)
}

theta_init = together_pred_theta

set.seed(12345)

SEED = 1
THIN = 1
CHAIN = 1
constant_IRF = 1
beta_prior_sds =  matrix(0.0, nrow = 2, ncol = ncol(data))
samples <- gpirtMCMC(data, SAMPLE_ITERS,BURNOUT_ITERS,
                     THIN, CHAIN, theta_init = theta_init,
                     beta_prior_sds = beta_prior_sds, theta_os = theta_os,
                     theta_ls = theta_ls, vote_codes = NULL, thresholds=NULL,
                     SEED=SEED, constant_IRF = constant_IRF)

samples = samples[[1]]
SAMPLE_ITERS = SAMPLE_ITERS/THIN

save.image(file='./results/gpirt_esm.RData')

# predicted ideology
pred_theta = matrix(0, nrow=nrow(data), ncol=dim(data)[3])
pred_theta_sd = matrix(0, nrow=nrow(data), ncol=dim(data)[3])
drop_wrong_signs = array(array(0, n*horizon*SAMPLE_ITERS), 
                         c(n, horizon, SAMPLE_ITERS))
for(h in 1:horizon){
  for (i in 1:n) {
    tmp = samples$theta[-1,i,h]
    drop_wrong_sign = (sign(cor(theta_init[,h],colMeans(samples$theta)[,h]))*tmp*theta_init[i,h]<0)
    drop_wrong_signs[i,h,] = drop_wrong_sign
    tmp = tmp[drop_wrong_sign==0]
    if(is.na(mean(tmp))){
      tmp = samples$theta[-1,i,h]
    }
    pred_theta[i,h] = mean(tmp)
    pred_theta_sd[i,h] = sd(tmp)
  }
}

sample_IRFs = vector(mode = "list", length = SAMPLE_ITERS)
for (iter in 1:SAMPLE_ITERS){
  # recover fstar
  # sample_IRFs[[iter]] = recover_fstar(BURNOUT_ITERS+iter-1,samples$f[[iter]],data,
  #                                     as.matrix(samples$theta[iter,,]), 
  #                                     samples$threshold[iter,],
  #                                     constant_IRF = constant_IRF)$fstar
  sample_IRFs[[iter]] = samples$fstar[[iter+1]]
}

# icc
xs = seq(-5,5,0.01)
source("getprob_gpirt.R")
idx = 301:701
gpirt_iccs = array(array(0, length(xs[idx])*m*horizon),
                   c(length(xs[idx]),m, horizon))

for (h in 1:horizon) {
  for (j in 1:m) {
    IRFs = matrix(0, nrow=SAMPLE_ITERS, ncol=length(idx))
    for(iter in 1:SAMPLE_ITERS){
      IRFs[iter, ] = sample_IRFs[[iter]][idx, j, h]
      IRFs[iter, ] = IRFs[iter, ] * sign(cor(IRFs[iter, ], grm_together_iccs[,j]))
    }
    probs = getprobs_gpirt(xs[idx], t(IRFs),
                           samples$threshold[1:SAMPLE_ITERS,])
    tmp = probs %>%
      group_by(xs) %>%
      summarize(icc=sum(order*p))
    gpirt_iccs[,j,h] = tmp$icc
  }
}

save.image(file='./results/gpirt_esm.RData')
