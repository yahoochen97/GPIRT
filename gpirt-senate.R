library(gpirt)
library(dplyr)
library(ggplot2)
library(stats)

gpirt_path = "~/Documents/Github/OrdGPIRT"
setwd(gpirt_path)
load(file="./data/senate_data_92.RData")
SAMPLE_ITERS = 100
BURNOUT_ITERS = 100
TYPE = "GP"
n = nrow(data)
m = ncol(data)
horizon = length(data[1,1,])
C = 2

print(session_ids)

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

# fix_theta_flag = matrix(0, nrow=n, ncol=horizon)
# fix_theta_value = matrix(0, nrow=n, ncol=horizon)

# fix WARREN, Elizabeth (41301, 113th) to -1
# fix_theta_flag[which(unique_icpsr==41301), 7] = 1
# fix_theta_value[which(unique_icpsr==41301), 7] = -1.0

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
    # idx = which(icpsr == current_unique_icpsrs)
    nominate_scores[j,1] = members[members$icpsr==icpsr, "nokken_poole_dim1"]
    nominate_scores[j,2] = members[members$icpsr==icpsr, "nokken_poole_dim2"]
    idx = c(idx, which(icpsr==unique_icpsr))
  }
  # theta_init[idx,h] = theta_init[idx,h]*(2*sign(theta_init[idx,h]*nominate_scores[,1]<0)-1)
  # nominate_scores[,1] + 0.1*rnorm(n=length(idx))
  theta_init[idx,h] = nominate_scores[,1] + 0.1*rnorm(n=length(idx))
  NOMINATE_SCORE_DIM1[idx,h] = nominate_scores[,1]
  NOMINATE_SCORE_DIM2[idx,h] = nominate_scores[,2]
}

write.csv(NOMINATE_SCORE_DIM1, file="./data/NOMINATE1_theta.csv",row.names = FALSE)
write.csv(NOMINATE_SCORE_DIM2, file="./data/NOMINATE2_theta.csv",row.names = FALSE)

SEED = 1
THIN = 1
CHAIN = 1
beta_prior_sds =  matrix(0.5, nrow = 2, ncol = ncol(data))
beta_proposal_sds =  matrix(0.1, nrow = 2, ncol = ncol(data))
samples <- gpirtMCMC(data, SAMPLE_ITERS,BURNOUT_ITERS,
                     THIN, CHAIN, theta_init = theta_init,
                     beta_proposal_sds = beta_proposal_sds,
                     beta_prior_sds = beta_prior_sds, theta_os = theta_os,
                     theta_ls = theta_ls, vote_codes = NULL, thresholds=NULL,
                     SEED=SEED, constant_IRF = 0)

samples = samples[[1]]
SAMPLE_ITERS = SAMPLE_ITERS/THIN

save.image(file='./results/gpirt_senate_92.RData')

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
  for (j in 1:length(current_unique_icpsrs)) {
    i = idx[j]
    tmp = samples$theta[-1,i,h]
    # if(sd(tmp)<1){
    #   # drop wrong sign
    #   drop_wrong_sign = (sign(cor(nominate_scores[,1],colMeans(samples$theta)[idx,h]))*tmp*nominate_scores[j,1]<0)
    #   tmp = tmp[drop_wrong_sign==0]
    #   if(is.na(mean(tmp))){
    #     tmp = samples$theta[-1,i,h]
    #   }
    # }
    
    pred_theta[i,h] = mean(tmp)
    pred_theta_sd[i,h] = sd(tmp)
  }
}

# # Mark Hatfield
# pred_theta[169,6] = mean(samples$theta[-1,169,6])
# pred_theta_sd[169,6] = sd(samples$theta[-1,169,6])
# 
# pred_theta[169,7] = mean(samples$theta[-1,169,7])
# pred_theta_sd[169,7] = sd(samples$theta[-1,169,7])
# 
# # RUSSELL, Donald Stuart
# pred_theta[163,5] = mean(samples$theta[-1,163,5])
# pred_theta_sd[163,5] = sd(samples$theta[-1,163,5])
# 
# # HOLLINGS, Ernest Frederick
# pred_theta[170,6] = mean(samples$theta[-1,170,6])
# pred_theta_sd[170,6] = sd(samples$theta[-1,170,6])
# 
# pred_theta[170,7] = mean(samples$theta[-1,170,7])
# pred_theta_sd[170,7] = sd(samples$theta[-1,170,7])
# 
# # EDWARDS, Elaine Schwartzenburg
# pred_theta[202,8] = mean(samples$theta[-1,202,8])
# pred_theta_sd[202,8] = sd(samples$theta[-1,202,8])

mask = (samples$theta[-1,60,7]>-3)
pred_theta[60,7] = mean(samples$theta[-1,60,7][mask])
pred_theta_sd[60,7] = sd(samples$theta[-1,60,7][mask])

mask = (samples$theta[-1,60,8]>-3)
pred_theta[60,8] = mean(samples$theta[-1,60,8][mask])
pred_theta_sd[60,8] = sd(samples$theta[-1,60,8][mask])

mask = (samples$theta[-1,60,9]>-3)
pred_theta[60,9] = mean(samples$theta[-1,60,9][mask])
pred_theta_sd[60,9] = sd(samples$theta[-1,60,9][mask])

mask = (samples$theta[-1,60,10]>-3)
pred_theta[60,10] = mean(samples$theta[-1,60,10][mask])
pred_theta_sd[60,10] = sd(samples$theta[-1,60,10][mask])

cor_theta = c()
pred_theta_ll = matrix(NA, nrow=nrow(data), ncol=dim(data)[3])
all_current_unique_icpsrs = list()
for(h in 1:length(session_ids)){
  session_id = session_ids[h]
  members = read.csv(paste("./data/S", session_id, "_members.csv", sep=""))
  members = members[members$chamber=="Senate", c("icpsr", 
                      "nokken_poole_dim1", "nokken_poole_dim2")]
  # exclude BARKLEY, Dean
  members = members[members$icpsr!=40106, ]
  # exclude RUSSELL, Richard Brevard, Jr. of GA
  members = members[members$icpsr!=8138, ]
  # exclude James Danforth
  members = members[members$icpsr!=14447, ]
  # exclude JOHNSTON, Olin DeWitt Talmadge
  members = members[members$icpsr!=5009, ]
  current_unique_icpsrs = unique(members$icpsr)
  all_current_unique_icpsrs[[h]] = current_unique_icpsrs
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
  pred_theta[idx,h] = sign(cor(pred_theta[idx, h],nominate_scores[,1]))*pred_theta[idx,h]
  pred_theta[idx,h] = (pred_theta[idx,h] - mean(pred_theta[idx,h]))
  cor_theta = c(cor_theta, cor(pred_theta[idx, h],nominate_scores[,1]))
  plot(pred_theta[idx,h], nominate_scores[,1])
  pred_theta_ll[idx,h] = log(dnorm(nominate_scores[,1],mean=pred_theta[idx,h],sd=pred_theta_sd[idx,h]))
}

# dynamic theta
all_service_senates = unique_icpsr
for(h in 1:horizon){
  all_service_senates = intersect(all_service_senates, all_current_unique_icpsrs[[h]])
}

for(icpsr in all_service_senates){
  idx = which(icpsr==unique_icpsr)
  session_id = session_ids[horizon]
  members = read.csv(paste("./data/S", session_id, "_members.csv", sep=""))
  bioname =toString(members$bioname[members$icpsr==icpsr])
  jpeg(file=paste('./results/', toString(icpsr) ,'.jpeg' , sep=''))
  plot(1:horizon, pred_theta[idx,])
  title(bioname)
  dev.off()
}

# MANCHIN, Joe, III
# all_service_senates = c(all_service_senates, 40915)

all_service_senate_data = data.frame(matrix(ncol =6, nrow = 0))
colnames(all_service_senate_data) <- c("session", "icpsr", "gpirt", "party", "name", "nominate")
for(icpsr in all_service_senates){
  idx = which(icpsr==unique_icpsr)
  # plot(1:horizon, pred_theta[idx,])
  for(h in 1:horizon){
    session_id = session_ids[h]
    members = read.csv(paste("./data/S", session_id, "_members.csv", sep=""))
    party = members$party_code[members$icpsr==icpsr]
    if(party==100){
      party = "Democratics"
    }else if(party==200){
      party = "Republicans"
    }else{
      party = "Independents"
    }
    bioname =toString(members$bioname[members$icpsr==icpsr])
    all_service_senate_data[nrow(all_service_senate_data)+1,] = c(session_id,
        icpsr, pred_theta[idx,h], party, bioname, members$nokken_poole_dim1[members$icpsr==icpsr])
  }
}

all_service_senate_data$gpirt = as.numeric(all_service_senate_data$gpirt)
all_service_senate_data$nominate = as.numeric(all_service_senate_data$nominate)
all_service_senate_data$norm_gpirt = (all_service_senate_data$gpirt-mean(all_service_senate_data$gpirt))/sd(all_service_senate_data$gpirt)

p = ggplot(all_service_senate_data, 
       aes(x=session,y=gpirt,group=factor(icpsr), color=factor(party)))+
  scale_y_continuous(name="GPIRT") + 
  geom_line() 

png("./results/senate_dynamic_scores_92.png")
print(p)
dev.off()

# ordinal_lls = function(f, thresholds){
#   result = c()
#   for (c in 1:(length(thresholds)-1)) {
#     z1 = thresholds[c] - f;
#     z2 = thresholds[c+1] - f;
#     result = c(result,log(pnorm(z2)-pnorm(z1)));
#   }
#   return(result)
# }
# 
# pred_lls = c()
# pred_acc = c()
# train_lls = c()
# train_acc = c()
# 
# sample_IRFs = vector(mode = "list", length = SAMPLE_ITERS)
# for (iter in 1:SAMPLE_ITERS){
#   sample_IRFs[[iter]] = recover_fstar(BURNOUT_ITERS+iter-1,samples$f[[iter]],data, 
#             as.matrix(samples$theta[iter,,]), samples$threshold[iter,])$fstar
# }
# 
# for (i in 1:n) {
#   for (j in 1:m) {
#     for(h in 1:horizon){
#       if(!is.na(data[[i,j,h]])){
#         lls = matrix(0,nrow=SAMPLE_ITERS, ncol = C)
#         y_pred = rep(0, SAMPLE_ITERS)
#         pred_idx = 1+as.integer((pred_theta[i,h]+5)*100)
#         for (iter in 1:SAMPLE_ITERS) {
#           f_pred = sample_IRFs[[iter]][pred_idx, j, h]
#           ll = ordinal_lls(f_pred, samples$threshold[iter,])
#           lls[iter,] = ll
#           y_pred[iter] =  which.max(ll)
#         }
#         ll = log(colMeans(exp(lls)))
#         y_pred = round(mean(y_pred))
#         if(train_idx[i,j,h]==0 & !is.na(data[[i,j,h]])){
#           pred_acc = c(pred_acc, y_pred==(data[[i,j, h]]))
#           pred_lls = c(pred_lls, ll[data[[i,j, h]]])
#         }else{
#           train_acc = c(train_acc, y_pred==(data[[i,j, h]]))
#           train_lls = c(train_lls, ll[data[[i,j, h]]])
#         }
#       }
#     }
#   }
# }
# 
# # train and preditive log likelihood and accuracies
# print(mean(train_lls[!is.infinite(train_lls)]))
# print(mean(train_acc[!is.infinite(train_lls)]))
# print(mean(pred_lls[!is.infinite(pred_lls)]))
# print(mean(pred_acc[!is.infinite(pred_lls)]))

# gp IRFs
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
    }
    probs = getprobs_gpirt(xs[idx], t(IRFs),
                           samples$threshold)
    tmp = probs %>%
      group_by(xs) %>%
      summarize(icc=sum(order*p))
    gpirt_iccs[,j,h] = tmp$icc
  }
}

save.image(file='./results/gpirt_senate_92.RData')
