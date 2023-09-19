library(gpirt)
library(dplyr)
library(ggplot2)
library(stats)

gpirt_path = "~/Documents/Github/OrdGPIRT"
setwd(gpirt_path)
load(file="./data/senate_data_90.RData")
SAMPLE_ITERS = 500
BURNOUT_ITERS = 500
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
THIN = 4
CHAIN = 1
beta_prior_sds =  matrix(1.0, nrow = 2, ncol = ncol(data))
beta_proposal_sds =  matrix(0.1, nrow = 2, ncol = ncol(data))
samples <- gpirtMCMC(data, SAMPLE_ITERS,BURNOUT_ITERS,
                     THIN, CHAIN, theta_init = theta_init,
                     beta_prior_sds = beta_prior_sds, theta_os = theta_os,
                     theta_ls = theta_ls, vote_codes = NULL, thresholds=NULL,
                     SEED=SEED, constant_IRF = 0)

samples = samples[[1]]
SAMPLE_ITERS = SAMPLE_ITERS/THIN

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

###################################
# effective sample size diagonis
library(coda)

ESS_theta = matrix(0, nrow = n, ncol=h)
for (h in 1:horizon) {
  for(i in 1:n){
    ESS_theta[i,h] = coda::effectiveSize(samples$theta[,i,h])
  }
}

#####################################

xs = seq(-5,5,0.01)
for(it in 1:SAMPLE_ITERS){
  for(h in 1:horizon){
    for(j in 1:m){
      samples$f[[it]][,j,h] = samples$f[[it]][,j,h] + samples$beta[[it]][1,j,h] + samples$beta[[it]][2,j,h]*samples$theta[it,,h]
      samples$fstar[[it]][,j,h] = samples$fstar[[it]][,j,h] + samples$beta[[it]][1,j,h] + samples$beta[[it]][2,j,h]*xs
    }
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

#####################################

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

# mask = (samples$theta[-1,60,7]>-3)
# pred_theta[60,7] = mean(samples$theta[-1,60,7][mask])
# pred_theta_sd[60,7] = sd(samples$theta[-1,60,7][mask])
# 
# mask = (samples$theta[-1,60,8]>-3)
# pred_theta[60,8] = mean(samples$theta[-1,60,8][mask])
# pred_theta_sd[60,8] = sd(samples$theta[-1,60,8][mask])
# 
# mask = (samples$theta[-1,60,9]>-3)
# pred_theta[60,9] = mean(samples$theta[-1,60,9][mask])
# pred_theta_sd[60,9] = sd(samples$theta[-1,60,9][mask])
# 
# mask = (samples$theta[-1,60,10]>-3)
# pred_theta[60,10] = mean(samples$theta[-1,60,10][mask])
# pred_theta_sd[60,10] = sd(samples$theta[-1,60,10][mask])

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

# for(icpsr in all_service_senates){
#   idx = which(icpsr==unique_icpsr)
#   session_id = session_ids[horizon]
#   members = read.csv(paste("./data/S", session_id, "_members.csv", sep=""))
#   bioname =toString(members$bioname[members$icpsr==icpsr])
#   jpeg(file=paste('./results/', toString(icpsr) ,'.jpeg' , sep=''))
#   plot(1:horizon, pred_theta[idx,])
#   title(bioname)
#   dev.off()
# }

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

png("./results/senate_dynamic_scores_90.png")
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

# read data
library(haven)
abortion_data = read_dta("./data/SenatePeriods.dta")
# remove president
abortion_data = abortion_data[abortion_data$name!="REAGAN",]
abortion_data = abortion_data[abortion_data$name!="BUSH",]
abortion_data = abortion_data[abortion_data$name!="CLINTON",]

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
  
  # abortion_rollcalls = unique(abortion_data[abortion_data$congress==session_id, c("rollcall")])
  # rollcalls = rollcalls[!(rollcalls$rollnumber %in% abortion_rollcalls),]
  rollcall_ids = unique(rollcalls$rollnumber)
  for (j in 1:length(rollcall_ids)) {
    print(j)
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
    gpirt_iccs[,j,h] = probs$p[probs$order==2]
  }
}

# plot irf
folder_path = "./figures/senate/"
dir.create(file.path(folder_path), showWarnings = FALSE)

for(h in 1:horizon){
  congress = session_ids[h]
  session_id = congress
  rollcalls = read.csv(paste("./data/S", session_id, "_rollcalls.csv", sep=""))
  rollcalls = rollcalls[,c("congress", "rollnumber", "yea_count","nay_count","date" )]
  rollcalls = rollcalls[(rollcalls$yea_count!=0)&(rollcalls$nay_count!=0),]
  
  # abortion_rollcalls = unique(abortion_data[abortion_data$congress==session_id, c("rollcall")])
  # rollcalls = rollcalls[!(rollcalls$rollnumber %in% abortion_rollcalls),]
  rollcall_ids = unique(rollcalls$rollnumber)
  
  members = read.csv(paste("./data/S", session_id, "_members.csv", sep=""))
  members = members[members$chamber=="Senate", c("icpsr", 
                                                 "nokken_poole_dim1", "nokken_poole_dim2")]
  
  senator_ids = unique(members$icpsr)
  subfolder = as.character(session_ids[h])
  dir.create(file.path(folder_path,subfolder), showWarnings = FALSE)
  for(j in 1:length(rollcall_ids)){
    rid = rollcall_ids[j]
    idx = c()
    for(i in 1:length(senator_ids)){
      icpsr = senator_ids[i]
      idx = c(idx, which(icpsr==unique_icpsr))
    }
    x = pred_theta[idx,h]
    response = data[idx,j,h] - 1
    irf_plot = data.frame(x,response)
    xs = seq(-5,5,0.01)
    idx = 201:801
    gpirt_plot = data.frame(xs[idx],gpirt_iccs[,j,h])
    colnames(gpirt_plot) = c("xs","icc")
    p = ggplot()+
      geom_point(data = na.omit(irf_plot), aes(x=x,y=response,color=factor(response)),
                 size=2, shape="|") +
      scale_color_manual(name='response',
                         labels=c('Nay', 'Yea'),
                         values=c('black', 'red'))+
      scale_y_continuous(name="P(yea)") +
      geom_line(data = gpirt_plot, aes(x=xs,y=icc))
    
    png(paste(folder_path, subfolder, "/", as.character(rid), ".png",sep = ""))
    print(p)
    dev.off()
  }
}

save.image(file='./results/gpirt_senate_90.RData')

# train/test statistics
# train
train_lls = c()
train_acc = c()
train_response = c()
train_prediction = c()

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
      # if (cor(samples$fstar[[1]][, j, h],samples$fstar[[iter]][, j, h])>0){
      #   mask[iter] = 1
      # }
    }
    # if(sum(mask)>=sum(!mask)){
    #   probs = getprobs_gpirt(xs, t(IRFs[mask==1,]), samples$threshold[mask==1,])
    # }
    # else{
    #   probs = getprobs_gpirt(xs, t(IRFs[mask==0,]), samples$threshold[mask==0,])
    # }
    probs = getprobs_gpirt(xs, t(IRFs), samples$threshold)
    # tmp = probs %>%
    #  group_by(xs) %>%
    #  summarize(icc=sum(order*p))
    # gpirt_iccs[,j,h] = probs$p[probs$order==2]
    
    for (i in 1:n) {
      if(!is.na(data[i,j,h])){
        # train
        pred_idx = 1+as.integer((pred_theta[i,h]+5)*100)
        ll = log(probs$p[probs$xs==xs[pred_idx]])
        y_pred = which.max(ll)
        train_acc = c(train_acc, y_pred==(data[i,j, h]))
        train_lls = c(train_lls, ll[data[i,j, h]])
        train_response = c(train_response, data[i,j, h])
        train_prediction = c(train_prediction, y_pred)
      }
    }
  }
}

results = data.frame(train_acc,train_lls, 
                     train_response, train_prediction)

colnames(results) = c("train_acc", "train_lls", "train_response", "train_prediction")
write.csv(results, "./results/gpirt_senate_fit.csv")

results_doirt = read.csv("./results/doirt_senate_fit.csv")

library(cvAUC)
print(paste("gpirt avg acc: ", round(mean(train_acc),3),sep=""))
print(paste("doirt avg acc: ", round(mean(results_doirt$train_acc),3),sep=""))
print(paste("t test p value: ", t.test(train_acc,results_doirt$train_acc)[["p.value"]],sep=""))

print(paste("gpirt avg ll: ", round(mean(train_lls),4),sep=""))
print(paste("doirt avg ll: ", round(mean(results_doirt$train_lls),4),sep=""))
print(paste("t test p value: ", t.test(train_lls,results_doirt$train_lls)[["p.value"]],sep=""))

print(paste("gpirt auc: ", round(AUC(train_prediction, train_response),3),sep=""))
print(paste("doirt auc: ", round(AUC(results_doirt$train_prediction, results_doirt$train_response),3),sep=""))
