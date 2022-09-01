gpirt_path = "~/Documents/Github/gpirt"
setwd(gpirt_path)
library(Rcpp)
Rcpp::compileAttributes()
install.packages(gpirt_path, type="source", repos = NULL)#,lib=R_path, INSTALL_opts = '--no-lock')
setwd("../OrdGPIRT")

library(gpirt)
library(dplyr)
library(ggplot2)
library(stats)
library(haven)

gpirt_path = "~/Documents/Github/OrdGPIRT"
setwd(gpirt_path)
TYPE = "GP"

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
  # votes = read.csv(paste("./data/S", toString(congress), "_votes.csv", sep=""))
  rollcall_ids = unique(data[data$congress==congress, "rollcall"]$rollcall)
  senator_ids = unique(data[data$congress==congress,"id"]$id)
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

# n by horizon
theta_init = matrix(0, nrow = n, ncol = horizon)
nominate_theta =  array(array(NA, n*horizon), c(n, horizon))
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
  theta_init[idx,h] = nominate_scores[,1] + 0.1*rnorm(length(idx))
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

SAMPLE_ITERS = 100
BURNOUT_ITERS = 100
SEED = 1
THIN = 1
CHAIN = 1
beta_prior_sds =  matrix(0.5, nrow = 2, ncol = ncol(rollcall_data))
beta_proposal_sds =  matrix(0.1, nrow = 2, ncol = ncol(rollcall_data))
samples <- gpirtMCMC(rollcall_data, SAMPLE_ITERS,BURNOUT_ITERS,
                     THIN, CHAIN, beta_prior_sds = beta_prior_sds,theta_init = theta_init,
                     beta_proposal_sds = beta_proposal_sds, theta_os = theta_os,
                     theta_ls = theta_ls, vote_codes = NULL, thresholds=NULL,
                     SEED=SEED, constant_IRF = 1)

samples = samples[[1]]
SAMPLE_ITERS = SAMPLE_ITERS/THIN

save.image(file='./results/gpirt_abortion.RData')

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

cor_theta_1 = c()
pred_theta_ll = matrix(NA, nrow=n, ncol=horizon)
for(h in 1:length(congresses)){
  session_id = congresses[h]
  members = read.csv(paste("./data/S", session_id, "_members.csv", sep=""))
  members = members[members$chamber=="Senate", c("icpsr", 
                                                 "nokken_poole_dim1", "nokken_poole_dim2")]
  senator_ids = unique(data[data$congress==session_id,"id"]$id)
  # nominate scores 
  nominate_scores = matrix(0, nrow=length(senator_ids), ncol=2)
  idx = c()
  for(j in 1:length(senator_ids)){
    icpsr = senator_ids[j]
    nominate_scores[j,1] = members[members$icpsr==icpsr, "nokken_poole_dim1"]
    nominate_scores[j,2] = members[members$icpsr==icpsr, "nokken_poole_dim2"]
    idx = c(idx, which(icpsr==all_senator_ids))
  }
  
  # pred_theta[idx,h] = sign(cor(pred_theta[idx, h],nominate_scores[,1]))*pred_theta[idx,h]
  # pred_theta[idx,h] = (pred_theta[idx,h] - mean(pred_theta[idx,h]))
  cor_theta_1 = c(cor_theta_1, cor(pred_theta[idx, h],nominate_scores[,1]))
  plot(pred_theta[idx,h], nominate_scores[,1])
  pred_theta_ll[idx,h] = log(dnorm(nominate_scores[,1],mean=pred_theta[idx,h],sd=pred_theta_sd[idx,h]))
}

cor_theta_2 = c()
pred_theta_ll = matrix(NA, nrow=n, ncol=horizon)
for(h in 1:length(congresses)){
  session_id = congresses[h]
  members = read.csv(paste("./data/S", session_id, "_members.csv", sep=""))
  members = members[members$chamber=="Senate", c("icpsr", 
                                                 "nokken_poole_dim1", "nokken_poole_dim2")]
  senator_ids = unique(data[data$congress==session_id,"id"]$id)
  # nominate scores 
  nominate_scores = matrix(0, nrow=length(senator_ids), ncol=2)
  idx = c()
  for(j in 1:length(senator_ids)){
    icpsr = senator_ids[j]
    nominate_scores[j,1] = members[members$icpsr==icpsr, "nokken_poole_dim1"]
    nominate_scores[j,2] = members[members$icpsr==icpsr, "nokken_poole_dim2"]
    idx = c(idx, which(icpsr==all_senator_ids))
  }
  
  # pred_theta[idx,h] = sign(cor(pred_theta[idx, h],nominate_scores[,2]))*pred_theta[idx,h]
  # pred_theta[idx,h] = (pred_theta[idx,h] - mean(pred_theta[idx,h]))
  cor_theta_2 = c(cor_theta_2, cor(pred_theta[idx, h],nominate_scores[,2]))
  plot(pred_theta[idx,h], nominate_scores[,2])
  pred_theta_ll[idx,h] = log(dnorm(nominate_scores[,2],mean=pred_theta[idx,h],sd=pred_theta_sd[idx,h]))
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

idx = 301:701
plot(xs[idx],gpirt_iccs[,1,1], ylim=c(1,2))
mask = is.na(rollcall_data[,1,1])
points(pred_theta[!mask,1], rollcall_data[!mask,1,1])

save.image(file='./results/gpirt_abortion.RData')

all_nominate_data = data.frame(matrix(ncol = 6, nrow = 0))
colnames(all_nominate_data) <- c("session", "gpirt", "nominate", "party", "icpsr", "bioname")
for(h in 1:length(congresses)){
  session_id = congresses[h]
  members = read.csv(paste("./data/S", session_id, "_members.csv", sep=""))
  members = members[members$chamber=="Senate",]
  # exlude BARKLEY, Dean
  members = members[members$icpsr!=40106, ]
  # excluse RUSSELL, Richard Brevard, Jr. of GA
  members = members[members$icpsr!=8138, ]
  # exclude JOHNSTON, Olin DeWitt Talmadge
  members = members[members$icpsr!=5009, ]
  # exclude James Danforth
  # members = members[members$icpsr!=14447, ]
  senator_ids = unique(data[data$congress==session_id,"id"]$id)
  # nominate scores 
  nominate_scores = matrix(0, nrow=length(senator_ids), ncol=2)
  idx = c()
  bionames = c()
  party_codes = c()
  for(j in 1:length(senator_ids)){
    icpsr = senator_ids[j]
    nominate_scores[j,1] = members[members$icpsr==icpsr, "nokken_poole_dim1"]
    nominate_scores[j,2] = members[members$icpsr==icpsr, "nokken_poole_dim2"]
    idx = c(idx, which(icpsr==all_senator_ids))
    bionames = c(bionames, toString(members[members$icpsr==icpsr, "bioname"]))
    party_codes = c(party_codes, members[members$icpsr==icpsr, "party_code"])
  }
  current_pred_theta = sign(cor(pred_theta[idx, h],nominate_scores[,1]))*pred_theta[idx, h]
  nominate_data = data.frame(current_pred_theta, nominate_scores[,1])
  colnames(nominate_data) = c("gpirt", "nominate")
  party_codes[(party_codes!=200)&(party_codes!=100)] = "Independents"
  party_codes[party_codes==100] = "Democrats"
  party_codes[party_codes==200] = "Republicans"
  nominate_data$party = party_codes
  nominate_data$session = session_id
  nominate_data$icpsr = senator_ids
  nominate_data$bioname = bionames
  all_nominate_data = rbind(all_nominate_data, nominate_data)
  # p[[h]] = ggplot(nominate_data, aes(y=x, x=y, colour = factor(party))) +
  #   geom_point(size=2, aes(shape=factor(party))) + 
  #   xlab("NOMINATE Dimension 1 Ideology") + ylab("GPIRT Ideology") + 
  #   labs(colour = "Party")
}

write.csv(all_nominate_data, file="./results/gpirt_abortion_results.csv")

all_service_senates = all_senator_ids
for(h in 1:horizon){
  senator_ids = unique(data[data$congress==congress,"id"]$id)
  all_service_senates = intersect(all_service_senates, senator_ids)
}

for(id in all_service_senates){
  plot(congresses, pred_theta[3,])
}
