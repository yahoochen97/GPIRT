# gpirt_path = "~/Documents/Github/gpirt"
# setwd(gpirt_path)
# library(Rcpp)
# Rcpp::compileAttributes()
# install.packages(gpirt_path, type="source", repos = NULL)#,lib=R_path, INSTALL_opts = '--no-lock')
# setwd("../OrdGPIRT")

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
# beta_prior_sds =  matrix(1.0, nrow = 3, ncol = ncol(rollcall_data))
# beta_proposal_sds =  matrix(0.1, nrow = 3, ncol = ncol(rollcall_data))
beta_prior_means = matrix(0, nrow = 3, ncol = m)
beta_prior_sds =  matrix(1.0, nrow = 3, ncol = m)
samples <- gpirtMCMC(rollcall_data, SAMPLE_ITERS,BURNOUT_ITERS,
                     THIN, CHAIN, beta_prior_means = beta_prior_means,theta_init = theta_init,
                     beta_prior_sds = beta_prior_sds, theta_os = theta_os,
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
  cor_theta_2 = c(cor_theta_2, cor(pred_theta[idx, h],nominate_scores[,2]))
  plot(pred_theta[idx,h], nominate_scores[,2])
  pred_theta_ll[idx,h] = log(dnorm(nominate_scores[,2],mean=pred_theta[idx,h],sd=pred_theta_sd[idx,h]))
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

xs = seq(-5,5,0.01)
idx = 201:801
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
    gpirt_iccs[,j,h] = probs$p[probs$order==2]
  }
}

plot(xs[idx],gpirt_iccs[,1,1], ylim=c(1,2))
mask = is.na(rollcall_data[,1,1])
points(pred_theta[!mask,1], rollcall_data[!mask,1,1])

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

# all_service_senates = all_senator_ids
# for(h in 1:horizon){
#   congress == congresses[h]
#   rollcall_ids = unique(data[data$congress==congress, "rollcall"]$rollcall)
#   senator_ids = unique(data[data$congress==congress,"id"]$id)
#   idx = c()
#   for(j in 1:length(senator_ids)){
#     icpsr = senator_ids[j]
#     i = which(icpsr==all_senator_ids)
#     mask = !is.na(rollcall_data[i,,h])
#     if(sum(mask)>0){
#       idx = c(idx, i)
#     }
#   }
#   all_service_senates = intersect(all_service_senates, all_senator_ids[idx])
# }


# 1252, 10818, 14101, 14503, 14511, 14713, 14904
# 14101, 14511, 12032, 14226,
all_service_senates = c(1366, 10808, 14101, 14503, 14511, 14713)

dynamic_score_data = data.frame(matrix(ncol = 5, nrow = 0))
colnames(dynamic_score_data) <- c("session","score", "type", "icpsr", "bioname")
for(h in 1:length(congresses)){
  session_id = congresses[h]
  members = read.csv(paste("./data/S", session_id, "_members.csv", sep=""))
  members = members[members$chamber=="Senate",]
  senator_ids = intersect(unique(members$icpsr), all_service_senates)
  tmp = nominate_theta[match(senator_ids,all_senator_ids),h]
  mask = !is.na(tmp)
  tmp = data.frame(tmp[mask])
  colnames(tmp) = c("score")
  tmp$type = "NOMINATE"
  tmp$session = congresses[h]
  tmp$icpsr = senator_ids[mask]
  tmp$bioname = as.character(members[members$icpsr %in% senator_ids[mask], "bioname"])
  
  for(k in 1:length(senator_ids[mask])){
    icpsr = senator_ids[mask][k]
    tmp$bioname[k] = as.character(members[members$icpsr==icpsr, "bioname"])
  }
  dynamic_score_data= rbind(dynamic_score_data, tmp)
  
  tmp$score = pred_theta[match(senator_ids[mask],all_senator_ids),h]/sd(pred_theta)*sd(nominate_theta[!is.na(nominate_theta)])
  tmp$type = "GPIRT"
  dynamic_score_data= rbind(dynamic_score_data, tmp)
}

write.csv(dynamic_score_data, file="./results/gpirt_abortion_dynamic.csv")


# for(id in all_service_senates){
#   plot(congresses, pred_theta[which(id==all_senator_ids),])
#   plot(congresses, nominate_theta[which(id==all_senator_ids),])
# }

xs = seq(-5,5,0.01)
idx = 201:801
plot(xs[idx],gpirt_iccs[,1,10], ylim=c(0,1))
mask = is.na(rollcall_data[,1,10])
points(pred_theta[!mask,10], rollcall_data[!mask,1,10]-1)

# plot irf
folder_path = "./figures/abortion/"
dir.create(file.path(folder_path), showWarnings = FALSE)

for(h in 1:horizon){
  congress = congresses[h]
  rollcall_ids = unique(data[data$congress==congress, "rollcall"]$rollcall)
  senator_ids = unique(data[data$congress==congress,"id"]$id)
  subfolder = as.character(congresses[h])
  dir.create(file.path(folder_path,subfolder), showWarnings = FALSE)
  for(j in 1:length(rollcall_ids)){
    rid = rollcall_ids[j]
    idx = c()
    for(i in 1:length(senator_ids)){
      icpsr = senator_ids[i]
      idx = c(idx, which(icpsr==all_senator_ids))
    }
    x = pred_theta[idx,h]
    response = rollcall_data[idx,j,h]
    irf_plot = data.frame(x,response-1)
    colnames(irf_plot) = c("x" , "response")
    xs = seq(-5,5,0.01)
    idx = 201:801
    gpirt_plot = data.frame(xs[idx],gpirt_iccs[,j,h])
    colnames(gpirt_plot) = c("xs","icc")
    p = ggplot()+
      geom_point(data = na.omit(irf_plot), aes(x=x,y=response,color=factor(response)),
                 size=4, shape="|") +
      scale_color_manual(name='response',
                         labels=c('Nay', 'Yea'),
                         values=c('black', 'red'))+
      scale_x_continuous(name="x",breaks = seq(-2, 2, by = 1)) + 
      # name=bquote(theta),
      scale_y_continuous(name="P(yea)", breaks=seq(0, 1, by = 0.2)) +
      geom_line(data = gpirt_plot, aes(x=xs,y=icc)) +
      theme(panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            axis.title.x=element_text(),
            axis.title.y=element_text(),
            legend.position = "none")
    
    ggsave(filename = paste(folder_path, subfolder, "/", as.character(rid), ".png",sep = ""),width = 3, height = 2, dpi = 300)
    # png(paste(folder_path, subfolder, "/", as.character(rid), ".png",sep = ""))
    # print(p)
    # dev.off()
  }
}

print(data[data$congress==100 & data$rollcall==653,"mynotes"][1,])

save.image(file='./results/gpirt_abortion.RData')

# GPIRT
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
      # tmp = sign(cor(samples$fstar[[iter]][, j, h],samples$fstar[[1]][, j, h]))
      IRFs[iter, ] = samples$fstar[[iter]][, j, h]# *tmp
    }
    probs = getprobs_gpirt(xs, t(IRFs), samples$threshold)
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
