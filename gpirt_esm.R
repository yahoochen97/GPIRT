R_path="~/R/x86_64-redhat-linux-gnu-library/4.0"
.libPaths(R_path)

# gpirt_path = "~/Documents/Github/gpirt"
# setwd(gpirt_path)
# library(Rcpp)
# Rcpp::compileAttributes()
# install.packages(gpirt_path, type="source", repos = NULL)
# setwd("../OrdGPIRT")

library(gpirt)
library(dplyr)
library(ggplot2)
library(stats)

SAMPLE_ITERS = 5000
BURNOUT_ITERS = 5000
TYPE = "GP"

source("2PL_esm.R")

# # esm.BFI09, esm.BFI19 and esm.BFI04 
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
# 
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
# 
# data = ems_data

if(TYPE=="GP"){
  theta_os = 1
  theta_ls = 8
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
THIN = 4
CHAIN = 3
constant_IRF = 0
beta_prior_sds =  matrix(1.0, nrow = 2, ncol = ncol(data))
beta_proposal_sds =  matrix(0.1, nrow = 2, ncol = ncol(data))
all_samples <- gpirtMCMC(data, SAMPLE_ITERS,BURNOUT_ITERS,
                     THIN, CHAIN,
                     beta_proposal_sds = beta_proposal_sds, theta_init = theta_init,
                     beta_prior_sds = beta_prior_sds, theta_os = theta_os,
                     theta_ls = theta_ls, vote_codes = NULL, thresholds=NULL,
                     SEED=SEED, constant_IRF = constant_IRF)

save.image(file='./results/gpirt_esm.RData')

SAMPLE_ITERS = SAMPLE_ITERS/THIN
library(rstan)
sims <- matrix(rnorm((SAMPLE_ITERS)*CHAIN), nrow = SAMPLE_ITERS, ncol = CHAIN)
theta_rhats = matrix(rnorm(n*horizon), nrow = n, ncol = horizon)
for(i in 1:n){
  for(h in 1:horizon){
    for(c in 1:CHAIN){
      sims[,c] = all_samples[[c]]$theta[,i,h]
    }
    theta_rhats[i, h] = Rhat(sims)
  }
}

samples = all_samples[[1]]



xs = seq(-5,5,0.01)
for(it in 1:SAMPLE_ITERS){
  for(h in 1:horizon){
    for(j in 1:m){
      samples$f[[it]][,j,h] = samples$f[[it]][,j,h] + samples$beta[[it]][1,j,h] + samples$beta[[it]][2,j,h]*samples$theta[it,,h]
      samples$fstar[[it]][,j,h] = samples$fstar[[it]][,j,h] + samples$beta[[it]][1,j,h] + samples$beta[[it]][2,j,h]*xs
    }
  }
}

# predicted ideology
pred_theta = matrix(0, nrow=nrow(data), ncol=dim(data)[3])
pred_theta_sd = matrix(0, nrow=nrow(data), ncol=dim(data)[3])

for(h in 1:horizon){
  for (i in 1:n) {
    tmp = samples$theta[-1,i,h]
    # mask = rep(1, SAMPLE_ITERS)
    # for(iter in 1:SAMPLE_ITERS){
    #   mask[iter, ] = samples$fstar[[it]][,j,h] * sign(cor(samples$fstar[[it]][,j,h], grm_together_iccs[,j]))
    # }
    pred_theta[i,h] = mean(tmp)
    pred_theta_sd[i,h] = sd(tmp)
  }
}

cor_theta = c()
for(h in 1:horizon){
  cor_theta = c(cor_theta, cor(pred_theta[,h],together_pred_theta[,h]))
}


tmp = c()
for(it in 1:SAMPLE_ITERS){
  tmp = c(tmp, samples$beta[[it]][2,1,1])
}
plot(1:SAMPLE_ITERS, tmp)

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
      IRFs[iter, ] = samples$fstar[[iter]][idx, j, h]
      IRFs[iter, ] = IRFs[iter, ] * sign(cor(IRFs[iter, ], grm_together_iccs[,j]))
    }
    thresholds = matrix(0,nrow=SAMPLE_ITERS,ncol=C+1)
    for(iter in 1:SAMPLE_ITERS){
      thresholds[iter, ] = samples$threshold[[iter]][j,,h]
    }
    probs = getprobs_gpirt(xs[idx], t(IRFs), thresholds)
    tmp = probs %>%
      group_by(xs) %>%
      summarize(icc=sum(order*p))
    gpirt_iccs[,j,h] = tmp$icc
  }
}

for(h in 1:horizon){
  all_C = TRUE
  for(j in 1:m){
    if(length(na.omit(unique(data[,j,h])))<5){
      all_C = FALSE
    }
  }
  if(all_C){
    print(h)
  }
}

folder_path = "./figures/esm/"
dir.create(file.path(folder_path), showWarnings = FALSE)

for(h in 1:horizon){
  subfolder = as.character(h)
  dir.create(file.path(folder_path,subfolder), showWarnings = FALSE)
  for(j in 1:m){
    x = pred_theta[,h]
    response = data[,j,h]
    irf_plot = data.frame(x,response)
    xs = seq(-5,5,0.01)
    idx = 401:601
    gpirt_plot = data.frame(xs[idx],gpirt_iccs[101:301,j,h])
    colnames(gpirt_plot) = c("xs","icc")
    p = ggplot()+
      geom_point(data = na.omit(irf_plot), aes(x=x,y=response,color=factor(response)),
                 size=4, shape="|") +
      scale_color_manual(name='Level',
                         labels=c("Strong Disagree","Disagree", "Neural", "Agree", "Stronly Dgree"),
                         values=c('red','red', 'blue', 'black', 'black'))+
      geom_line(data = gpirt_plot, aes(x=xs,y=icc), size=1)+
      scale_x_continuous(name=bquote(theta), breaks = seq(-2, 2, by = 1)) + 
      scale_y_discrete(name=NULL, limits=c("Strong Disagree","Disagree", "Neural", "Agree", "Stronly Agree")) +
      theme(panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=2),
            legend.position = "none",
            axis.text.y = element_text(size=12,colour = "black"),
            axis.text.x = element_text(size=12,colour = "black"),
            axis.title.x=element_text(size=20,face="bold",colour = "black"),
            plot.title = element_text(hjust = 0.5)) +
      ggtitle(irf_names[j])
    ggsave(filename = paste(folder_path, subfolder, "/", as.character(j), ".png",sep = ""),width = 4, height = 3, dpi = 300)
  }
}

idx = 301:701
for(j in 1:m){
  IRFs = matrix(0, nrow=SAMPLE_ITERS, ncol=length(idx))
  for(iter in 1:SAMPLE_ITERS){
    IRFs[iter, ] = samples$fstar[[iter]][idx, j, 1]
    IRFs[iter, ] = IRFs[iter, ] * sign(cor(IRFs[iter, ], grm_together_iccs[,j]))
  }
  thresholds = matrix(0,nrow=SAMPLE_ITERS,ncol=C+1)
  for(iter in 1:SAMPLE_ITERS){
    thresholds[iter, ] = samples$threshold[[iter]][j,,h]
  }
  probs = getprobs_gpirt(xs[idx], t(IRFs),thresholds)
  
q = ggplot(probs, aes(x=xs, y=p, group=order, color=factor(order))) +
  geom_line(size=2) +ggtitle(paste("IRT q",j, sep="")) +
  theme(plot.title = element_text(hjust = 0.5))
print(q)
# ggsave(paste("./figures/trueirfq",l,".pdf", sep=""), plot=q, width = 7, height = 4, units = "in")
}

save.image(file='./results/gpirt_esm.RData')
