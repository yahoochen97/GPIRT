gpirt_path = "~/Documents/Github/OrdGPIRT"
setwd(gpirt_path)
# library(Rcpp)
# Rcpp::compileAttributes()
# install.packages(gpirt_path, type="source", repos = NULL)# ,lib=R_path, INSTALL_opts = '--no-lock')
setwd("./TAPS")

library(gpirt)
library(dplyr)
library(ggplot2)
library(stats)
library(haven)
library(ltm)

source("load_TAPS.R")

theta_os = 1
theta_ls = 12 # length scale is set to a year

SEED = 1
SAMPLE_ITERS = 100
BURNOUT_ITERS = 100
THIN = 1
CHAIN = 1
beta_prior_sds =  matrix(3.0, nrow = 3, ncol = ncol(gpirt_data))
theta_prior_sds =  matrix(1.0, nrow = 2, ncol = nrow(gpirt_data))
theta_prior_sds[2,] = 0
beta_prior_sds[1,] = 0.0
beta_prior_sds[3,] = 0.0
samples_all <- gpirtMCMC(gpirt_data, SAMPLE_ITERS,BURNOUT_ITERS,
                         THIN, CHAIN, theta_init = theta_init,
                         beta_prior_sds = beta_prior_sds,
                         theta_prior_sds = theta_prior_sds,
                         theta_os = theta_os, theta_ls = theta_ls,
                         vote_codes = NULL, thresholds=NULL,
                         SEED=SEED, constant_IRF = 1)

SAMPLE_ITERS = SAMPLE_ITERS/THIN
samples = samples_all[[1]]

save.image(file='gpirt_TAPS_2014.RData')

###################################
# effective sample size diagonis
library(coda)

ESS_theta = matrix(0, nrow = n, ncol=h)
for (h in 1:horizon) {
  mask = rep(0, SAMPLE_ITERS)
  for(iter in 1:SAMPLE_ITERS){
    if(cor(samples$theta[1,,h],  samples$theta[iter,,h])>0){
      mask[iter] = 1
    }
  }
  for(i in 1:n){
    if(sum(mask)>=sum(!mask)){
      tmp = samples$theta[mask==1,i,h]
    }
    else{
      tmp = samples$theta[mask==0,i,h]
    }
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

pred_theta = matrix(0, nrow=n, ncol=horizon)
pred_theta_sd = matrix(0, nrow=n, ncol=horizon)
for (h in 1:horizon) {
  mask = rep(0, SAMPLE_ITERS)
  for(iter in 1:SAMPLE_ITERS){
    if(cor(samples$theta[1,,h],  samples$theta[iter,,h])>0){
      mask[iter] = 1
    }
  }
  for(i in 1:n){
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

xs = seq(-5,5,0.01)
source("../getprob_gpirt.R")
gpirt_iccs = array(array(0, length(xs)*m*horizon),
                   c(length(xs),m, horizon))

C=5

for (h in 1:horizon) {
  wave = unique(data$wave)[h]
  for (j in 1:m) {
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
    gpirt_iccs[,j,h] = tmp$icc
  }
}

folder_path = "../figures/TAPS2014/"
dir.create(file.path(folder_path), showWarnings = FALSE)

idx = sum(xs<=min(pred_theta)):sum(xs<=max(pred_theta))
for(h in 1:horizon){
  wave = unique(data$wave)[h]
  subfolder = as.character(wave)
  dir.create(file.path(folder_path,subfolder), showWarnings = FALSE)
  for(j in 1:m){
    caseName = questions[j]
    x = -pred_theta[,h] 
    response = gpirt_data[,j,h]
    irf_plot = data.frame(x,response) # negate x 
    gpirt_plot = data.frame(-xs[idx],gpirt_iccs[idx,j,h]) # negate xs
    colnames(gpirt_plot) = c("xs","icc")
    p = ggplot()+
      geom_point(data = na.omit(irf_plot), aes(x=x,y=response,color=factor(response)),
                 size=4, shape="|") +
      scale_color_manual(name='Level',
                         labels=c("Strongly\ndisagree", "Disagree", "Neutral", "Agree", "Strongly\nagree"),
                         values=c('black', 'black','blue', 'red','red'))+
      geom_line(data = gpirt_plot, aes(x=xs,y=icc), size=1)+
      scale_x_continuous(name="x", breaks = seq(-5, 5, by = 1)) +
      scale_y_continuous(name = "Response",
                         breaks = NULL,
                         labels = NULL,
                         limits = c(1,5.2),
                         sec.axis = sec_axis(~.,
                                             breaks = 1:5,
                                             labels = c("Strongly\ndisagree", "Disagree", "Neutral", "Agree", "Strongly\nagree")))+
      theme(panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=2),
            legend.position = "none",
            axis.text.y = element_text(size=20,colour = "black"),
            axis.text.x = element_text(size=20,colour = "black"),
            axis.title.x=element_text(size=20),
            axis.title.y=element_text(size=20),
            title = element_blank())
    ggsave(filename = paste(folder_path, subfolder, "/",
                            as.character(j), ".png",sep = ""),
           width = 8, height = 4, dpi = 300)
  }
}

# plot item response function
# probs = getprobs_gpirt(xs[idx], irfs, matrix(thresholds,nrow=1))
# q = ggplot(probs, aes(x=xs, y=p, group=order, color=factor(order))) +
#   geom_line(size=2) +ggtitle(paste("True GP IRT q",l, sep="")) +
#   theme(plot.title = element_text(hjust = 0.5))

# save estimated scores
dynamic_score_data = data.frame(matrix(ncol = 3, nrow = 0))
colnames(dynamic_score_data) <- c("wave","score", "wustlid")
for(h in 1:horizon){
  tmp = -pred_theta[,h] # negate x
  tmp = data.frame(tmp)
  colnames(tmp) = c("score")
  tmp$wave = h
  tmp$WUSTLID = all_ids
  dynamic_score_data= rbind(dynamic_score_data, tmp)
}

write.csv(dynamic_score_data, file="../results/gpirt_TAPS2014_dynamic.csv")

# train/test statistics
# train
train_lls = c()
train_acc = c()
train_response = c()
train_prediction = c()

for (h in 1:horizon) {
  wave = unique(data$wave)[h]
  for (j in 1:m) {
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
    gpirt_iccs[,j,h] = tmp$icc
    
    # test/train statistic
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

results = data.frame(train_lls,
                     train_acc,
                     train_response,
                     train_prediction)

write.csv(results, "./gpirt_TAPS_2014_train.csv")

save.image(file='gpirt_TAPS_2014.RData')