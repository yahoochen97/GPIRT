R_path="~/R/x86_64-redhat-linux-gnu-library/4.0"
.libPaths(R_path)

library(gpirt)
library(dplyr)
library(ggplot2)
library(stats)

# gpirt_path = "~/Documents/Github/OrdGPIRT/VAA"
# setwd(gpirt_path)
SAMPLE_ITERS = 500
BURNOUT_ITERS = 500
THIN = 4
CHAIN = 1
TYPE = "GP"
SEED=12345
set.seed(SEED)

# load data
data = readRDS("VAA_Data.rds")
item_texts = data$item_text
item_types = data$item_types
data = data$Finland_Dataset

# drop nominal and open-end questions
drop_items = item_types %in% c("Categorical","Open-Ended" )
item_types = item_types[!drop_items]
item_texts = item_texts[!drop_items]

years = unique(data$year)
horizon = length(years)

# iterate items
item_names = names(item_texts)
MIN_ITEMS = 1000 # 1000 # filter items with too few responses
num_response_items = rep(0, length(item_names))
for(i in 1:length(item_names)){
  num_response_items[i] = sum(!is.na(data[,item_names[i]]))
}
item_names = item_names[num_response_items>=MIN_ITEMS]

# iterate years
num_item_years = rep(0, horizon)
for(h in 1:horizon){
  mask = unlist(lapply(item_names, FUN=function(x) 
        substr(x,2,3)==substr(toString(years[h]),3,4)))
  num_item_years[h] = sum(mask)
}

# iterate respondent
cids = levels(data$c_id)
num_item_cids = rep(0, length(cids))
for(i in 1:length(cids)){
  num_item_cids[i] = sum(!is.na(data[data$c_id==cids[i],item_names]))
} 


MIN_CIDS = 100 # 100 # filter respondents with too few responses
n = sum(num_item_cids>=MIN_CIDS)
cids = cids[num_item_cids>=MIN_CIDS]
m = max(num_item_years)
C = 5

# iterate responses
unique_responses = c()
response_levels = list()
for(j in 1:length(item_names)){
  tmp = pull(unique(data[data$c_id %in% cids,item_names[j]]))
  response_levels[[item_names[j]]] = tmp
  unique_responses = unique(c(unique_responses, tmp))
}

# build response matrix
responses = array(rep(NA,n*m*horizon),c(n,m,horizon))
for(h in 1:horizon){
  mask = unlist(lapply(item_names, FUN=function(x) 
    substr(x,2,3)==substr(toString(years[h]),3,4)))
  for(i in 1:n){
    for(j in 1:num_item_years[h]){
      tmp = data[data$c_id==cids[i] & data$year==years[h], item_names[mask][j]]
      if(nrow(tmp)==0){next}
      tmp = pull(tmp)[[1]]
      if(is.na(tmp)){next}
      if (tmp=="Totally Agree" | tmp=="Strongly agree" | tmp == "Yes" | tmp=="Agree" | tmp=="Mentioned"){
        responses[i,j,h] = 5
      }
      if (tmp=="Totally Disagree" | tmp=="Strongly disagree" | tmp == "No" | tmp=="Disagree" | tmp=="Not mentioned"){
        responses[i,j,h] = 1
      }
      if (tmp=="Agree to some extent" | tmp=="Somewhat Agree"){
        responses[i,j,h] = 4
      }
      if (tmp=="Disagree to some extent" | tmp=="Somewhat Disagree"){
        responses[i,j,h] = 2
      }
      if (tmp=="I Can't Say" | tmp=="Can't say" | tmp=="Cannot say" | tmp=="Pass" | tmp=="Caused equally by the Finns and immigrants themselves"){
        responses[i,j,h] = 3
      }
      if (tmp=="Definitely more by the Finns"){
        responses[i,j,h] = 5
      }
      if (tmp=="Definitely more by immigrants themselves"){
        responses[i,j,h] = 1
      }
    }
  }
}

data = responses

if(TYPE=="GP"){
  theta_os = 1
  theta_ls = 2
}else if(TYPE=="CST"){
  theta_os = 0
  theta_ls = -1
}else{
  theta_os = 0
  theta_ls = 1
}


beta_prior_means = matrix(0, nrow = 3, ncol = m)
beta_prior_sds =  matrix(1.0, nrow = 3, ncol = m)
beta_prior_sds[3,] = 0
theta_prior_means = matrix(0, nrow = 2, ncol = n)
theta_prior_sds =  matrix(0.0, nrow = 2, ncol = n)
theta_init = matrix(0, nrow = n, ncol = horizon)
theta_init[,1] = rnorm(n)
for (h in 2:horizon) {
 theta_init[,h] = theta_init[,1]
}

all_samples <- gpirtMCMC(data, SAMPLE_ITERS,BURNOUT_ITERS,
                         THIN=THIN, CHAIN=CHAIN, vote_codes = NULL,
                         beta_prior_means = beta_prior_means,
                         beta_prior_sds = beta_prior_sds, 
                         theta_prior_means = theta_prior_means,
                         theta_prior_sds = theta_prior_sds,
                         theta_os = theta_os, theta_ls = theta_ls, 
                         theta_init = theta_init, KERNEL = "RBF",
                         thresholds=NULL, SEED=SEED, constant_IRF = 0)


SAMPLE_ITERS = SAMPLE_ITERS/THIN
# library(rstan)
# sims <- matrix(rnorm((SAMPLE_ITERS)*CHAIN), nrow = SAMPLE_ITERS, ncol = CHAIN)
# theta_rhats = matrix(rnorm(n*horizon), nrow = n, ncol = horizon)
# for(i in 1:n){
#   for(h in 1:horizon){
#     for(c in 1:CHAIN){
#       pred_theta = colMeans(all_samples[[c]]$theta)
#       sims[,c] = all_samples[[c]]$theta[,i,h]
#     }
#     theta_rhats[i, h] = Rhat(sims)
#   }
# }

samples = all_samples[[1]]

# posterior mean and std of latent positions
xs = seq(-5,5,0.01)
pred_theta = matrix(0, nrow=nrow(data), ncol=dim(data)[3])
pred_theta_sd = matrix(0, nrow=nrow(data), ncol=dim(data)[3])
for(i in 1:n){
  for (h in 1:horizon) {
    tmp = samples$theta[-1,i,h]
    pred_theta[i,h] = mean(tmp)
    pred_theta_sd[i,h] = sd(tmp)
  }
}

# sum over linear and non-lineaer components in IRF
for(it in 1:SAMPLE_ITERS){
  for(h in 1:horizon){
    for(j in 1:num_item_years[h]){
      samples$f[[it]][,j,h] = samples$f[[it]][,j,h] + samples$beta[[it]][1,j,h] + samples$beta[[it]][2,j,h]*samples$theta[it,,h]
      samples$fstar[[it]][,j,h] = samples$fstar[[it]][,j,h] + samples$beta[[it]][1,j,h] + samples$beta[[it]][2,j,h]*xs
    }
  }
}

# posterior of IRFs
folder_path = "./results/figures/"
dir.create(file.path(folder_path), showWarnings = FALSE)
xs = seq(-5,5,0.01)
idx = 1:1001
gpirt_iccs = array(array(0, length(xs[idx])*m*horizon),
                   c(length(xs[idx]),m, horizon))

source("../getprob_gpirt.R")
for (h in 1:horizon){
  # iterate year
  mask = unlist(lapply(item_names, FUN=function(x) 
    substr(x,2,3)==substr(toString(years[h]),3,4)))
  for (j in 1:num_item_years[h]) {
    # iterate each item
    IRFs = matrix(0, nrow=SAMPLE_ITERS, ncol=length(idx))
    for(iter in 1:SAMPLE_ITERS){
      IRFs[iter, ] = samples$fstar[[iter]][idx, j, h]
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

save.image("./results/gpirt_VAA.RData")

# plot posterior of IRFs, estimated latent positions and actual responses
for (h in 1:horizon){
  # iterate year / generate subfolders
  subfolder = as.character(years[h])
  dir.create(file.path(folder_path,subfolder), showWarnings = FALSE)
  mask = unlist(lapply(item_names, FUN=function(x) 
    substr(x,2,3)==substr(toString(years[h]),3,4)))
  for (j in 1:num_item_years[h]) {
    # iterate each item
    item_name = item_names[mask][j]
    item_x = pred_theta[,h]
    item_response = data[,j,h]
    irf_plot = data.frame(item_x,item_response)
    idx = floor((min(item_x)+5)*100):ceiling((max(item_x)+5)*100)
    gpirt_plot = data.frame(xs[idx],gpirt_iccs[idx,j,h])
    colnames(gpirt_plot) = c("xs","icc")
    response_labels = response_levels[[item_name]][!is.na(response_levels[[item_name]])]
    if ("Caused equally by the Finns and immigrants themselves" %in% response_labels){
      response_labels = response_labels[response_labels!="Can't say"]
    }
    if(length(response_labels)==4 & !("Can't say" %in% response_labels)){
      response_labels = c(response_labels, "Can't say")
    }
    if (length(response_labels)==3){
      response_colors = c('black', 'blue', 'red')
      response_breaks = c(1,3,5)
    }
    if (length(response_labels)==2){
      response_colors = c('black', 'red')
      response_breaks = c(1,5)
    }
    if (length(response_labels)==5){
      response_colors = c('black','black', 'blue', 'red', 'red')
      response_breaks = 1:5
    }
    p = ggplot()+
      geom_point(data = na.omit(irf_plot), aes(x=item_x,y=item_response,color=factor(item_response)),
                 size=4, shape="|") +
      # geom_text(data = na.omit(irf_plot), 
      #           aes(x=item_x,y=item_response, label=all_justice_names),
      #           vjust=1.5*VAJUST+0.5,hjust=0, angle=ANGLE) +
      scale_color_manual(name='Level',
                         labels=response_labels,
                         values=response_colors)+
      geom_line(data = gpirt_plot, aes(x=xs,y=icc), size=1)+
      scale_x_continuous(name=bquote(theta), breaks = seq(-4, 4, by = 1)) +
      scale_y_continuous(name = "Expected Responses",
                         breaks = NULL,
                         labels = NULL,
                         limits = c(1,5),
                         sec.axis = sec_axis(~.,
                                             breaks = response_breaks,
                                             labels = response_labels))+
      theme(panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=2),
            legend.position = "none",
            axis.text.y = element_text(size=12,colour = "black"),
            axis.text.x = element_text(size=12,colour = "black"),
            axis.title.x=element_text(size=20,face="bold",colour = "black"),
            plot.title = element_text(hjust = 0.5))
      # ggtitle(caseName)
    ggsave(filename = paste(folder_path, subfolder, "/",
                            item_name, ".png",sep = ""),
           width = 8, height = 6, dpi = 300)
  }
}

