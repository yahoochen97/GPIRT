#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

R_path="~/R/x86_64-redhat-linux-gnu-library/4.0"
.libPaths(R_path)
library(rstan)

if (length(args)==0) {
  TRAIN_START_YEAR = 1981
  TRAIN_END_YEAR = 1990
  TEST_YEAR = 1995
}

if (length(args)==3){
  TRAIN_START_YEAR = as.integer(args[1])
  TRAIN_END_YEAR = as.integer(args[2])
  TEST_YEAR = as.integer(args[3])
}

library(gpirt)
library(dplyr)
library(ggplot2)
library(stats)
# gpirt_path = "~/Documents/Github/OrdGPIRT"
# setwd(gpirt_path)
TYPE = "GP"

data = read.csv("./data/dataverse_files/HumanRightsProtectionScores_v4.01.csv")

DOIRT_theta = read.csv("./data/dataverse_files/DOIRT_PhysintLatentVariableEsimates_1981_2010.csv")

# ordered variable for disappearances/extra-judical/political imprisonment/torture
# in the CIRI dataset
data = data[data$YEAR>=1981 & data$YEAR<=2010,
            c("YEAR", "CIRI", "DISAP", "KILL", 
              "POLPRIS", "TORT", "theta_mean",
              "theta_sd", "country_name")]

unique_sessions = unique(data$YEAR)
horizon = length(unique_sessions)
m = 4 # 4 items

# impute missing value in CIRI by country_name if possible
country_names = unique(data$country_name)
for(country_name in country_names){
  tmp = data[data$country_name==country_name & !is.na(data$CIRI), c("YEAR","CIRI", "country_name")]
  if(nrow(tmp)){
    id = unique(tmp$CIRI)
    tmp2 = data[data$country_name==country_name & is.na(data$CIRI), ]
    if(nrow(tmp2)){
      data[data$country_name==country_name & is.na(data$CIRI), "CIRI"] = id
    }
  }else{
    print(country_name)
  }
}

# remove  German Democratic Republic/Yemen People's Republic/Kosovo/South Sudan
data = data[!is.na(data$CIRI), ]

# Yugoslavia has 4 CIRIs: 689(1981-1991) 563(1992-1999, 2003-2005) 692(2000-2002) 560(2006-2011)
# Czechia has 2 CIRIs: 242(1981-1992) and 239(1993-2011)
# Russia has 2 CIRIs: 590(1981-1991) and 530(1992-2011)
# CIRI 299 represents German Federal Republic(1981-1990) and Germany(1990-2011)

# combine Yugoslavia
data[data$country_name=="Yugoslavia", "CIRI"] = 689

# combine Czechia
data[data$CIRI==239,"CIRI"] = 242

# combine Russian
data[data$CIRI==530, "CIRI"] = 590

# combine Germany
data = data[data$YEAR!=1990 | data$country_name!="German Federal Republic", ]
data[data$country_name=="German Federal Republic", "country_name"] = "Germany"

country_names = unique(data$country_name)
CIRIs = unique(data$CIRI)
n = length(CIRIs)

CIRI_data = array(array(NA, n*m*horizon), c(n, m, horizon))
C = 3

questions = c("DISAP", "KILL", "POLPRIS", "TORT")
for (h in 1:length(unique(data$YEAR))) {
  year = unique(data$YEAR)[h]
  for(j in 1:m){
    for(i in 1:length(CIRIs)){
      tmp = data[data$YEAR==year & data$CIRI==CIRIs[i], c("YEAR", "country_name","theta_mean","theta_sd", questions[j])]
      if(nrow(tmp)){
        CIRI_data[i,j,h] = tmp[[questions[j]]] + 1
      }
    }
  }
}

CIRI_theta = array(array(NA, n*horizon), c(n, horizon))
for (h in 1:length(unique(data$YEAR))) {
  year = unique(data$YEAR)[h]
  for(j in 1:m){
    for(i in 1:length(CIRIs)){
      tmp = DOIRT_theta[DOIRT_theta$YEAR==year & DOIRT_theta$CIRI==CIRIs[i],
                        c("YEAR", "CTRY","latentmean","latentsd", questions[j])]
      if(nrow(tmp)){
        CIRI_theta[i,h] = tmp[["latentmean"]]
      }
    }
  }
}

if(TYPE=="GP"){
  theta_os = 1
  theta_ls = 6
}else if(TYPE=="CST"){
  theta_os = 1
  theta_ls = 10*horizon
}else{
  theta_os = 1
  theta_ls = 0.1
}


theta_init = array(array(0, n*horizon), c(n, horizon))
theta_init[!is.na(CIRI_theta)] = CIRI_theta[!is.na(CIRI_theta)]
for(h in 1:horizon){
  theta_init[,h] = theta_init[,h] + 0.1*rnorm(n)
}

theta_init = theta_init[,(TRAIN_START_YEAR-1980):(TEST_YEAR-1980)]

CIRI_data_train = CIRI_data
CIRI_data_train[,,(TRAIN_END_YEAR-1979):(TEST_YEAR-1980)] = NA
CIRI_data_train = CIRI_data_train[,,(TRAIN_START_YEAR-1980):(TEST_YEAR-1980)]

SEED = 1
SAMPLE_ITERS = 500
BURNOUT_ITERS = 500
THIN = 4
CHAIN = 1

CIRI_data_train[is.na(CIRI_data_train)] = 0
stan_data <- list(n=n,
                  m=m,
                  horizon=ncol(CIRI_data_train[1,,]),
                  K=C,
                  sigma=0.1,
                  y=CIRI_data_train)
CIRI_data_train[CIRI_data_train==0] = NA

# train stan model
fit <- stan(file = "doirt-synthetic.stan",
            data = stan_data, 
            warmup = BURNOUT_ITERS, 
            iter = BURNOUT_ITERS + SAMPLE_ITERS, 
            chains = CHAIN, 
            cores = 1, 
            thin = THIN,
            control=list(adapt_delta=.98, max_treedepth = 15),
            seed = SEED,
            refresh=1
)

fit_params <- as.data.frame(fit)

SAMPLE_ITERS = SAMPLE_ITERS / THIN

samples = list()
samples[["theta"]] = array(array(0, SAMPLE_ITERS*n*ncol(CIRI_data_train[1,,])), 
                           c(SAMPLE_ITERS,n, ncol(CIRI_data_train[1,,])))
samples[["threshold"]] = array(array(0, SAMPLE_ITERS*(C+1)), 
                               c(SAMPLE_ITERS,(C+1)))
xs = seq(-5,5,0.01)
idx = 1:1001
for (it in 1:SAMPLE_ITERS) {
  samples[["fstar"]][[it]] =  array(array(0, length(xs[idx])*m*ncol(CIRI_data_train[1,,])), 
                                    c(length(xs[idx]),m, ncol(CIRI_data_train[1,,])))
  for(j in 1:m){
    for(h in 1:ncol(CIRI_data_train[1,,])){
      samples[["fstar"]][[it]][,j,h] = xs[idx]*fit_params[[paste("beta[",j,",",h,"]",sep="")]][it]
    }
  }
  for(i in 1:n){
    for(h in 1:ncol(CIRI_data_train[1,,])){
      samples[["theta"]][it,i,h] = fit_params[[paste("theta[",i,",",h,"]",sep="")]][it]
    }
  }
}

# train/test statistics
# train
train_lls = c()
train_acc = c()
train_response = c()
train_prediction = c()

# test
test_lls = c()
test_acc = c()
test_response = c()
test_prediction = c()

source("getprob_gpirt.R")

for(h in 1:ncol(CIRI_data_train[1,,])) {
  for (j in 1:m){
    IRFs = matrix(0, nrow=SAMPLE_ITERS, ncol=length(xs))
    for(iter in 1:SAMPLE_ITERS){
      IRFs[iter, ] = samples$fstar[[iter]][, j, h]
    }
    probs = getprobs_gpirt(xs, t(IRFs), samples$threshold)
    for (i in 1:n) {
      if(h<(TRAIN_END_YEAR-1979) & !is.na(CIRI_data[[i,j,h]]) & !is.na(CIRI_theta[i,h]) ){
        # train
        pred_idx = 1+as.integer((pred_theta[i,h]+5)*100)
        ll = log(probs$p[probs$xs==xs[pred_idx]])
        y_pred = which.max(ll)
        train_acc = c(train_acc, y_pred==(CIRI_data[[i,j, h]]))
        train_lls = c(train_lls, ll[CIRI_data[[i,j, h]]])
        train_response = c(train_response, CIRI_data[[i,j, h]])
        train_prediction = c(train_prediction, y_pred)
      }
      if(h==(TEST_YEAR-1980) & !is.na(CIRI_data[[i,j,h]]) & !is.na(CIRI_theta[i,h]) ){
        # test
        pred_idx = 1+as.integer((CIRI_theta[i,h]+5)*100)
        ll = log(probs$p[probs$xs==xs[pred_idx]])
        y_pred = which.max(ll)
        test_acc = c(test_acc, y_pred==(CIRI_data[[i,j, h]]))
        test_lls = c(test_lls, ll[CIRI_data[[i,j, h]]])
        test_response = c(test_response, CIRI_data[[i,j, h]])
        test_prediction = c(test_prediction, y_pred)
      }
    }
  }
}

x = TRAIN_END_YEAR - TRAIN_START_YEAR + 1
i = TEST_YEAR - TRAIN_END_YEAR
HYP = paste("_year_", TRAIN_START_YEAR, '_x_', x, '_i_', i, sep="")

# save.image(file=paste("./results/doirt_CIRI", HYP, ".RData" , sep=""))