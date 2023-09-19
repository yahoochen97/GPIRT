#!/usr/bin/env Rscript

gpirt_path = "~/Documents/Github/OrdGPIRT"
setwd(gpirt_path)
setwd("./TAPS")

R_path="~/R/x86_64-redhat-linux-gnu-library/4.0"
.libPaths(R_path)
library(rstan)
rstan_options(auto_write = TRUE)
args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

if (length(args)==0) {
  TRAIN_START_YEAR = 1
  TRAIN_END_YEAR = 42
  TEST_YEAR = 41
  DROP_RATIO = 0
  SEED = 12345
}

if (length(args)==5){
  TRAIN_START_YEAR = as.integer(args[1])
  TRAIN_END_YEAR = as.integer(args[2])
  TEST_YEAR = as.integer(args[3])
  SEED = as.integer(args[4])
  DROP_RATIO = as.integer(args[5])
}

library(gpirt)
library(dplyr)
library(ggplot2)
library(stats)
library(ltm)

source("load_TAPS.R")

# set random seed
set.seed(SEED)

# get icc
C=5

# train/test statistics
# train
train_lls = c()
train_acc = c()
train_response = c()
train_prediction = c()

for (h in 1:horizon) {
  wave = unique(data$wave)[h]
  for (j in 1:m) {
    
    # test/train statistic
    for (i in 1:n) {
      if(!is.na(gpirt_data_train[i,j,h]) & !is.na(gpirt_data[i,j,h])){
        # train
        pred_idx = 1+as.integer((pred_theta[i,h]+5)*100)
        ll = log(probs$p[probs$xs==xs[pred_idx]])
        y_pred = which.max(ll)
        train_acc = c(train_acc, y_pred==(gpirt_data[i,j, h]))
        train_lls = c(train_lls, ll[gpirt_data[i,j, h]])
        train_response = c(train_response, gpirt_data[i,j, h])
        train_prediction = c(train_prediction, y_pred)
      }
      if(is.na(gpirt_data_train[i,j,h]) & !is.na(gpirt_data[i,j,h])){
        # test
        pred_idx = 1+as.integer((pred_theta[i,h]+5)*100)
        ll = log(probs$p[probs$xs==xs[pred_idx]])
        y_pred = which.max(ll)
        test_acc[[h_]] = c(test_acc[[h_]], y_pred==(gpirt_data[i,j,h]))
        test_lls[[h_]] = c(test_lls[[h_]], ll[gpirt_data[i,j, h]])
        test_response[[h_]] = c(test_response[[h_]],gpirt_data[i,j, h])
        test_prediction[[h_]] = c(test_prediction[[h_]], y_pred)
      }
    }
  }
}

results = data.frame(train_lls,
                     train_acc,
                     train_response,
                     train_prediction)

# write.csv(results, "./doirt_TAPS_2014_train.csv")
