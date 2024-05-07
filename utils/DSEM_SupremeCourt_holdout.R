#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

if (length(args)==0) {
  TRAIN_START_YEAR = 2000
  TRAIN_END_YEAR = 2015
  TEST_YEAR = 2020
  SEED = 1
  TYPE = "DSEM"
}

if (length(args)==5){
  TRAIN_START_YEAR = as.integer(args[1])
  TRAIN_END_YEAR = as.integer(args[2])
  TEST_YEAR = as.integer(args[3])
  SEED = as.integer(args[4])
  TYPE = args[5]
}

R_path="~/R/x86_64-redhat-linux-gnu-library/4.0"
.libPaths(R_path)

library(dplyr)
library(ggplot2)
library(stats)
library("lavaan")

# gpirt_path = "~/Documents/Github/OrdGPIRT"
# setwd(gpirt_path)
print("loading data...")
source("supremecourt_data.R")
C = 3
gpirt_data_train = gpirt_data

# set random seed
set.seed(SEED)

print("setting up training data...")
# random drop one observation for every case from TRAIN_END_YEAR+1 to TEST_YEAR
for(h in (TRAIN_END_YEAR-2000+1+1):(TEST_YEAR-2000+1)){
  year = unique(data$term)[h]
  caseIds = as.character(unique(data[data$term==year,"caseId"]))
  for (j in 1:length(caseIds)){
    mask = !is.na(gpirt_data[,j,h])
    drop_justice = sample(which(mask==1), 1)
    gpirt_data_train[drop_justice, j, h] = NA
  }
}

# build n*horizon*m 3d matrix
NA_MASK = is.na(gpirt_data_train)
gpirt_data_train[NA_MASK] = 2 
DSEM_data = aperm(gpirt_data_train, c(1,3,2))

# reshape to 2d matrix
C = length(unique(array(DSEM_data, c(n*horizon*m))))
dim(DSEM_data) = c(n*horizon,m)

model_data = as.data.frame(DSEM_data)
colnames(model_data) = unlist(lapply(1:(m),function(i) paste("y",as.character(i), sep="")))

myModel <- '
  l1 =~ y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8 + y9 + y10 + 
 y11 + y12 + y13 + y14 + y15 + y16 + y17 + y18 + y19 + y20 +
 y21 + y22 + y23 + y24 + y25 + y26 + y27 + y28 + y29 + y30 + 
 y31 + y32 + y33 + y34 + y35 + y36 + y37 + y38 + y39 + y40 +
 y41 + y42 + y43 + y44 + y45 + y46 + y47 + y48 + y49 + y50 + 
 y51 + y52 + y53 + y54 + y55 + y56 + y57 + y58 + y59 + y60 + 
 y61 + y62 + y63 + y64 + y65 + y66 + y67 + y68 + y69 + y70 + 
 y71 + y72 + y73 + y74 + y75
'

fit <- sem(model = myModel, 
           data = model_data) 

gpirt_data_train[NA_MASK] = NA
loadings = parameterEstimates(fit)[1:m,"est"]
log_lik = fitMeasures(fit)[["logl"]]
BIC = BIC(fit)
pred_theta = predict(fit, newdata = model_data)
dim(pred_theta)=c(n,horizon)

print("analyzing results...")
# train/test statistics
# train
train_lls = c()
train_acc = c()
train_response = c()
train_prediction = c()

# test
test_lls = vector(mode='list', length=(TEST_YEAR-TRAIN_END_YEAR))
test_acc = vector(mode='list', length=(TEST_YEAR-TRAIN_END_YEAR))
test_response = vector(mode='list', length=(TEST_YEAR-TRAIN_END_YEAR))
test_prediction = vector(mode='list', length=(TEST_YEAR-TRAIN_END_YEAR))

source("getprob_gpirt.R")

for(h in (2+TRAIN_END_YEAR-TRAIN_START_YEAR):horizon) {
  h_ = h+1999-TRAIN_END_YEAR
  year = unique(data$term)[h]
  caseIds = as.character(unique(data[data$term==year,"caseId"]))
  for (j in 1:length(caseIds)){
    for (i in 1:n) {
      if(!is.na(gpirt_data_train[i,j,h]) & !is.na(gpirt_data[i,j,h])){
        # train
        y_pred = sum(pred_theta[i, h_] * loadings)
        train_lls = c(train_lls, log_lik / sum(!is.na(gpirt_data_train)))
        train_response = c(train_response, gpirt_data[i,j, h])
        train_prediction = c(train_prediction, y_pred)
      }
      if(is.na(gpirt_data_train[i,j,h]) & !is.na(gpirt_data[i,j,h])){
        # test
        y_pred = sum(pred_theta[i, h_] * loadings)
        test_response[[h_]] = c(test_response[[h_]],gpirt_data[i,j, h])
        test_prediction[[h_]] = c(test_prediction[[h_]], y_pred)
      }
    }
  }
  test_prediction[[h_]]  = 1 + (test_prediction[[h_]]-min(test_prediction[[h_]] )+1)/(max(test_prediction[[h_]] )-min(test_prediction[[h_]] ))*(C-1)
  tmp = pnorm(abs(test_prediction[[h_]]-test_response[[h_]]))
  test_lls[[h_]] = log((1-tmp)/(tmp))
  test_prediction[[h_]] = round(test_prediction[[h_]], digits = 0)
  test_acc[[h_]] = test_prediction[[h_]]==test_response[[h_]]
}

train_prediction  = 1 +(train_prediction -min(train_prediction )+1)/(max(train_prediction )-min(train_prediction ))*(C-1)
train_prediction = round(train_prediction, digits = 0)

train_acc = train_prediction==train_response


if(TYPE=="DSEM"){
  file_name = paste("./results/DSEM_SupremeCourt_holdout_SEED_", SEED, ".RData" , sep="")
}

save(gpirt_data_train, gpirt_data, pred_theta,pred_theta_sd,train_lls,
     train_acc, train_response, train_prediction,test_lls,
     test_acc, test_response, test_prediction,
     file=file_name)
