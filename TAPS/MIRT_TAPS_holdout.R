#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)


TYPEs = c("graded_uni", "gpcm_uni", "sequential_uni" ,"ggum_uni")

if (length(args)==0) {
  TRAIN_START_YEAR = 1
  TRAIN_END_YEAR = 31
  TEST_YEAR = 41
  DROP_RATIO = 10
  SEED = 1
  TYPE = "graded_uni"
}


if (length(args)==6){
  TRAIN_START_YEAR = as.integer(args[1])
  TRAIN_END_YEAR = as.integer(args[2])
  TEST_YEAR = as.integer(args[3])
  SEED = as.integer(args[4])
  DROP_RATIO = as.integer(args[5])
  TYPE = args[6]
}

R_path="~/R/x86_64-redhat-linux-gnu-library/4.0"
.libPaths(R_path)

library(mirt)
library(dplyr)
library(ggplot2)
library(stats)
library(ltm)

print("loading data...")
source("load_TAPS.R")
C = 5

gpirt_data_train = gpirt_data

# set random seed
set.seed(SEED)

print("setting up training data...")
# random drop one observation for every case from TRAIN_END_YEAR+1 to TEST_YEAR
for(h in (TRAIN_END_YEAR+1):(TEST_YEAR)){
  for (j in 1:m){
    mask = !is.na(gpirt_data[,j,h])
    if (sum(mask==1)>0){
      drop_unit = sample(which(mask==1), as.integer(DROP_RATIO/100*length(mask)))
    }
    gpirt_data_train[drop_unit, j, h] = NA
  }
}

# build n*horizon*m 3d matrix
NA_MASK = is.na(gpirt_data_train)
gpirt_data_train[NA_MASK] = 2 
DSEM_data = aperm(gpirt_data_train, c(1,3,2))

# reshape to 2d matrix
C = length(unique(array(DSEM_data, c(n*horizon*m))))
dim(DSEM_data) = c(n*horizon,m)

# define mirt model
RANK = 1
factor_strings = c()
for (r in 1:RANK){
  factor_strings = c(factor_strings, paste('F',r, ' = ', m/RANK*(r-1)+1, '-', m/RANK*r, sep=''))
}

s = paste(factor_strings, collapse="\n")
factor_model <- mirt.model(s)

MODEL_NAME = unlist(strsplit(TYPE, "_"))[1]
UNI = unlist(strsplit(TYPE, "_"))[2]
EM_method = "QMCEM"
if(UNI=="uni"){
  factor_model = 1
  EM_method="EM"
}

# fit mirt model
mirt_fit <- mirt(data = data.frame(DSEM_data), 
                 model = factor_model,
                 itemtype = MODEL_NAME,
                 method = EM_method,
                 optimizer = "nlminb",
                 verbose = FALSE)

gpirt_data_train[NA_MASK] = NA

if(MODEL_NAME=="sequential"){
  coefs = coef(mirt_fit, simplify = TRUE)$items
} else if(MODEL_NAME=="ggum"){
  coefs = coef(mirt_fit, simplify = TRUE)$items
} else {
  coefs = coef(mirt_fit, IRTpars = TRUE, simplify = TRUE)$items
}

if(UNI=="uni"){
  loadings = matrix(as.vector(coefs[,1]))
} else{
  loadings = matrix(as.vector(coefs[,1:RANK]), nrow=m)
}

correlation_matrix = loadings %*% t(loadings)
log_lik = mirt_fit@Fit$logLik
BIC = mirt_fit@Fit$BIC
if(UNI=="uni"){
  pred_theta = array(fscores(mirt_fit), c(n,horizon))
} else{
  pred_theta = array(as.vector(fscores(mirt_fit)), c(n,horizon, RANK))
}

get_latent_f = function(as, theta, bs){
  # set na as very extreme number
  if(MODEL_NAME=="sequential"){
    bs[is.na(bs)] = -1000
  } else{
    if(as>0){
      bs[is.na(bs)] = 1000
    } else{
      bs[is.na(bs)] = -1000
    }
  }
  
  # compute latent f
  if(MODEL_NAME=="graded"){
    f = as*(theta-bs)
  } else if(MODEL_NAME=="gpcm"){
    f = as*(theta-bs)
    f = cumsum(f)
  } else if(MODEL_NAME=="sequential"){
    f = as*theta-bs
  } else if(MODEL_NAME=="ggum"){
    f = exp(as*(theta-sum(bs)))
  }
  return(f)
}

dim(DSEM_data) = c(n, horizon, m)

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

source("../getprob_gpirt.R")

for(h in (TRAIN_END_YEAR+1):horizon) {
  wave = unique(data$wave)[h]
  h_ = h-TRAIN_END_YEAR
  for (j in 1:m){
    for (i in 1:n) {
      tmp = get_latent_f(coefs[j],pred_theta[i,h],coefs[j,2:C])
      if( MODEL_NAME=="graded"){
        if(C==2){
          ps = c(1-pnorm(tmp), pnorm(tmp))
        } else {
          tmp = exp(tmp)/sum(exp(tmp))
          ps = c(1-tmp[1])
          for(c in 1:(C-2)){
            ps = c(ps, tmp[c]-tmp[c+1])
          }
          ps = as.vector(c(ps, tmp[C-1]))
        }
      } else if (MODEL_NAME=="gpcm"){
        tmp = c(0, tmp)
        ps = as.vector(exp(tmp)/sum(exp(tmp)))
      } else if (MODEL_NAME=="sequential"){
        tmp = plogis(tmp)
        ps = c(tmp,1)*c(1,cumprod(1-tmp))
      } else if (MODEL_NAME=="ggum"){
        A = rep(0, C)
        for(w in 1:C){
          tmp = c(0)
          if (w>=2){ tmp = c(tmp, coefs[j,3:(w+1)])}
          A[w] = A[w] + get_latent_f(coefs[j],w*(pred_theta[i,h]-coefs[j,2]),tmp)
          A[w] = A[w] + get_latent_f(coefs[j],(2*C-1-w)*(pred_theta[i,h]-coefs[j,2]),tmp)
        }
        ps = A/sum(A)
      }
      pred_y = which.max(ps)
      
      if(!is.na(gpirt_data_train[i,j,h]) & !is.na(gpirt_data[i,j,h])){
        # train
        train_lls = c(train_lls, log(ps[gpirt_data[i,j, h]]))
        train_response = c(train_response, gpirt_data[i,j, h])
        train_prediction = c(train_prediction, pred_y)
      }
      if(is.na(gpirt_data_train[i,j,h]) & !is.na(gpirt_data[i,j,h])){
        # test
        test_lls[[h_]] = c(test_lls[[h_]],  log(ps[gpirt_data[i,j, h]]))
        test_response[[h_]] = c(test_response[[h_]],gpirt_data[i,j, h])
        test_prediction[[h_]] = c(test_prediction[[h_]], pred_y)
      }
    }
  }
  test_acc[[h_]] = test_prediction[[h_]]==test_response[[h_]]
}

train_acc = train_prediction==train_response

file_name = paste("./results/", MODEL_NAME, "_TAPS_holdout_", "DR_", DROP_RATIO, "_SEED_", SEED, ".RData" , sep="")

save(gpirt_data_train, gpirt_data, pred_theta,train_lls,
     train_acc, train_response, train_prediction,test_lls,
     test_acc, test_response, test_prediction,
     file=file_name)