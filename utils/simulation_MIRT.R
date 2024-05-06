args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

# gpirt_path = "~/Documents/Github/OrdGPIRT"
# setwd(gpirt_path)
TYPEs = c("graded_uni", "gpcm_uni", "sequential_uni" ,"ggum_uni")

if (length(args)==0) {
  SEED = 1
  C = 2
  n = 100
  m = 10
  horizon = 10
  TYPE = "graded_uni"
  CONSTANT_IRF = 0
  DATA_TYPE = "GP"
}

if (length(args)==8){
  SEED = as.integer(args[1])
  C = as.integer(args[2])
  n = as.integer(args[3])
  m = as.integer(args[4])
  horizon = as.integer(args[5])
  TYPE = args[6]
  CONSTANT_IRF = as.integer(args[7])
  DATA_TYPE = args[8]
}

R_path="~/R/x86_64-redhat-linux-gnu-library/4.0"
.libPaths(R_path)
library(mirt)
library(dplyr)
set.seed(SEED)

# load data
source("getprob_gpirt.R")
HYP = paste(DATA_TYPE, "_C_", C, '_n_', n, '_m_', m, '_h_', horizon,'_CSTIRF_', CONSTANT_IRF , '_SEED_', SEED, sep="")
print(HYP)
load(file=paste("./data/", HYP, ".RData" , sep=""))

# build n*horizon*m 3d matrix
mirt_data = aperm(data, c(1,3,2))
train_data = aperm(data_train, c(1,3,2))

# reshape to 2d matrix
C = length(unique(array(mirt_data, c(n*horizon*m))))
dim(mirt_data) = c(n*horizon,m)
dim(train_data) = c(n*horizon,m)

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
mirt_fit <- mirt(data = data.frame(train_data), 
                   model = factor_model,
                   itemtype = MODEL_NAME,
                   method = EM_method,
                   optimizer = "nlminb",
                   verbose = FALSE)

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
pred_theta_sd = matrix(0.1, nrow=n, ncol = horizon)
pred_theta_ll = matrix(0, nrow=n, ncol=horizon)

for(i in 1:n){
  for (h in 1:horizon) {
    pred_theta_ll[i,h] = log(dnorm(theta[i,h], 
                                   mean=pred_theta[i,h], 
                                   sd=pred_theta_sd[i,h]))
  }
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

# predict test observations and likelihood
train_acc = c()
train_ll = c()
test_acc = c()
test_ll = c()
dim(train_data) = c(n, horizon, m)
for(i in 1:n){
  for(j in 1:m){
    for(h in 1:horizon){
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
      if(is.na(data_train[i,j,h])){
        test_acc = c(test_acc, pred_y==data[i,j,h])
        if(ps[data[i,j,h]]<=-1e-6){
          print(paste(i,"_",h,"_",j,sep=""))
        }
        test_ll = c(test_ll, log(1e-6+ps[data[i,j,h]]))
      }else{
        train_acc = c(train_acc, pred_y==data[i,j,h])
        if(ps[data[i,j,h]]<=-1e-6){
          print(paste(i,"_",h,"_",j,sep=""))
        }
        train_ll = c(train_ll, log(1e-6+ps[data[i,j,h]]))
      }
    }
  }
}

train_acc = mean(train_acc)
train_lls = mean(train_ll[!is.na(train_ll)])
pred_acc = mean(test_acc)
pred_lls = mean(test_ll[!is.na(test_ll)])


# get cor of icc
xs = seq(-5,5,0.01)
idx = 301:701
gpirt_iccs = array(array(0, length(xs[idx])*m*horizon), 
                   c(length(xs[idx]),m, horizon))
true_iccs = array(array(0, length(xs[idx])*m*horizon), 
                  c(length(xs[idx]),m, horizon))
cor_icc = matrix(0, nrow=m, ncol=horizon)
rmse_icc = matrix(0, nrow=m, ncol=horizon)

for (h in 1:horizon) {
  for (j in 1:m) {
    source("true_irf.R")
    probs = getprobs_gpirt(xs[idx], irfs, matrix(thresholds,nrow=1))
    tmp = probs %>% 
      group_by(xs) %>%
      summarize(icc=sum(order*p))
    true_iccs[,j,h] = tmp$icc
    if (MODEL_NAME=="ggum"){
      for(iter in 1:length(idx)){
        A = rep(0, C)
        for(w in 1:C){
          tmp = c(0)
          if (w>=2){ tmp = c(tmp, coefs[j,3:(w+1)])}
          A[w] = A[w] + get_latent_f(coefs[j],w*(xs[idx][iter]-coefs[j,2]),tmp)
          A[w] = A[w] + get_latent_f(coefs[j],(2*C-1-w)*(xs[idx][iter]-coefs[j,2]),tmp)
        }
        ps = A/sum(A)
        gpirt_iccs[iter,j,h] = sum(ps*(1:C))
      }
    } else{
      IRFs = matrix(0, nrow=1, ncol=length(idx))
      for(iter in 1:length(idx)){
        tmp = get_latent_f(coefs[j],xs[idx][iter],coefs[j,2:C])
        IRFs[1,iter] = tmp[ceiling(iter/length(idx)*(C-1))]
      }
      probs = getprobs_gpirt(xs[idx], t(IRFs), matrix(thresholds,nrow=1))
      
      tmp = probs %>% 
        group_by(xs) %>%
        summarize(icc=sum(order*p))
      gpirt_iccs[,j,h] = tmp$icc
    }
    cor_icc[j,h] = cor(gpirt_iccs[,j,h], true_iccs[,j,h])
    rmse_icc[j,h] = sqrt(mean((gpirt_iccs[,j,h]-true_iccs[,j,h])^2))
  }
}

HYP = paste(MODEL_NAME, "_C_", C, '_n_', n, '_m_', m, '_h_', horizon,'_CSTIRF_', CONSTANT_IRF , '_SEED_', SEED, sep="")

save(theta, pred_theta,pred_theta_ll,pred_theta_sd,
     train_lls,train_acc, pred_lls, pred_acc,
     cor_icc, rmse_icc, gpirt_iccs, true_iccs,
     file=paste("./results/gpirt_", HYP, ".RData" , sep=""))

# quit()