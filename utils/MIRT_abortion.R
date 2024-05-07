# gpirt_path = "~/Documents/Github/gpirt"
# setwd(gpirt_path)
# library(Rcpp)
# Rcpp::compileAttributes()
# install.packages(gpirt_path, type="source", repos = NULL)#,lib=R_path, INSTALL_opts = '--no-lock')
# setwd("../OrdGPIRT")

library(dplyr)
library(ggplot2)
library(stats)
library(haven)
library(mirt)



gpirt_path = "~/Documents/Github/OrdGPIRT"
setwd(gpirt_path)

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

# build n*horizon*m 3d matrix
MIRT_data = aperm(rollcall_data, c(1,3,2))
NA_MASK = is.na(MIRT_data)
# DSEM_data[NA_MASK] = 2

# reshape to 2d matrix
C = length(unique(array(MIRT_data, c(n*horizon*m))))
dim(MIRT_data) = c(n*horizon,m)

# define mirt model
RANK = 1
factor_strings = c()
for (r in 1:RANK){
  factor_strings = c(factor_strings, paste('F',r, ' = ', m/RANK*(r-1)+1, '-', m/RANK*r, sep=''))
}

s = paste(factor_strings, collapse="\n")
factor_model <- mirt.model(s)

TYPEs = c("graded_uni", "ggum_uni")
TYPEs = c("gpcm_uni", "sequential_uni")

for (TYPE in TYPEs){
MODEL_NAME = unlist(strsplit(TYPE, "_"))[1]
UNI = unlist(strsplit(TYPE, "_"))[2]
EM_method = "QMCEM"
if(UNI=="uni"){
  factor_model = 1
  EM_method="EM"
}

# fit mirt model
mirt_fit <- mirt(data = data.frame(MIRT_data), 
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

loadings = matrix(as.vector(coefs[,1]))

log_lik = mirt_fit@Fit$logLik
BIC = mirt_fit@Fit$BIC
pred_theta = array(fscores(mirt_fit), c(n,horizon))

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
C = 2
train_acc = c()
train_ll = c()
dim(MIRT_data) = c(n, horizon, m)
pred_ys = array(rep(0, n*horizon*m), c(n, horizon, m))
for(i in 1:n){
  for(j in 1:m){
    for(h in 1:horizon){
      if( MODEL_NAME=="graded"){
        tmp = get_latent_f(coefs[j],pred_theta[i,h],coefs[j,1:2])
      }else{
        tmp = get_latent_f(coefs[j],pred_theta[i,h],coefs[j,2:C])
      }
      
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
          if (w>=2){ tmp = c(tmp, coefs[j,2:(w+1)])}
          A[w] = A[w] + get_latent_f(coefs[j],w*(pred_theta[i,h]-coefs[j,2]),tmp)
          A[w] = A[w] + get_latent_f(coefs[j],(2*C-1-w)*(pred_theta[i,h]-coefs[j,2]),tmp)
        }
        ps = A/sum(A)
      }
      pred_y = which.max(ps)
      if(!is.na(MIRT_data[i,h, j])){
        pred_ys[i,h,j] = pred_y
        train_acc = c(train_acc, pred_y==MIRT_data[i,h, j])
        train_ll = c(train_ll, log(1e-6+ps[MIRT_data[i,h, j]]))
      }
    }
  }
}

train_acc = mean(train_acc)
train_lls = mean(train_ll[!is.na(train_ll)])
pred_acc = mean(test_acc)
pred_lls = mean(test_ll[!is.na(test_ll)])


library(pROC)
print(TYPE)
print(train_acc)
print(train_lls)
print(auc(MIRT_data[!NA_MASK], pred_ys[!NA_MASK]))
}