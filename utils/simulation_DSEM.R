# install.packages("dynr")
args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

if (length(args)==0) {
  SEED = 24
  C = 5
  n = 100
  m = 10
  horizon = 10
  TYPE = "DSEM"
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
library("lavaan")
set.seed(SEED)

# load data
source("getprob_gpirt.R")
HYP = paste(DATA_TYPE, "_C_", C, '_n_', n, '_m_', m, '_h_', horizon,'_CSTIRF_', CONSTANT_IRF , '_SEED_', SEED, sep="")
print(HYP)
load(file=paste("./data/", HYP, ".RData" , sep=""))

# build n*horizon*m 3d matrix
DSEM_data = aperm(data, c(1,3,2))
train_data = aperm(data_train, c(1,3,2))

# reshape to 2d matrix
C = length(unique(array(DSEM_data, c(n*horizon*m))))
dim(DSEM_data) = c(n*horizon,m)
dim(train_data) = c(n*horizon,m)

model_data = as.data.frame(train_data)
colnames(model_data) = unlist(lapply(1:m,function(i) paste("y",as.character(i), sep="")))

myModel <- ' 
 # latent variables 
   l1 =~ y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8 + y9 + y10
'
fit <- sem(model = myModel, 
           data = model_data) 
correlation_matrix = matrix(0, nrow = m, ncol=m)
for (i in 1:m){
  for (j in 1:m){
    correlation_matrix[i,j] = fitted(fit)$cov[i,j]
    correlation_matrix[j,i] = fitted(fit)$cov[j,i]
  }
}

loadings = parameterEstimates(fit)[1:m,"est"]

log_lik = fitMeasures(fit)[["logl"]]
BIC = BIC(fit)
pred_theta = predict(fit, newdata = model_data)
dim(pred_theta)=c(n,horizon)

pred_theta_sd = matrix(0.1, nrow=n, ncol = horizon)
pred_theta_ll = matrix(0, nrow=n, ncol=horizon)

for(i in 1:n){
  for (h in 1:horizon) {
    pred_theta_ll[i,h] = log(dnorm(theta[i,h], 
                                   mean=pred_theta[i,h], 
                                   sd=pred_theta_sd[i,h]))
  }
}

pred_y = matrix(0, nrow = n*horizon, ncol=m)
dim(pred_theta)=c(n*horizon)
for (i in 1:(n*horizon)){
  pred_y[i,] = rep(pred_theta[i], each = m) * loadings
}
dim(pred_theta)=c(n,horizon)

pred_y = (pred_y-min(pred_y)+1)/(max(pred_y)-min(pred_y))*(C-1)
pred_y = round(pred_y, digits = 0)

train_acc = mean(pred_y==model_data)
train_lls = log_lik / n / m / horizon
pred_acc = NULL
pred_lls = NULL

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
    IRFs = matrix(0, nrow=1, ncol=length(idx))
    for(iter in 1:length(idx)){
      tmp = rep(xs[idx][iter], each = m) * loadings
      IRFs[1,iter] = tmp[ceiling(iter/length(idx)*(C-1))]
    }
    probs = getprobs_gpirt(xs[idx], t(IRFs), matrix(thresholds,nrow=1))
    
    tmp = probs %>% 
      group_by(xs) %>%
      summarize(icc=sum(order*p))
    gpirt_iccs[,j,h] = tmp$icc
    cor_icc[j,h] = cor(gpirt_iccs[,j,h], true_iccs[,j,h])
    rmse_icc[j,h] = sqrt(mean((gpirt_iccs[,j,h]-true_iccs[,j,h])^2))
  }
}

HYP = paste(TYPE, "_C_", C, '_n_', n, '_m_', m, '_h_', horizon,'_CSTIRF_', CONSTANT_IRF , '_SEED_', SEED, sep="")

save(theta, pred_theta,pred_theta_ll,pred_theta_sd,
     train_lls,train_acc, pred_lls, pred_acc,
     cor_icc, rmse_icc, gpirt_iccs, true_iccs,
     file=paste("./results/gpirt_", HYP, ".RData" , sep=""))
