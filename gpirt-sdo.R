#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
    SEED = 1234
    C = 5
    n = 100
    m = 20
    PLOT = 0
}
if (length(args)==1){
    SEED = as.integer(args[1])
    C = 2
    n = 50
    m = 20
}
if (length(args)==4){
    SEED = as.integer(args[1])
    C = as.integer(args[2])
    n = as.integer(args[3])
    m = as.integer(args[4])
}

R_path="~/R/x86_64-redhat-linux-gnu-library/4.0"
.libPaths(R_path)
gpirt_path = "../gpirt"
setwd(gpirt_path)
remove.packages("gpirt")
library(Rcpp)
Rcpp::compileAttributes()
install.packages(gpirt_path, type="source", repos = NULL, lib=R_path)
library(gpirt)
library(catSurv)
library(ggplot2)

setwd("../OrdGPIRT")

thresholds <- c(-Inf)
# data = data.matrix(SDO)[1:500,]
# unique_ys = unique(as.vector(data))
# C = length(unique(unique_ys[!is.na(unique_ys)]))
for(i in 1:(C-1)){
    thresholds = c(thresholds, qnorm(i/C, 0, 1, 1, 0))
}
thresholds = c(thresholds, Inf)

thresholds = c(-Inf,-3,0,1,2,Inf)

gen_responses <- function(theta, alpha, beta, thresholds) {
  # ordinal regression
  C <- length(thresholds) - 1
  n <- length(theta)
  m <- length(alpha)
  responses <- matrix(0, n, m)
  for ( j in 1:m ) {
    for ( i in 1:n ) {
      f = alpha[j] + beta[j] * theta[i]
      ps = rep(0, C)
      for (c in 1:C) {
        z1 = thresholds[c] - f
        z2 = thresholds[c+1] -f
        ps[c] = pnorm(z2) - pnorm(z1)
      }
      responses[i, j] <- sample(1:C, 1, prob = ps)
    }
  }
  return(responses)
}
set.seed(SEED)
theta <- runif(n, -2, 2) # Respondent ability parameters
alpha <- runif(m, -2, 2) # Item difficulty parameters
beta  <- runif(m, -2, 2) # Item discrimination parameters
data <- gen_responses(theta, alpha, beta, thresholds)

# split into train and test data
N = nrow(data)*ncol(data)
na_mask = matrix(is.na(data), nrow = N)
train_idx = rep(0,N)
train_idx[sample((1:N)[!na_mask],as.integer(0.8*N),replace=FALSE)] = 1
train_idx = matrix(train_idx, nrow = nrow(data))
data_train = data

for (i in 1:nrow(data)) {
    for (j in 1:ncol(data)) {
        if(train_idx[i,j]==0){
            data_train[i,j] = NA
        }
    }
}

SAMPLE_ITERS = 100
BURNOUT_ITERS = 100
samples <- gpirtMCMC(data_train, SAMPLE_ITERS,BURNOUT_ITERS,
                     beta_prior_sds = matrix(0.5, nrow = 2, ncol = ncol(data_train)),
                     vote_codes = NULL, thresholds=NULL)
# save(samples, file = "vignettes/sdo.RData")
# load(file = "vignettes/sdo.RData")

xs = seq(-5,5,0.01)
pred_theta = colMeans(samples$theta)

ordinal_lls = function(f, thresholds){
    result = c()
    for (c in 1:(length(thresholds)-1)) {
        z1 = thresholds[c] - f;
        z2 = thresholds[c+1] - f;
        result = c(result,log(pnorm(z2)-pnorm(z1)));
    }
    return(result)
}

pred_lls = c()
pred_acc = c()
train_lls = c()
train_acc = c()
for (i in 1:nrow(data)) {
    for (j in 1:ncol(data)) {
        if(!is.na(data[[i,j]])){
          lls = matrix(0,nrow=SAMPLE_ITERS, ncol = C)
          y_pred = rep(0, SAMPLE_ITERS)
          for (iter in 1:SAMPLE_ITERS) {
            f_pred = samples$IRFs[as.integer((pred_theta[i]+5)*100), j, iter]
            ll = ordinal_lls(f_pred, samples$threshold[iter,])
            lls[iter,] = ll
            y_pred[iter] =  which.max(ll)
          }
          ll = colMeans(lls)
          y_pred = round(mean(y_pred))
            if(ll[data[[i,j]]]<(-16.6)){
                print(i)
                print(j)
            }
            if(train_idx[i,j]==0){
                pred_acc = c(pred_acc, y_pred==(data[[i,j]]))
                pred_lls = c(pred_lls, ll[data[[i,j]]])
            }else{
                train_acc = c(train_acc, y_pred==(data[[i,j]]))
                train_lls = c(train_lls, ll[data[[i,j]]])
            }
        }
    }
}

print(cor(theta,pred_theta))
print(mean(train_lls))
print(mean(train_acc))
print(mean(pred_lls))
print(mean(pred_acc))

# data("sdo_cat")
getprobs_catSurv = function(xs, discrimination, difficulties){
  C = length(difficulties) + 1
  probs = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(probs) = c("order", "xs","p")
  for (i in 1:length(xs)) {
    ps = c(rep(0, C),1)
    for (c in 1:(C-1)){
      ps[1+c] = exp(difficulties[c]-discrimination*xs[i])
      ps[1+c] = ps[1+c]/(1+ps[1+c])
    }
    for (c in 1:C){
      probs[nrow(probs) + 1,] = c(c, xs[i],  ps[c+1] - ps[c])
    }
  }
  return(probs)
}

getprobs_gpirt = function(xs, irfs, thresholds){
  C = ncol(thresholds) - 1
  probs = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(probs) = c("order", "xs","p")
  for (i in 1:length(xs)) {
    ps = matrix(0, nrow=ncol(irfs), ncol=C)
    for (c in 1:C){
      for (iter in 1:ncol(irfs)){
        z1 = thresholds[iter, c] - irfs[i,iter]
        z2 = thresholds[iter, c+1] - irfs[i,iter]
        ps[iter, c] = pnorm(z2)-pnorm(z1)
      }
    }
    ps = colMeans(ps)
    for (c in 1:C){
      probs[nrow(probs) + 1,] = c(c, xs[i],  ps[c])
    }
  }
  return(probs)
}

if (PLOT){
idx = (as.integer(min(pred_theta)*100+500)):(as.integer(max(pred_theta)*100+500))
for (j in 1:m) {
  # probs = getprobs_catSurv(xs[idx], sdo_cat@discrimination[[paste("q",j, sep="")]],
  #                           sdo_cat@difficulty[[paste("q",j, sep="")]])
  probs = getprobs_gpirt(xs[idx], matrix(alpha[j]+beta[j]*xs[idx], ncol=1),
                         matrix(thresholds,nrow=1))
  p = ggplot(probs, aes(x=xs, y=p, group=order, color=factor(order))) +
    geom_line(size=2) +ggtitle(paste("True IRT q",j, sep="")) +
      theme(plot.title = element_text(hjust = 0.5))
  print(p)
  # ggsave(paste("/figures/trueirtq",i,".pdf", sep=""), plot=p, width = 7, height = 4, units = "in")
  probs2 = getprobs_gpirt(xs[idx], samples$IRFs[idx,j,],
                          samples$threshold)
  q = ggplot(probs2, aes(x=xs, y=p, group=order, color=factor(order))) +
    geom_line(size=2) +ggtitle(paste("GPIRT q",j, sep="")) +
      theme(plot.title = element_text(hjust = 0.5))
  print(q)
  # ggsave(paste("figures/gpirtq",i,".pdf", sep=""), plot=q, width = 7, height = 4, units = "in")
  # plot(pred_theta, data[,9], pch=4, ylim=c(-2,5))
  # points(pred_theta,rowMeans(samples$f[,9,]))
  # for (c in 2:C) {
  #   abline(h = mean(samples$threshold[,c]), col="red", lwd=1, lty=1)
  # }
}
}
# hist(colMeans(samples$theta))
# cor(theta,colMeans(samples$theta))
# for(i in 1:16){plot(xs[idx], log(samples$IRFs[idx,i]/(1-samples$IRFs[idx,i])))}
