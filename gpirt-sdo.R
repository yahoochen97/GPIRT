#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
    SEED = 1234
    C = 2
    n = 50
    m = 20
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

gpirt_path = "~/Documents/GitHub/gpirtr"
setwd(gpirt_path)
remove.packages("gpirt")
library(Rcpp)
Rcpp::compileAttributes()
install.packages(gpirt_path, type="source", repos = NULL)
library(gpirt)
library(catSurv)
library(ggplot2)

setwd("~/Documents/GitHub/GPIRT")

data = data.matrix(SDO)[1:500,]
thresholds <- c(-Inf)
unique_ys = unique(as.vector(data))
C = length(unique(unique_ys[!is.na(unique_ys)]))
for(i in 1:(C-1)){
    thresholds = c(thresholds, qnorm(i/C, 0, 1, 1, 0))
}
thresholds = c(thresholds, Inf)

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
# set.seed(SEED)
# theta <- seq(-2, 2, length.out = n) # Respondent ability parameters
# alpha <- seq(-1, 1, length.out = m) # Item difficulty parameters
# beta  <- seq(-3, 3, length.out = m) # Item discrimination parameters
# data <- gen_responses(theta, alpha, beta, thresholds)

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

samples <- gpirtMCMC(data_train, 100,100,
                     beta_prior_sds = matrix(0.5, nrow = 2, ncol = ncol(data_train)),
                     vote_codes = NULL, thresholds=thresholds)
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
            f_pred = samples$IRFs[as.integer((pred_theta[i]+5)*100), j]
            f_pred = log(f_pred/(1-f_pred))
            ll = ordinal_lls(f_pred, thresholds)
            y_pred = which.max(ll)
            if(ll[data[[i,j]]]<(-8.8)){
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
  C = length(thresholds) - 1
  probs = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(probs) = c("order", "xs","p")
  for (i in 1:length(xs)) {
    ps = c(rep(0, C))
    for (c in 1:C){
      z1 = thresholds[c] - irfs[i]
      z2 = thresholds[c+1] - irfs[i]
      ps[c] = pnorm(z2)-pnorm(z1)
      probs[nrow(probs) + 1,] = c(c, xs[i],  ps[c])
    }
  }
  return(probs)
}

idx = 250:650
for (i in 1:m) {
  probs = getprobs_catSurv(xs[idx], sdo_cat@discrimination[[paste("q",i, sep="")]], sdo_cat@difficulty[[paste("q",i, sep="")]])
  # probs = getprobs_gpirt(xs[idx], alpha[i]+beta[i]*xs[idx],thresholds)
  p = ggplot(probs, aes(x=xs, y=p, group=order, color=factor(order))) +
    geom_line(size=2) +ggtitle(paste("True IRT q",i, sep="")) +
      theme(plot.title = element_text(hjust = 0.5))
  print(p)
  # ggsave(paste("~/Documents/GitHub/GPIRT/figures/trueirtq",i,".pdf", sep=""), plot=p, width = 7, height = 4, units = "in")
  probs2 = getprobs_gpirt(-xs[idx], log(samples$IRFs[idx,i]/(1-samples$IRFs[idx,i])),
                          thresholds)
  q = ggplot(probs2, aes(x=xs, y=p, group=order, color=factor(order))) +
    geom_line(size=2) +ggtitle(paste("GPIRT q",i, sep="")) +
      theme(plot.title = element_text(hjust = 0.5))
  print(q)
  # ggsave(paste("~/Documents/GitHub/GPIRT/figures/gpirtq",i,".pdf", sep=""), plot=q, width = 7, height = 4, units = "in")
  # ggsave("mtcars.pdf", width = 20, height = 20, units = "cm")
  # plot(xs[idx], log(samples$IRFs[idx,i]/(1-samples$IRFs[idx,i])))
  # lines(xs[idx], alpha[i]+beta[i]*xs[idx])
}

# hist(colMeans(samples$theta))
# cor(theta,colMeans(samples$theta))
# for(i in 1:16){plot(xs[idx], log(samples$IRFs[idx,i]/(1-samples$IRFs[idx,i])))}
