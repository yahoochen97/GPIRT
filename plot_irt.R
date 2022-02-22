args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  SEED = 1234
  C = 5
  n = 1000
  m = 10
  TYPE = "GP"
}
if (length(args)==5){
  SEED = as.integer(args[1])
  C = as.integer(args[2])
  n = as.integer(args[3])
  m = as.integer(args[4])
  TYPE = args[5]
}

library(ggplot2)
library(KRLS)
library(dplyr)
HYP = paste(TYPE, "_C_", C, '_n_', n, '_m_', m, '_SEED_', SEED, sep="")
load(file=paste("./data/", HYP, ".RData" , sep=""))

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

xs = seq(-5,5,0.01)
idx = 1:1001
if(TYPE=="2PL"){
  for(j in 1:m){
    probs = getprobs_gpirt(xs[idx], matrix(alpha[j]+beta[j]*xs[idx], ncol=1),
                           matrix(thresholds[j,],nrow=1))
    p = ggplot(probs, aes(x=xs, y=p, group=order, color=factor(order))) +
      geom_line(size=2) +ggtitle(paste("True 2PL IRT q",j, sep="")) +
      theme(plot.title = element_text(hjust = 0.5))
    print(p)
    
    tmp = probs %>% 
      group_by(xs) %>%
      summarize(icc=sum(order*p))
    q = plot(xs, tmp$icc)
    print(q)
  }
}
if(TYPE=="GP"){
  for(j in 1:m){
    K = gausskernel(anchor_xs[j,], sigma=SIGMA)
    inv_K = ginv(K)
    K1 = matrix(0, nrow=length(xs), ncol=ncol(anchor_xs))
    for ( i in 1:length(xs) ) {
      K1[i,] = dnorm(xs[i]-anchor_xs[j,], sd=SIGMA)/dnorm(0, sd=SIGMA)
    }
    irfs = K1 %*% inv_K %*% anchor_ys[j,]
    probs = getprobs_gpirt(xs[idx], irfs, matrix(thresholds[j,],nrow=1))
    q = ggplot(probs, aes(x=xs, y=p, group=order, color=factor(order))) +
      geom_line(size=2) +ggtitle(paste("True GP IRT q",j, sep="")) +
      theme(plot.title = element_text(hjust = 0.5))
    
    tmp = probs %>% 
          group_by(xs) %>%
          summarize(icc=sum(order*p))
    pdf(file=paste("./figures/trueiccq",j,".pdf", sep=""))
    plot(xs, tmp$icc)
    dev.off()
    # ggsave(paste("./figures/trueirfq",j,".pdf", sep=""), plot=q, width = 7, height = 4, units = "in")
  }
}

# plot(pred_theta, data[,9], pch=4, ylim=c(-2,5))
# points(pred_theta,rowMeans(samples$f[,9,]))
# for (c in 2:C) {
#   abline(h = mean(samples$threshold[,c]), col="red", lwd=1, lty=1)
# }

# hist(colMeans(samples$theta))
# cor(theta,colMeans(samples$theta))
# for(i in 1:16){plot(xs[idx], log(samples$IRFs[idx,i]/(1-samples$IRFs[idx,i])))}

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

# if (PLOT){
#   idx = (as.integer(min(pred_theta)*100+500)):(as.integer(max(pred_theta)*100+500))
#   for (j in 1:m) {
#     probs = getprobs_catSurv(xs[idx], sdo_cat@discrimination[[paste("q",j, sep="")]],
#                               sdo_cat@difficulty[[paste("q",j, sep="")]])
#   }
# }
