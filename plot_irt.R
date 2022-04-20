args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

if (length(args)==0) {
  SEED = 1
  C = 5
  n = 50
  m = 20
  horizon = 1
  TYPE = "GP"
}
if (length(args)==6){
  SEED = as.integer(args[1])
  C = as.integer(args[2])
  n = as.integer(args[3])
  m = as.integer(args[4])
  horizon = as.integer(args[5])
  TYPE = args[6]
}
library(ggplot2)
library(KRLS)
library(dplyr)
source("getprob_gpirt.R")
HYP = paste(TYPE, "_C_", C, '_n_', n, '_m_', m, '_h_', horizon, '_SEED_', SEED, sep="")
load(file=paste("./data/", HYP, ".RData" , sep=""))

xs = seq(-5,5,0.01)
idx = (as.integer(min(theta)*100+500)):(as.integer(max(theta)*100+500))

for(l in 1:5){
  # plot latent f
  source("true_irf.R")
  pdf(file=paste("./figures/truefq",l,".pdf", sep=""))
  plot(xs[idx], irfs)
  dev.off()
  
  # plot item response function
  probs = getprobs_gpirt(xs[idx], irfs, matrix(thresholds,nrow=1))
  q = ggplot(probs, aes(x=xs, y=p, group=order, color=factor(order))) +
    geom_line(size=2) +ggtitle(paste("True GP IRT q",l, sep="")) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste("./figures/trueirfq",l,".pdf", sep=""), plot=q, width = 7, height = 4, units = "in")
  
  # plot item characteristics curve
  tmp = probs %>% 
    group_by(xs) %>%
    summarize(icc=sum(order*p))
  
  pdf(file=paste("./figures/trueiccq",l,".pdf", sep=""))
  plot(xs[idx], tmp$icc, col=2, lty=1)
  # source('~/Documents/GitHub/OrdGPIRT/gpirt-sdo.R')
  # lines(xs[idx], gpirt_iccs[,l],col=3)
  # source('~/Documents/GitHub/OrdGPIRT/2PL.R')
  # lines(xs[idx], grm_iccs[,l],col=4)
  # source('~/Documents/GitHub/OrdGPIRT/bgrm_logit.R')
  # lines(xs[idx], bgrm_iccs[,l],col=1)
  # legend(x = "topright", 
  #        legend=c("truth", "gpirt", "grm", "bgrm"),
  #        lty = c(1, 1, 1, 1),           # Line types
  #        col = c(2, 3, 4, 1),           # Line colors
  #        lwd = 2)
  dev.off()
}
  
  # for(j in 1:5){
  #   probs = getprobs_2PL(xs[idx], betas[[paste('Item ', j, sep='')]])
  # 
  #   tmp = probs %>%
  #     group_by(xs) %>%
  #     summarize(icc=sum(order*p))
  #   plot(xs[idx], tmp$icc)
  # }

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


