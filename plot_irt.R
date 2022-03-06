args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

if (length(args)==0) {
  SEED = 1
  C = 5
  n = 50
  m = 20
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
source("getprob_gpirt.R")
HYP = paste(TYPE, "_C_", C, '_n_', n, '_m_', m, '_SEED_', SEED, sep="")
load(file=paste("./data/", HYP, ".RData" , sep=""))

xs = seq(-5,5,0.01)
idx = (as.integer(min(theta)*100+500)):(as.integer(max(theta)*100+500))
if(TYPE=="GP"){
  for(j in 1:5){
    source("true_irf.R")
    pdf(file=paste("./figures/truefq",j,".pdf", sep=""))
    plot(xs[idx], irfs)
    points(anchor_xs[j,],anchor_ys[j,])
    dev.off()
    
    probs = getprobs_gpirt(xs[idx], irfs, matrix(thresholds,nrow=1))
    q = ggplot(probs, aes(x=xs, y=p, group=order, color=factor(order))) +
      geom_line(size=2) +ggtitle(paste("True GP IRT q",j, sep="")) +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(paste("./figures/trueirfq",j,".pdf", sep=""), plot=q, width = 7, height = 4, units = "in")
    
    tmp = probs %>% 
          group_by(xs) %>%
          summarize(icc=sum(order*p))
    pdf(file=paste("./figures/trueiccq",j,".pdf", sep=""))
    plot(xs[idx], tmp$icc)
    dev.off()
  }
  
  for(j in 1:5){
    probs = getprobs_2PL(xs[idx], betas[[paste('Item ', j, sep='')]])
    # p = ggplot(probs, aes(x=xs, y=p, group=order, color=factor(order))) +
    #   geom_line(size=2) +ggtitle(paste("2PL IRT q",j, sep="")) +
    #   theme(plot.title = element_text(hjust = 0.5))
    # print(p)

    tmp = probs %>%
      group_by(xs) %>%
      summarize(icc=sum(order*p))
    plot(xs[idx], tmp$icc)
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
