R_path="~/R/x86_64-redhat-linux-gnu-library/4.0"
.libPaths(R_path)

library(dplyr)
library(ltm)

gpirt_path = "~/Documents/Github/OrdGPIRT/VAA"
setwd(gpirt_path)

# initialize with grm
theta_init = matrix(0, nrow = n, ncol = horizon)

data_reshape = matrix(0, nrow=n, ncol=0)
for(h in 1:horizon){
  mask = unlist(lapply(item_names, FUN=function(x) 
    substr(x,2,3)==substr(toString(years[h]),3,4)))
  for(j in 1:num_item_years[h]){
    tmp = data[,j,h]
    if(sum(!is.na(unique(tmp)))==C){
      data_reshape = cbind(data_reshape, tmp)
    }
  }
}

results = grm(data = data_reshape, na.action = NULL)
pred_theta = factor.scores(results, resp.patterns = data_reshape)
pred_theta = pred_theta$score.dat[,"z1"]

for(h in 1:horizon){
  theta_init[,h] = pred_theta
}