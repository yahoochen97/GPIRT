library(dplyr)
library(ggplot2)
library(stats)
library(ltm)

gpirt_path = "~/Documents/Github/OrdGPIRT"
setwd(gpirt_path)

# esm.BFI21 and esm.BFI36
data = read.csv("./data/esm_w1_redacted.csv")
data = data[,c("esm.IDnum.w1","esm.PRO01.w1","esm.PRO03.w1","esm.PRO04.w1",
               "esm.BFI09.w1","esm.BFI04.w1","esm.BFI19.w1")]
colnames(data) = c("SID", "freq", "hourblock","day","N_relax", "N_depress", "N_worried")
data = data[!is.na(data$N_relax),]
data = data[!is.na(data$N_depress),]
data = data[!is.na(data$N_worried),]

# data = data[(data$day<50) & (data$SID<5000), ]

num_days = length(unique(data$day))
num_hours = length(unique(data$hourblock))
unique_ids = unique(data$SID)
n = length(unique_ids)
m = 3 # three questions across all time periods
horizon = num_days*num_hours
C = 5

ems_data = array(array(NA, n*m*horizon), c(n,m, horizon))
for(iter in 1:nrow(data)){
  id = data$SID[iter]
  hourblock = data$hourblock[iter]
  day = data$day[iter]
  N_relax = data$N_relax[iter]
  N_depress = data$N_depress[iter]
  N_worried = data$N_worried[iter]
  i = which(id==unique_ids)
  h = (day-1)*num_hours + hourblock
  ems_data[i,1,h] = N_relax
  ems_data[i,2,h] = N_depress
  ems_data[i,2,h] = N_worried
}

data = ems_data

results = grm(data = data, na.action = NULL)
betas = results$coefficients
pred_theta = factor.scores(results, resp.patterns = data_train)
pred_theta = pred_theta$score.dat[,m+3]

xs = seq(-5,5,0.01)
idx = 301:701
gpirt_iccs = array(array(0, length(xs[idx])*m*horizon),
                   c(length(xs[idx]),m, horizon))

for(h in 1:horizon){
  for (j in 1:m) {
    tmp = c()
    for (i in 1:length(xs[idx])) {
      ps = rep(0, C+1)
      ps[C+1] = 1
      beta_coef = betas[[paste('Item ', j, sep='')]]
      for (c in 1:(C-1)){
        lp = beta_coef[[c]] - xs[idx][i]*beta_coef[[C]]
        ps[c+1] = 1 / (1+ exp(-lp))
      }
      
      p = rep(0, C)
      for (c in 1:C) {
        p[c] = ps[c+1] - ps[c]
      }
      tmp = c(tmp, sum(p*(1:C)))
    }
    grm_iccs[,j] = tmp
  }
}
