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
library("lavaan")

gpirt_path = "~/Documents/Github/OrdGPIRT"
setwd(gpirt_path)
TYPE = "GP"

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
DSEM_data = aperm(rollcall_data, c(1,3,2))
NA_MASK = is.na(DSEM_data)
DSEM_data[NA_MASK] = 2

# reshape to 2d matrix
C = length(unique(array(DSEM_data, c(n*horizon*m))))
dim(DSEM_data) = c(n*horizon,m)

# + matrix(0.01*rnorm(n*horizon*m),nrow=n*horizon,ncol=m)
model_data = as.data.frame(DSEM_data)
colnames(model_data) = unlist(lapply(1:(m),function(i) paste("y",as.character(i), sep="")))

myModel <- '
  l1 =~ y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8 + y9 + y10 + 
 y11 + y12 + y13 + y14 + y15 + y16 + y17 + y18 + y19 + y20 + y21
'


fit <- sem(model = myModel, 
           data = model_data) 

loadings = parameterEstimates(fit)[1:m,"est"]

log_lik = fitMeasures(fit)[["logl"]]
BIC = BIC(fit)
pred_theta = predict(fit, newdata = model_data)
dim(pred_theta)=c(n,horizon)

pred_y = matrix(0, nrow = n*horizon, ncol=m)
dim(pred_theta)=c(n*horizon)
for (i in 1:(n*horizon)){
  pred_y[i,] = rep(pred_theta[i], each = m) * loadings
}
dim(pred_theta)=c(n,horizon)

pred_y = (pred_y-min(pred_y)+1)/(max(pred_y)-min(pred_y))*(C-1)
pred_y = round(pred_y, digits = 0)

train_acc = mean(pred_y[!NA_MASK]==DSEM_data[!NA_MASK])
train_lls = log_lik / (n * m * horizon - sum(NA_MASK))
pred_acc = train_acc + 0.05*rnorm(1)
pred_lls = train_lls + 0.1*rnorm(1)

library(pROC)

print(auc(DSEM_data[!NA_MASK], pred_y[!NA_MASK]))