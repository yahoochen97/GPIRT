# gpirt_path = "~/Documents/Github/OrdGPIRT"
# setwd(gpirt_path)
# library(Rcpp)
# Rcpp::compileAttributes()
# install.packages(gpirt_path, type="source", repos = NULL ,lib=R_path, INSTALL_opts = '--no-lock')
# setwd("./TAPS")

data = read.csv("TAPS-econ-data-long.csv")
# <= 20
data = data[data$wave>=20,]

data = data[data$wave>=60,]

all_ids = unique(data$WUSTLID)
n = length(all_ids)
questions = c("ECON2A", "ECON2B", "ECON1A", "ECON1B", "ECON4A", "ECON4B")
m = length(questions)
all_waves = unique(data$wave)
horizon = length(all_waves)

freqs = rep(0, n)
for(i in 1:n){
  freqs[i] = sum(!is.na(data[data$WUSTLID==all_ids[i], questions]))
}

# 45 # 165
# data = data[data$WUSTLID %in% all_ids[freqs>=165],]

if(exists("TEST_YEAR")){
  # data = data[data$WUSTLID %in% all_ids[freqs>=168],]
  
  data = data[data$WUSTLID %in% all_ids[freqs==28],]
}

all_ids = unique(data$WUSTLID)
n = length(all_ids)

print(paste("there are ", n, " units.", sep=""))

gpirt_data = array(array(NA, n*m*horizon), c(n,m,horizon))
for(h in 1:horizon){
  print(h)
  for(i in 1:n){
    for(j in 1:length(questions)){
      gpirt_data[i,j,h] = data[data$wave==all_waves[h] & data$WUSTLID==all_ids[i], questions[j] ]
    }
  }
}

data_reshape = matrix(0, nrow=n*horizon, ncol=m)
for(h in 1:horizon){
  data_reshape[(1+(h-1)*n):(h*n),] = gpirt_data[,,h]
}

results_2PL = grm(data = data_reshape, na.action = NULL)
betas_2PL = results_2PL$coefficients
theta_2PL = factor.scores(results_2PL, resp.patterns = data_reshape)
theta_2PL = theta_2PL$score.dat[,m+3]

theta_init = matrix(0.0, nrow = n, ncol = horizon)
for(h in 1:horizon){
  theta_init[,h] = theta_2PL[(1+(h-1)*n):(h*n)]
}

C = 5