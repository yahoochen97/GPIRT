args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

if (length(args)==0) {
  SEED = 1
  C = 2
  n = 100
  m = 20
  horizon = 10
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

source("getprob_gpirt.R")
HYP = paste(TYPE, "_C_", C, '_n_', n, '_m_', m, '_h_', horizon, '_SEED_', SEED, sep="")
load(file=paste("./data/", HYP, ".RData" , sep=""))

MODELS = c("gpirt","grm", "bgrm")
MODELS = c("gpirt")

results = matrix(0, nrow = 7, ncol = length(MODELS))
for(i in 1:length(MODELS)){
  load(file=paste("./results/", MODELS[i], "_", HYP, ".RData" , sep=""))
  cor_theta = c()
  for (h in 1:horizon) {
    cor_theta = c(cor_theta, cor(theta[,h], pred_theta[,h]))
  }
  results[1,i] = mean(cor_theta)
  results[2,i] = mean(train_lls)
  results[3,i] = mean(train_acc)
  results[4,i] = mean(pred_lls)
  results[5,i] = mean(pred_acc)
  results[6,i] = mean(cor_icc)
  results[7,i] = mean(rmse_icc)
}

write.csv(results, file=paste("./results/", HYP, ".csv" , sep=""))
