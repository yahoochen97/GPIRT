args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

if (length(args)==0) {
  SEED = 1
  C = 5
  n = 100
  m = 10
  horizon = 10
  TYPE = "GPDM"
  CONSTANT_IRF = 0
  DATA_TYPE = "GP"
}
if (length(args)==8){
  SEED = as.integer(args[1])
  C = as.integer(args[2])
  n = as.integer(args[3])
  m = as.integer(args[4])
  horizon = as.integer(args[5])
  TYPE = args[6]
  CONSTANT_IRF = as.integer(args[7])
  DATA_TYPE = args[8]
}

source("getprob_gpirt.R")
HYP = paste(DATA_TYPE, "_C_", C, '_n_', n, '_m_', m, '_h_', horizon,'_CSTIRF_', CONSTANT_IRF , '_SEED_', SEED, sep="")
load(file=paste("./data/", HYP, ".RData" , sep=""))
HYP = paste(TYPE, "_C_", C, '_n_', n, '_m_', m, '_h_', horizon,'_CSTIRF_', CONSTANT_IRF , '_SEED_', SEED, sep="")

MODELS = c("gpirt")

results = matrix(0, nrow = 10, ncol = length(MODELS))
for(i in 1:length(MODELS)){
  load(file=paste("./results/", MODELS[i], "_", HYP, ".RData" , sep=""))
  cor_theta = c()
  for (h in 1:horizon) {
    cor_theta = c(cor_theta, cor(theta[,h], pred_theta[,h]))
  }
  theta_rhats = c(1)
  results[1,i] = mean(abs(cor_theta), na.rm = TRUE)
  results[2,i] = mean(pred_theta_sd[!is.infinite(pred_theta_sd)], na.rm = TRUE)
  results[3,i] = mean(pred_theta_ll[!is.infinite(pred_theta_ll)], na.rm = TRUE)
  results[4,i] = mean(train_lls[!is.infinite(train_lls)], na.rm = TRUE)
  results[5,i] = mean(train_acc, na.rm = TRUE)
  results[6,i] = mean(pred_lls[!is.infinite(pred_lls)], na.rm = TRUE)
  results[7,i] = mean(pred_acc, na.rm = TRUE)
  results[8,i] = mean(abs(cor_icc), na.rm = TRUE)
  results[9,i] = mean(rmse_icc, na.rm = TRUE)
  results[10,i] = mean(theta_rhats, na.rm = TRUE)
}

write.csv(results, file=paste("./results/", HYP, ".csv" , sep=""))
