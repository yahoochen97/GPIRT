args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

if (length(args)==0) {
  SEED = 1
  C = 2
  n = 50
  m = 10
  horizon = 20
  TYPE = "GP"
  CONSTANT_IRF = 1
}
if (length(args)==7){
  SEED = as.integer(args[1])
  C = as.integer(args[2])
  n = as.integer(args[3])
  m = as.integer(args[4])
  horizon = as.integer(args[5])
  TYPE = args[6]
  CONSTANT_IRF = as.integer(args[7])
}

source("getprob_gpirt.R")
HYP = paste("GP_C_", C, '_n_', n, '_m_', m, '_h_', horizon,'_CSTIRF_', CONSTANT_IRF , '_SEED_', SEED, sep="")
load(file=paste("./data/", HYP, ".RData" , sep=""))
HYP = paste(TYPE, "_C_", C, '_n_', n, '_m_', m, '_h_', horizon,'_CSTIRF_', CONSTANT_IRF , '_SEED_', SEED, sep="")

MODELS = c("gpirt","grm", "bgrm")
MODELS = c("gpirt")

# icc -1 to 1
# idx = 101:301
# cor_icc = matrix(0, nrow=m, ncol=horizon)
# rmse_icc = matrix(0, nrow=m, ncol=horizon)
# 
# for (h in 1:horizon) {
#   for (j in 1:m) {
#     cor_icc[j,h] = cor(gpirt_iccs[idx,j,h], true_iccs[idx,j,h])
#     rmse_icc[j,h] = sqrt(mean((gpirt_iccs[idx,j,h]-true_iccs[idx,j,h])^2))
#   }
# }

results = matrix(0, nrow = 9, ncol = length(MODELS))
for(i in 1:length(MODELS)){
  load(file=paste("./results/", MODELS[i], "_", HYP, ".RData" , sep=""))
  cor_theta = c()
  for (h in 1:horizon) {
    cor_theta = c(cor_theta, cor(theta[,h], pred_theta[,h]))
  }
  results[1,i] = mean(abs(cor_theta), na.rm = TRUE)
  results[2,i] = mean(pred_theta_sd[!is.infinite(pred_theta_sd)], na.rm = TRUE)
  results[3,i] = mean(pred_theta_ll[!is.infinite(pred_theta_ll)], na.rm = TRUE)
  results[4,i] = mean(train_lls, na.rm = TRUE)
  results[5,i] = mean(train_acc, na.rm = TRUE)
  results[6,i] = mean(pred_lls, na.rm = TRUE)
  results[7,i] = mean(pred_acc, na.rm = TRUE)
  results[8,i] = mean(abs(cor_icc), na.rm = TRUE)
  results[9,i] = mean(rmse_icc, na.rm = TRUE)
}

write.csv(results, file=paste("./results/", HYP, ".csv" , sep=""))
