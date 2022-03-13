args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

if (length(args)==0) {
  SEED = 1
  C = 5
  n = 50
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

source("getprob_gpirt.R")
HYP = paste(TYPE, "_C_", C, '_n_', n, '_m_', m, '_SEED_', SEED, sep="")
load(file=paste("./data/", HYP, ".RData" , sep=""))

MODELS = c("gpirt","grm")#, "bgrm")

results = matrix(0, nrow = 7, ncol = length(MODELS))
for(i in 1:length(MODELS)){
  load(file=paste("./results/", MODELS[i], "_", HYP, ".RData" , sep=""))
  results[1,i] = cor(theta, pred_theta)
  results[2,i] = mean(train_lls)
  results[3,i] = mean(train_acc)
  results[4,i] = mean(pred_lls)
  results[5,i] = mean(pred_acc)
  results[6,i] = mean(cor_icc)
  results[7,i] = mean(rmse_icc)
}

write.csv(results, file=paste("./results/", HYP, ".csv" , sep=""))
