args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

if (length(args)==0) {
  MAXSEED = 1
  C = 5
  n = 50
  m = 10
  TYPE = "GP"
}
if (length(args)==5){
  MAXSEED = as.integer(args[1])
  C = as.integer(args[2])
  n = as.integer(args[3])
  m = as.integer(args[4])
  TYPE = args[5]
}

source("getprob_gpirt.R")
load(file=paste("./data/", HYP, ".RData" , sep=""))

MODELS = c("gpirt","grm", "bgrm")

cor_theta = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
train_ll = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
pred_ll = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
train_acc = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
pred_acc = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
cor_icc = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
rmse_icc = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
for(i in 1:length(MODELS)){
  for (SEED in 1:MAXSEED) {
    HYP = paste(TYPE, "_C_", C, '_n_', n, '_m_', m, '_SEED_', SEED, sep="")
    result = read.csv(file=paste("./results/", MODELS[i], "_", HYP, ".csv" , sep=""))
    cor_theta[SEED, i] = result[1,i]
    train_ll[SEED, i] = result[2,i]
    train_acc[SEED, i] = result[3,i]
    pred_ll[SEED, i] = result[4,i]
    pred_acc[SEED, i] = result[5,i]
    cor_icc[SEED, i] = result[6,i]
    rmse_icc[SEED, i] = result[7,i]
    }
}
# write.csv(results, file=paste("./results/compare", HYP, ".csv" , sep=""))
