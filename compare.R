args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

if (length(args)==0) {
  MAXSEED = 100
  C = 2
  n = 100
  m = 20
  horizon = 10
  TYPE = "GP"
}
if (length(args)==6){
  MAXSEED = as.integer(args[1])
  C = as.integer(args[2])
  n = as.integer(args[3])
  m = as.integer(args[4])
  horizon = as.integer(args[5])
  TYPE = args[6]
}

MODELS = c("gpirt","grm", "bgrm")
MODELS = c("GP", "CST", "RDM")

cor_theta = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
train_ll = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
pred_ll = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
train_acc = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
pred_acc = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
cor_icc = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
rmse_icc = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
for(i in 1:length(MODELS)){
  for (SEED in 1:MAXSEED) {
    HYP = paste(MODELS[i], "_C_", C, '_n_', n, '_m_', m, '_h_', horizon, '_SEED_', SEED, sep="")
    result = read.csv(file=paste("./results/", HYP, ".csv" , sep=""))
    cor_theta[SEED, i] = result[1,i+1]
    train_ll[SEED, i] = result[2,i+1]
    train_acc[SEED, i] = result[3,i+1]
    pred_ll[SEED, i] = result[4,i+1]
    pred_acc[SEED, i] = result[5,i+1]
    cor_icc[SEED, i] = result[6,i+1]
    rmse_icc[SEED, i] = result[7,i+1]
    }
}
print(colMeans(abs(cor_theta)))
print(colMeans((train_ll)))
print(colMeans(abs(train_acc)))
print(colMeans((pred_ll)))
print(colMeans(abs(pred_acc)))
print(colMeans(abs(cor_icc)))
print(colMeans(rmse_icc))

HYP = paste("_C_", C, '_n_', n, '_m_', m, '_h_', horizon, '_SEED_', SEED, sep="")
write.csv(results, file=paste("./results/compare", HYP, ".csv" , sep=""))
