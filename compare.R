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
    # cor_theta[SEED, i] = result[1,i+1]
    # train_ll[SEED, i] = result[2,i+1]
    # train_acc[SEED, i] = result[3,i+1]
    # pred_ll[SEED, i] = result[4,i+1]
    # pred_acc[SEED, i] = result[5,i+1]
    # cor_icc[SEED, i] = result[6,i+1]
    # rmse_icc[SEED, i] = result[7,i+1]
    cor_theta[SEED, i] = result[1,2]
    train_ll[SEED, i] = result[2,2]
    train_acc[SEED, i] = result[3,2]
    pred_ll[SEED, i] = result[4,2]
    pred_acc[SEED, i] = result[5,2]
    cor_icc[SEED, i] = result[6,2]
    rmse_icc[SEED, i] = result[7,2]
    }
}

COR_THETA_ERR = apply(abs(cor_theta), 2, sd)/sqrt(MAXSEED)
TRAIN_LL_ERR =  apply(train_ll, 2, sd)/sqrt(MAXSEED)
TRAIN_ACC_ERR =  apply(train_acc, 2, sd)/sqrt(MAXSEED)
PRED_LL_ERR =  apply(pred_ll, 2, sd)/sqrt(MAXSEED)
PRED_ACC_ERR = apply(pred_acc, 2, sd)/sqrt(MAXSEED)
COR_ICC_ERR =  apply(cor_icc, 2, sd)/sqrt(MAXSEED)
RMSE_ICC_ERR =  apply(rmse_icc, 2, sd)/sqrt(MAXSEED)

print(colMeans(abs(cor_theta)))
print(COR_THETA_ERR)
print(colMeans((train_ll)))
print(TRAIN_LL_ERR)
print(colMeans(abs(train_acc)))
print(TRAIN_ACC_ERR)
print(colMeans((pred_ll)))
print(PRED_LL_ERR)
print(colMeans(abs(pred_acc)))
print(PRED_ACC_ERR)
print(colMeans(abs(cor_icc)))
print(COR_ICC_ERR)
print(colMeans(rmse_icc))
print(RMSE_ICC_ERR)

HYP = paste("_C_", C, '_n_', n, '_m_', m, '_h_', horizon, '_SEED_', SEED, sep="")
# write.csv(results, file=paste("./results/compare", HYP, ".csv" , sep=""))
