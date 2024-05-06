args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

if (length(args)==0) {
  MAXSEED = 1
  C = 5
  n = 100
  m = 10
  horizon = 10
  TYPE = "DSEM"
  CONSTANT_IRF = 0
  DATA_TYPE = "GP"
}
if (length(args)==8){
  MAXSEED = as.integer(args[1])
  C = as.integer(args[2])
  n = as.integer(args[3])
  m = as.integer(args[4])
  horizon = as.integer(args[5])
  TYPE = args[6]
  CONSTANT_IRF = as.integer(args[7])
  DATA_TYPE = args[8]
}

MODELS = c("gpirt","grm", "bgrm")
MODELS = c("GP", "CST", "RDM", "BRW")
MODELS = c("ggum", "graded", "gpcm", "sequential")
MODELS = c("DSEM")

cor_theta = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
sd_theta = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
ll_theta = matrix(0, nrow = MAXSEED, ncol = length(MODELS))

train_llss = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
pred_llss = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
train_accss = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
pred_accss = matrix(0, nrow = MAXSEED, ncol = length(MODELS))

cor_iccss = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
rmse_iccss = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
theta_rhatsss = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
for(i in 1:length(MODELS)){
  for (SEED in 1:MAXSEED) {
    tryCatch(
      {
        HYP = paste(MODELS[i], "_C_", C, '_n_', n, '_m_', m, '_h_', horizon,'_CSTIRF_', CONSTANT_IRF , '_SEED_', SEED, sep="")
        result = read.csv(file=paste("./results/", HYP, ".csv" , sep=""))
      }, error = function(cond){
        HYP = paste(MODELS[i], "_C_", C, '_n_', n, '_m_', m, '_h_', horizon,'_CSTIRF_', CONSTANT_IRF , '_SEED_', 1, sep="")
        result = read.csv(file=paste("./results/", HYP, ".csv" , sep=""))
      }
    )
    cor_theta[SEED, i] = result[1,2]
    sd_theta[SEED, i] = result[2,2]
    ll_theta[SEED, i] = result[3,2]
    train_llss[SEED, i] = result[4,2]
    train_accss[SEED, i] = result[5,2]
    pred_llss[SEED, i] = result[6,2]
    pred_accss[SEED, i] = result[7,2]
    cor_iccss[SEED, i] = result[8,2]
    rmse_iccss[SEED, i] = result[9,2]
    theta_rhatsss[SEED, i] = result[10,2]
    # load(paste("./results/gpirt_", HYP, ".RData" , sep=""))
    # train_llss[SEED, i] = mean(train_lls[!is.na(train_lls)])
    # train_accss[SEED, i] = mean(train_acc[!is.na(train_acc)])
    # pred_llss[SEED, i] = mean(pred_lls[!is.na(pred_lls)])
    # pred_accss[SEED, i] = mean(pred_acc[!is.na(pred_acc)])
    }
}

train_llss[is.infinite(train_llss)] = mean(train_llss[!is.infinite(train_llss)])
pred_llss[is.infinite(pred_llss)] = mean(pred_llss[!is.infinite(pred_llss)])

COR_THETA_ERR = apply(abs(cor_theta), 2, sd)/sqrt(MAXSEED)
TRAIN_LL_ERR =  apply(train_llss, 2, sd)/sqrt(MAXSEED)
TRAIN_ACC_ERR =  apply(train_accss, 2, sd)/sqrt(MAXSEED)
PRED_LL_ERR =  apply(pred_llss, 2, sd)/sqrt(MAXSEED)
PRED_ACC_ERR = apply(pred_accss, 2, sd)/sqrt(MAXSEED)
COR_ICC_ERR =  apply(cor_iccss, 2, sd)/sqrt(MAXSEED)
RMSE_ICC_ERR =  apply(rmse_iccss, 2, sd)/sqrt(MAXSEED)

print(HYP)
print("cor theta")
print(colMeans(abs(cor_theta)))
print(COR_THETA_ERR)
print("mean theta rhat")
print(colMeans(theta_rhatsss))
print("ll theta")
print(colMeans(ll_theta))
print(apply(ll_theta, 2, sd)/sqrt(MAXSEED))
print("cor icc")
print(colMeans(abs(cor_iccss)))
print(COR_ICC_ERR)
print("rmse icc")
print(colMeans(rmse_iccss))
print(RMSE_ICC_ERR)
print("\n")

print("train ll")
print(colMeans((train_llss)))
print(TRAIN_LL_ERR)
print("train acc")
print(colMeans(abs(train_accss)))
print(TRAIN_ACC_ERR)
print("pred ll")
print(colMeans((pred_llss)))
print(PRED_LL_ERR)
print("pred acc")
print(colMeans(abs(pred_accss)))
print(PRED_ACC_ERR)
print("\n")

HYP = paste("_C_", C, '_n_', n, '_m_', m, '_h_', horizon,'_CSTIRF_', CONSTANT_IRF , '_SEED_', SEED, sep="")
# write.csv(results, file=paste("./results/compare", HYP, ".csv" , sep=""))
