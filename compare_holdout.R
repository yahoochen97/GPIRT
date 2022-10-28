args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

if (length(args)==0) {
  TRAIN_START_YEAR = 2000
  TRAIN_END_YEAR = 2015
  TEST_YEAR = 2020
  MAXSEED = 25
  DATANAME = 'SupremeCourt'
}
if (length(args)==5){
  TRAIN_START_YEAR = as.integer(args[1])
  TRAIN_END_YEAR = as.integer(args[2])
  TEST_YEAR = as.integer(args[3])
  MAXSEED = as.integer(args[4])
  DATANAME = args[5]
}

MODELS = c("gpirt","doirt")

cor_theta = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
for(i in 1:length(MODELS)){
  for (SEED in 1:MAXSEED) {
    HYP = paste(MODELS[i], "_", DATANAME, "_holdout_SEED_", SEED , sep="")
    load(file=paste("./results/", HYP, ".RData", sep=""))
    cor_theta[SEED, i] = result[1,2]
  }
}

