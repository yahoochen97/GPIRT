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
horizon = TEST_YEAR - TRAIN_START_YEAR + 1

ttest_acc = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
ttest_lls = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
ttest_response = matrix(0, nrow = MAXSEED, ncol = length(MODELS))
ttest_prediction = matrix(0, nrow = MAXSEED, ncol = length(MODELS))

for(h in 1:horizon) {
  h_ = h+1999-TRAIN_END_YEAR
  if(h_>0){
    for (SEED in 1:MAXSEED){
      for(i in 1:length(MODELS)){
        HYP = paste(MODELS[1], "_", DATANAME, "_holdout_SEED_", SEED , sep="")
        load(file=paste("./results/", HYP, ".RData", sep=""))
        train_acc
        train_lls
        train_response 
        train_prediction
        test_acc[[h_]]
        test_lls[[h_]]
        test_response[[h_]]
        test_prediction[[h_]]
      }
    }
  } 
}

