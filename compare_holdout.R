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

results = data.frame(matrix(ncol = 5, nrow = 0))
colnames(results) = c("type","measure","SEED","pvalue","diff_mean")
# train 0/test h, acc, seed 1, gpirt, 0.5, E[gpirt]-E[doirt]

# ttest_acc = matrix(0, nrow = MAXSEED, ncol = TEST_YEAR-TRAIN_END_YEAR+1)
# ttest_lls = matrix(0, nrow = MAXSEED, ncol = TEST_YEAR-TRAIN_END_YEAR+1)
# ttest_response = matrix(0, nrow = MAXSEED, ncol = TEST_YEAR-TRAIN_END_YEAR+1)
# ttest_prediction = matrix(0, nrow = MAXSEED, ncol = TEST_YEAR-TRAIN_END_YEAR+1)

for (SEED in 1:MAXSEED){
  HYP = paste(MODELS[1], "_", DATANAME, "_holdout_SEED_", SEED , sep="")
  load(file=paste("./results/", HYP, ".RData", sep=""))
  train_acc1 = train_acc
  train_lls1 = train_lls
  train_response1 = train_response
  train_prediction1 = train_prediction 
  HYP = paste(MODELS[2], "_", DATANAME, "_holdout_SEED_", SEED , sep="")
  load(file=paste("./results/", HYP, ".RData", sep=""))
  train_acc2 = train_acc
  train_lls2 = train_lls
  train_response2 = train_response
  train_prediction2 = train_prediction 
  
  tmp = t.test(train_acc1, train_acc2)
  results[nrow(results)+1,] = c(0, "acc", SEED, tmp[["p.value"]],
                     tmp[['estimate']][[1]]-tmp[['estimate']][[2]])
  tmp = t.test(train_lls1, train_lls2)
  results[nrow(results)+1,] = c(0, "ll", SEED, tmp[["p.value"]],
                                tmp[['estimate']][[1]]-tmp[['estimate']][[2]])
  tmp = t.test(train_lls1, train_lls2)
  results[nrow(results)+1,] = c(0, "ll", SEED, tmp[["p.value"]],
                                tmp[['estimate']][[1]]-tmp[['estimate']][[2]])
        
  for(h in 1:horizon) {
    h_ = h+1999-TRAIN_END_YEAR
    if(h_>0){
      HYP = paste(MODELS[1], "_", DATANAME, "_holdout_SEED_", SEED , sep="")
      load(file=paste("./results/", HYP, ".RData", sep=""))
      test_acc1 = test_acc[[h_]]
      test_lls1 = test_lls[[h_]]
      test_response1 = test_response[[h_]]
      test_prediction1 = test_prediction[[h_]]
      
      HYP = paste(MODELS[2], "_", DATANAME, "_holdout_SEED_", SEED , sep="")
      load(file=paste("./results/", HYP, ".RData", sep=""))
      test_acc2 = test_acc[[h_]]
      test_lls2 = test_lls[[h_]]
      test_response2 = test_response[[h_]]
      test_prediction2 = test_prediction[[h_]]
    }
  }
}

