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

ttest_acc = array(array(0, MAXSEED*(TEST_YEAR-TRAIN_END_YEAR+1)*length(MODELS)),
                  c(MAXSEED,(TEST_YEAR-TRAIN_END_YEAR+1),length(MODELS)))
ttest_lls = array(array(0, MAXSEED*(TEST_YEAR-TRAIN_END_YEAR+1)*length(MODELS)),
                  c(MAXSEED,(TEST_YEAR-TRAIN_END_YEAR+1),length(MODELS)))

for (SEED in 1:MAXSEED){
  for(k in 1:length(MODELS)){
    HYP = paste(MODELS[k], "_", DATANAME, "_holdout_SEED_", SEED , sep="")
    load(file=paste("./results/", HYP, ".RData", sep=""))
    ttest_acc[SEED,1,k] = mean(train_acc)
    ttest_lls[SEED,1,k] = mean(train_lls)
    
    for(h in 1:horizon) {
      h_ = h+1999-TRAIN_END_YEAR
      if(h_>0){
        ttest_acc[SEED,1+h_,k] = mean(test_acc[[h_]])
        ttest_lls[SEED,1+h_,k] = mean(test_lls[[h_]])

      }
    }
    # ttest_acc[SEED,2,k] = mean(unlist(test_acc))
    # ttest_lls[SEED,2,k] = mean(unlist(test_lls))
  }
}

# t tests
results = data.frame(matrix(ncol = 4, nrow = 0))
colnames(results) = c("type","measure","pvalue","diff_mean")
# train 0/test h, acc, pvalue0.5, E[gpirt]-E[doirt]


for(h in 1:horizon) {
  h_ = h+1999-TRAIN_END_YEAR
  if(h_>=0){
    tmp = t.test(ttest_acc[,h_+1,1], ttest_acc[,h_+1,2], paired = TRUE)
    results[nrow(results)+1,] = c(h_, "acc", tmp[["p.value"]],
                       tmp[['estimate']][[1]])
    # -tmp[['estimate']][[2]]
    tmp = t.test(ttest_lls[,h_+1,1], ttest_lls[,h_+1,2], paired = TRUE)
    results[nrow(results)+1,] = c(h_, "ll", tmp[["p.value"]],
                                  tmp[['estimate']][[1]])
    # -tmp[['estimate']][[2]]
  }
}

write.csv(results, paste("./results/",DATANAME,"_holdout.csv", sep=""))