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

MODELS = c("gpirt","SEirt","doirt")
horizon = TEST_YEAR - TRAIN_START_YEAR + 1

ttest_acc = array(array(0, MAXSEED*(TEST_YEAR-TRAIN_END_YEAR+1)*length(MODELS)),
                  c(MAXSEED,(TEST_YEAR-TRAIN_END_YEAR+1),length(MODELS)))
ttest_lls = array(array(0, MAXSEED*(TEST_YEAR-TRAIN_END_YEAR+1)*length(MODELS)),
                  c(MAXSEED,(TEST_YEAR-TRAIN_END_YEAR+1),length(MODELS)))
ttest_sds = array(array(0, MAXSEED*(TEST_YEAR-TRAIN_END_YEAR+1)*length(MODELS)),
                  c(MAXSEED,(TEST_YEAR-TRAIN_END_YEAR+1),length(MODELS)))

for (SEED in 1:MAXSEED){
  for(k in 1:length(MODELS)){
    HYP = paste(MODELS[k], "_", DATANAME, "_holdout_SEED_", SEED , sep="")
    load(file=paste("./results/", HYP, ".RData", sep=""))
    ttest_acc[SEED,1,k] = mean(train_acc)
    ttest_lls[SEED,1,k] = mean(train_lls)
    
    TEST_MASK = (!is.na(gpirt_data)) & (is.na(gpirt_data_train))
    TRAIN_MASK = (!is.na(gpirt_data)) & (!is.na(gpirt_data_train))
    pred_theta_sd_ = array(array(0, nrow(pred_theta_sd)*ncol(gpirt_data)*ncol(pred_theta_sd)),
                          c(nrow(pred_theta_sd),ncol(gpirt_data),ncol(pred_theta_sd)))
    for(j in 1:ncol(gpirt_data)){
      pred_theta_sd_[,j,]=pred_theta_sd
    }
    ttest_sds[SEED,1,k] = mean(unique(pred_theta_sd_[TRAIN_MASK]))
    
    for(h in 1:horizon) {
      h_ = h+1999-TRAIN_END_YEAR
      if(h_>0){
        ttest_acc[SEED,1+h_,k] = mean(test_acc[[h_]])
        ttest_lls[SEED,1+h_,k] = mean(test_lls[[h_]])
        ttest_sds[SEED,1+h_,k] = mean(unique(pred_theta_sd_[,,h][TEST_MASK[,,h]]))
      }
    }
    # ttest_acc[SEED,2,k] = mean(unlist(test_acc))
    # ttest_lls[SEED,2,k] = mean(unlist(test_lls))
  }
}

# t tests
results = data.frame(matrix(ncol = 5, nrow = 0))
colnames(results) = c("type","measure","pvalue","diff_mean","compare")
# train 0/test h, acc, pvalue0.5, E[gpirt]-E[doirt], gpirt vs 


for(h in 1:horizon) {
  h_ = h+1999-TRAIN_END_YEAR
  for(k in 2:length(MODELS)){
    MODEL_COMPARE = paste(MODELS[1], " vs ", MODELS[k], sep="")
    if(h_>=0){
      tmp = t.test(ttest_acc[,h_+1,1], ttest_acc[,h_+1,2])
      results[nrow(results)+1,] = c(h_, "acc", tmp[["p.value"]],
                mean(ttest_acc[,h_+1,1])-mean(ttest_acc[,h_+1,2]), MODEL_COMPARE)
      # -tmp[['estimate']][[2]]
      tmp = t.test(ttest_lls[,h_+1,1], ttest_lls[,h_+1,2])
      results[nrow(results)+1,] = c(h_, "ll", tmp[["p.value"]],
                mean(ttest_lls[,h_+1,1])-mean(ttest_lls[,h_+1,2]),MODEL_COMPARE)
      
      tmp = t.test(ttest_sds[,h_+1,1], ttest_sds[,h_+1,2])
      results[nrow(results)+1,] = c(h_, "sd", tmp[["p.value"]],
                 mean(ttest_sds[,h_+1,1])-mean(ttest_sds[,h_+1,2]),MODEL_COMPARE)
    }
  }
  
}

write.csv(results, paste("./results/",DATANAME,"_holdout.csv", sep=""))