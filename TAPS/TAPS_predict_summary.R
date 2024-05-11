R_path="~/R/x86_64-redhat-linux-gnu-library/4.0"
.libPaths(R_path)
library(rstan)
rstan_options(auto_write = TRUE)
args = commandArgs(trailingOnly=TRUE)
options(show.error.locations = TRUE)

if (length(args)==0) {
  TRAIN_START_YEAR = 1
  TRAIN_END_YEAR = 31
  TEST_YEAR = 41
}

MAXSEED = 25
DRS=c(5, 10, 15, 20, 25)
MODELS = c("DOIRT", "GPIRT")

DRS=c(20)
MODELS = c("NIRT", "LIRT")
MODELS = c("ggum", "graded", "gpcm", "DSEM", "sequential")
MODELS = c("GPDM")

TRAIN_ACC = array(rep(0,MAXSEED*length(DRS)*length(MODELS)), 
                  c(length(MODELS),length(DRS),MAXSEED)) 
TRAIN_LL = array(rep(0,MAXSEED*length(DRS)*length(MODELS)), 
                  c(length(MODELS),length(DRS),MAXSEED)) 
TEST_ACC = array(rep(0,MAXSEED*length(DRS)*10*length(MODELS)), 
                  c(length(MODELS),length(DRS),10,MAXSEED)) 
TEST_LL = array(rep(0,MAXSEED*length(DRS)*10*length(MODELS)), 
                 c(length(MODELS),length(DRS),10,MAXSEED))

RESULTS = data.frame(matrix(ncol = 6, nrow = 0),stringsAsFactors=FALSE)
colnames(RESULTS) = c("v","metric","model","dropratio","horizon","seed")

for(SEED in 1:MAXSEED){
  for(i in 1:length(DRS)){
    DR = DRS[i]
    for(k in 1:length(MODELS)){
      MODEL = MODELS[k]
      if(MODEL=="DOIRT"){
        file_name = paste("./results/doirt_TAPS_holdout_DR_", DR, "_SEED_", SEED, ".RData", sep="")
      }else if(MODEL=="GPIRT"){
        file_name = paste("./results/gpirt_TAPS_holdout_SEED_DR_", DR, "_", SEED, ".RData", sep="")
      }else if(MODEL=="NIRT"){
        file_name = paste("./results/Nirt_TAPS_holdout_DR_", DR, "_SEED_", SEED, ".RData", sep="")
      }else if(MODEL=="LIRT"){
        file_name = paste("./results/Lirt_TAPS_holdout_DR_", DR, "_SEED_", SEED, ".RData", sep="")
      } else {
        file_name = paste("./results/", MODEL,"_TAPS_holdout_DR_", DR, "_SEED_", SEED, ".RData", sep="")
      }

      if (file.exists(file_name)){
      load(file_name)
      
      TRAIN_ACC[k,i,SEED] = mean(train_acc)
      TRAIN_LL[k,i,SEED] = mean(train_lls)
      for(h in 1:10){
        TEST_ACC[k,i,h,SEED] = mean(test_acc[[h]])
        TEST_LL[k,i,h,SEED] = mean(test_lls[[h]])
      }
      
      RESULTS[nrow(RESULTS)+1,] = list(mean(train_acc),"train_acc",MODEL,DR,0,SEED)
      RESULTS[nrow(RESULTS)+1,] = list(mean(train_lls),"train_lls",MODEL,DR,0,SEED)
      for(h in 1:10){
        RESULTS[nrow(RESULTS)+1,] = list(mean(test_acc[[h]]),"test_acc",MODEL,DR,h,SEED)
        RESULTS[nrow(RESULTS)+1,] = list(mean(test_lls[[h]]),"test_lls",MODEL,DR,h,SEED)
      }
      
      # for(h in 1:10){
      #   RESULTS[nrow(RESULTS)+1,] = list(mean(pred_theta_sd[TRAIN_START_YEAR+h]),"test_sds",MODEL,DR,h,SEED)
      # }
      }
    }
  }
}

# write.csv(RESULTS, "./results/TAPS_holdout_summary.csv")
write.csv(RESULTS, "./results/TAPS_holdout_summary_rebuttal.csv")