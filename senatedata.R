gpirt_path = "~/Documents/Github/OrdGPIRT"
gpirt_path = "../gpirt"
setwd(gpirt_path)

# get member icpsr for all sessions
session_ids = 107:116
unique_icpsr = c()
for(session_id in session_ids){
  votes = read.csv(paste("./data/S", session_id, "_votes.csv", sep=""))
  unique_icpsr = unique(c(unique_icpsr, unique(votes$icpsr)))
}

# exclude Senator BARKLEY, Dean Independent of Minnesota 
# Served in Senate 2001-2002
# who have not completed a minimum number of rollcall votes.
unique_icpsr = setdiff(unique_icpsr, c(40106))

# number of bills per session
num_bill_per_session = 703
data = array(array(NA, length(unique_icpsr)*num_bill_per_session*length(session_ids)), 
                          c(length(unique_icpsr),num_bill_per_session,length(session_ids)))

for(h in 1:length(session_ids)){
  session_id = session_ids[h]
  rollcalls = read.csv(paste("./data/S", session_id, "_rollcalls.csv", sep=""))
  
  rollcalls = rollcalls[,c("congress", "rollnumber", "yea_count","nay_count","date" )]
  
  rollcalls = rollcalls[(rollcalls$yea_count!=0)&(rollcalls$nay_count!=0),]
  print(length(unique(rollcalls$rollnumber)))
  # rollcalls = rollcalls[1:num_bill_per_session, ]
  rollnumbers = unique(rollcalls$rollnumber)
  votes = read.csv(paste("./data/S", session_id, "_votes.csv", sep=""))
  
  for(j in 1:length(rollnumbers)){
    rollnumber = rollnumbers[j]
    for(i in 1:length(unique_icpsr)){
      icpsr = unique_icpsr[i]
      cast_code = votes[(votes$icpsr==icpsr)&(votes$rollnumber==rollnumber),]
      
      if(nrow(cast_code)!=0){
        cast_code = cast_code$cast_code
        if(cast_code==1){
          data[i, j, h] = 2
        }else if(cast_code==6){
          data[i, j, h] = 1
        }
      }
    }
  }
}

n = nrow(data)
m = ncol(data)
horizon = length(data[1,1,])
C = 2
train_ratio = 0.8

# split into train and test data
N = n*m*horizon
na_mask = matrix(is.na(data), nrow = N)
train_idx = rep(0,N)
train_idx[sample((1:N)[!na_mask],as.integer(train_ratio*sum(na_mask==0)),replace=FALSE)] = 1
train_idx = array(train_idx, c(n,m,horizon))
data_train = data

for (i in 1:nrow(data)) {
  for (j in 1:ncol(data)) {
    for (h in 1:horizon) {
      if(train_idx[i,j,h]==0){
        data_train[i,j,h] = NA
      }
    }
  }
}

save(data,data_train,train_idx,unique_icpsr,session_ids, file="./data/senate_data.RData")
