library(haven)
gpirt_path = "~/Documents/Github/OrdGPIRT"
setwd(gpirt_path)

# read data
abortion_data = read_dta("./data/SenatePeriods.dta")
# remove president
abortion_data = abortion_data[abortion_data$name!="REAGAN",]
abortion_data = abortion_data[abortion_data$name!="BUSH",]
abortion_data = abortion_data[abortion_data$name!="CLINTON",]

# get member icpsr for all sessions
session_ids = 92:101
unique_icpsr = c()
for(session_id in session_ids){
  votes = read.csv(paste("./data/S", session_id, "_votes.csv", sep=""))
  unique_icpsr = unique(c(unique_icpsr, unique(votes$icpsr)))
}

# exclude Senator BARKLEY, Dean Independent of Minnesota 
# Served in Senate 2001-2002
# who have not completed a minimum number of rollcall votes.
unique_icpsr = setdiff(unique_icpsr, c(40106))

# excluse RUSSELL, Richard Brevard, Jr. of GA
# served in 92th congress
unique_icpsr = setdiff(unique_icpsr, c(8138))

# number of bills per session
num_bill_per_session = 1202 # 703
data = array(array(NA, length(unique_icpsr)*num_bill_per_session*length(session_ids)), 
                          c(length(unique_icpsr),num_bill_per_session,length(session_ids)))

for(h in 1:length(session_ids)){
  session_id = session_ids[h]
  rollcalls = read.csv(paste("./data/S", session_id, "_rollcalls.csv", sep=""))
  rollcalls = rollcalls[,c("congress", "rollnumber", "yea_count","nay_count","date" )]
  rollcalls = rollcalls[(rollcalls$yea_count!=0)&(rollcalls$nay_count!=0),]
  
  abortion_rollcalls = unique(abortion_data[abortion_data$congress==session_id, c("rollcall")])
  rollcalls = rollcalls[!(rollcalls$rollnumber %in% abortion_rollcalls),]
  
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
        if(cast_code>=1 & cast_code<=3){
          data[i, j, h] = 2
        }else if(cast_code<=6 & cast_code>=4){
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

save(data,unique_icpsr,session_ids, file="./data/senate_data_92.RData")
