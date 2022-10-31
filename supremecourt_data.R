data = read.csv("./data/SCDB_2021_01_justiceCentered_Citation.csv")

data = data[data$term>=2000,]
data = data[, c("caseId", "term", "caseName", "justice",
                "justiceName", "vote", "caseDisposition")]
data = na.omit(data)
data = data[data$vote %in% c(1,2,3,4,5,7), ]

# data$caseId_noyear = as.integer(data.frame(do.call("rbind", strsplit(as.character(data$caseId), "-", fixed = TRUE)))$X2)
# data = data[data$caseId_noyear<=10,]

all_justices = unique(data[,c("justiceName", "justice")])
all_justice_names = as.character(all_justices$justiceName)
all_justice_ids = all_justices$justice
all_case_ids = as.character(unique(data$caseId))

n = length(unique(all_justice_ids))
m = 1
for(year in unique(data$term)){
  if(length(unique(data[data$term==year,"caseId"]))>m){
    m = length(unique(data[data$term==year,"caseId"]))
  }
}
C = 3
horizon = length(unique(data$term))
gpirt_data = array(array(NA, n*m*horizon), c(n,m, horizon))

for (i in 1:n) {
  justice_id = all_justice_ids[i]
  for(h in 1:horizon){
    year = unique(data$term)[h]
    caseIds = as.character(unique(data[data$term==year,"caseId"]))
    for(j in 1:length(caseIds)){
      caseId = caseIds[j]
      tmp = data[data$caseId==caseId & data$justice==justice_id,]
      if(nrow(tmp)){
        vote = tmp$vote
        if(tmp$vote %in% c(1,5)){
          # voted with majority/judgement of court
          gpirt_data[i,j,h] = 3
        }
        if(tmp$vote %in% c(3,4)){
          # concurrence
          gpirt_data[i,j,h] = 2
        }
        if(tmp$vote %in% c(2,7)){
          # dissent
          gpirt_data[i,j,h] = 1
        }
      }
    }
  }
}

# drop unanimous votes
m = 1
for(h in 1:horizon){
  year = unique(data$term)[h]
  caseIds = unique(data[data$term==year,"caseId"])
  mask = rep(0, length(caseIds))
  for(j in 1:length(caseIds)){
    if(length(unique(gpirt_data[!is.na(gpirt_data[,j,h]),j,h]))!=1){
      mask[j] = 1
    }else{
      data = data[data$caseId!=caseIds[j],]
    }
  }
  gpirt_data[,1:sum(mask),h] = gpirt_data[,(1:length(caseIds))[which(mask==1)],h]
  if(sum(mask)>m){
    m = sum(mask)
  }
}
gpirt_data = gpirt_data[,1:m,]