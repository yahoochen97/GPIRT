h = 1
j = 205

# gp IRFs
xs = seq(-5,5,0.01)
idx = 301:701
gpirt_iccs = rep(0, length(xs[idx]))

source("getprob_gpirt.R")
IRFs = matrix(0, nrow=SAMPLE_ITERS, ncol=length(idx))
for(iter in 1:SAMPLE_ITERS){
  IRFs[iter, ] = sample_IRFs[[iter]][idx, j, h]
}
probs = getprobs_gpirt(xs[idx], t(IRFs),
                       samples$threshold)
tmp = probs %>%
  group_by(xs) %>%
  summarize(icc=sum(order*p))

gpirt_iccs = tmp$icc - 1
plot(xs[idx],gpirt_iccs, ylim=c(0,1))

# party
session_id = session_ids[h]
members = read.csv(paste("./data/S", session_id, "_members.csv", sep=""))
members = members[members$chamber=="Senate",]
# exlude BARKLEY, Dean
members = members[members$icpsr!=40106, ]
# excluse RUSSELL, Richard Brevard, Jr. of GA
members = members[members$icpsr!=8138, ]
current_unique_icpsrs = unique(members$icpsr)
# nominate scores 
nominate_scores = matrix(0, nrow=length(current_unique_icpsrs), ncol=2)
idx = c()
bionames = c()
for(j in 1:length(current_unique_icpsrs)){
  icpsr = current_unique_icpsrs[j]
  # idx = which(icpsr == current_unique_icpsrs)
  nominate_scores[j,1] = members[members$icpsr==icpsr, "nominate_dim1"]
  nominate_scores[j,2] = members[members$icpsr==icpsr, "nominate_dim2"]
  idx = c(idx, which(icpsr==unique_icpsr))
  bionames = c(bionames, toString(members[members$icpsr==icpsr, "bioname"]))
}
current_pred_theta = sign(cor(pred_theta[idx, h],nominate_scores[,1]))*pred_theta[idx, h]
nominate_data = data.frame(current_pred_theta, nominate_scores[,1])
colnames(nominate_data) = c("gpirt", "nominate")
party_code = members$party_code
party_code[(party_code!=200)&(party_code!=100)] = "Independents"
party_code[party_code==100] = "Democrats"
party_code[party_code==200] = "Republicans"

for(j in 1:length(current_unique_icpsrs)){
  i = idx[j]
  if(!is.na(data[i,j,h])){
    color = "blue"
    if(party_code[j]=="Republicans"){
      color = "red"
    }
    points(pred_theta[i,h], data[i,j,h]-1, pch=19, cex=0.25, col=color)
  }
}
