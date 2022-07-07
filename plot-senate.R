library(dplyr)
library(ggplot2)

gpirt_path = "~/Documents/Github/OrdGPIRT"
setwd(gpirt_path)
load(file="./data/senate_data_85.RData")
load(file='./results/gpirt_senate_85.RData')

# plot three ideology scores
# 2001-2003, 2009-2011 , 2019-2021

all_nominate_data = data.frame(matrix(ncol = 6, nrow = 0))
colnames(all_nominate_data) <- c("session", "gpirt", "nominate", "party", "icpsr", "bioname")
for(h in 1:length(session_ids)){
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
  nominate_data$party = party_code
  nominate_data$session = session_id
  nominate_data$icpsr = current_unique_icpsrs
  nominate_data$bioname = bionames
  all_nominate_data = rbind(all_nominate_data, nominate_data)
  # p[[h]] = ggplot(nominate_data, aes(y=x, x=y, colour = factor(party))) +
  #   geom_point(size=2, aes(shape=factor(party))) + 
  #   xlab("NOMINATE Dimension 1 Ideology") + ylab("GPIRT Ideology") + 
  #   labs(colour = "Party")
}

write.csv(all_nominate_data, file="./results/all_nominate_data_85.csv")

# plot six IRFs
# 2001-2003, 2009-2011 , 2019-2021

# re-generate IRFs
# sample_IRFs = vector(mode = "list", length = SAMPLE_ITERS)
# for (iter in 1:SAMPLE_ITERS){
#   # recover fstar
#   sample_IRFs[[iter]] = recover_fstar(BURNOUT_ITERS+iter-1,samples$f[[iter]],data, 
#                                       as.matrix(samples$theta[iter,,]), samples$threshold[iter,])$fstar
# }

# xs = seq(-5,5,0.01)
# gpirt_iccs = array(array(0, length(xs)*m*horizon), 
#                    c(length(xs),m, horizon))
# 
# source("getprob_gpirt.R")
# for (h in 1:horizon) {
#   for (j in 1:m) {
#     IRFs = matrix(0, nrow=SAMPLE_ITERS, ncol=length(xs))
#     for(iter in 1:SAMPLE_ITERS){
#       IRFs[iter, ] = sample_IRFs[[iter]][, j, h]
#     } 
#     probs = getprobs_gpirt(xs, t(IRFs), 
#                            samples$threshold)
#     tmp = probs %>% 
#       group_by(xs) %>%
#       summarize(icc=sum(order*p))
#     gpirt_iccs[,j,h] = tmp$icc
#   }
# }

# To provide additional amounts for low-income home energy assistance
# h = 2, j = 9, rollnumber = 10

# To increase veterans medical care by $410,000,000 in fiscal year 2006
# h = 3, j = 48, rollnumber = 54

# To create a Rural Policing Institute as part of the Federal Law Enforcement Training Center
# h = 4, j = 38, rollnumber = 58

# To specify the criminal offenses that disqualify an applicant from the receipt of a transportation security card
# h = 4, j = 36, rollnumber = 55

# A bill to amend the Fair Labor Standards Act of 1938 to provide for an increase in the Federal minimum wage
# h = 4, j = 28, rollnumber = 42

# A bill to provide greater transparency in the legislative process
# h = 4, j = 15, rollnumber = 19

# To establish a Senate Office of Public Integrity
# h = 4, j = 14, rollnumber = 18

# To prohibit lobbyists and entities that retain or employ lobbyists from throwing lavish parties honoring Members at party conventions.
# h = 4, j = 9, rollnumber = 13

# To prevent government shutdowns
# h = 4, j = 5, rollnumber = 6

# To provide assistance for States with percentages of children with no health insurance coverage above the national average.
# h = 5, j = 30, rollnumber = 30

# To require the use of competitive procedures to award contracts, grants, and cooperative agreements funded under this Act
# h = 5, j = 50, rollnumber = 50

# To prohibit the regulation of greenhouse gases from certain sources
# To suspend, for 2 years, any Environmental Protection Agency enforcement of greenhouse gas regulations, to exempt American agriculture from greenhouse gas regulations, and to increase the number of companies eligible to participate in the successful Advanced Energy Manufacturing Tax Credit Program
# To suspend, until the end of the 2-year period beginning on the date of enactment of this Act, any Environmental Protection Agency action under the Clean Air Act with respect to carbon dioxide or methane pursuant to certain proceedings, other than with respect to motor vehicle emissions.
# h = 6, j = 36-38, rollnumber = 51-53

# A joint resolution making further continuing appropriations for fiscal year 2011, and for other purposes.
# h = 6, j = 23, rollnumber = 29

# A bill to reauthorize the Violence Against Women Act of 1994
# h = 7, j = 18, rollnumber = 19

# A bill to extend the termination date of the Terrorism Insurance Program established under the Terrorism Risk Insurance Act of 2002
# h = 8, j = 2, rollnumber = 2

hs = c(2, 3, 4, 5, 6, 7)
js = c(9, 48, 36, 49, 35, 17)

all_iccs = data.frame(matrix(ncol = 4, nrow = 0))
colnames(all_iccs) <- c("session", "lik", "nominate", "desc")

xs = seq(-5,5,0.01)
idx = 301:701
xs = xs[idx]
  
for(k in 1:length(hs)){
  h = hs[k]
  j = js[k]
  
  # read roll call datga
  session_id = session_ids[h]
  rollcalls = read.csv(paste("./data/S", session_id, "_rollcalls.csv", sep=""))
  
  rollcalls = rollcalls[rollcalls$session==1,
                        c("congress", "rollnumber", "yea_count","nay_count" )]
  
  rollcalls = rollcalls[(rollcalls$yea_count!=0)&(rollcalls$nay_count!=0),]
  rollcalls = rollcalls[1:num_bill_per_session, ]
  rollnumbers = unique(rollcalls$rollnumber)
  
  rollcalls = read.csv(paste("./data/S", session_id, "_rollcalls.csv", sep=""))
  
  # load gpirt icc
  for(idx in 1:length(xs)){
    all_iccs[nrow(all_iccs)+1,] = c(session_id, gpirt_iccs[idx,j, h]-1, xs[idx], 
                     toString(rollcalls[rollcalls$rollnumber==rollnumbers[j], "vote_desc"]))
  }
  
}

write.csv(all_iccs, file="./results/gpirt_icc_result.csv")