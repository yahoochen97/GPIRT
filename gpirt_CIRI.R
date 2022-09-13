# gpirt_path = "~/Documents/Github/gpirt"
# setwd(gpirt_path)
# library(Rcpp)
# Rcpp::compileAttributes()
# install.packages(gpirt_path, type="source", repos = NULL)#,lib=R_path, INSTALL_opts = '--no-lock')
# setwd("../OrdGPIRT")

library(gpirt)
library(dplyr)
library(ggplot2)
library(stats)
library(haven)

gpirt_path = "~/Documents/Github/OrdGPIRT"
setwd(gpirt_path)
TYPE = "GP"

data = read.csv("./data/dataverse_files/HumanRightsProtectionScores_v4.01.csv")

DOIRT_theta = read.csv("./data/dataverse_files/DOIRT_PhysintLatentVariableEsimates_1981_2010.csv")

# ordered variable for disappearances/extra-judical/political imprisonment/torture
# in the CIRI dataset
data = data[data$YEAR>=1981 & data$YEAR<=2010,
            c("YEAR", "CIRI", "DISAP", "KILL", 
              "POLPRIS", "TORT", "theta_mean",
              "theta_sd", "country_name")]

unique_sessions = unique(data$YEAR)
horizon = length(unique_sessions)
m = 4 # 4 items

# impute missing value in CIRI by country_name if possible
country_names = unique(data$country_name)
for(country_name in country_names){
  tmp = data[data$country_name==country_name & !is.na(data$CIRI), c("YEAR","CIRI", "country_name")]
  if(nrow(tmp)){
    id = unique(tmp$CIRI)
    tmp2 = data[data$country_name==country_name & is.na(data$CIRI), ]
    if(nrow(tmp2)){
      data[data$country_name==country_name & is.na(data$CIRI), "CIRI"] = id
    }
  }else{
    print(country_name)
  }
}

# remove  German Democratic Republic/Yemen People's Republic/Kosovo/South Sudan
data = data[!is.na(data$CIRI), ]

# Yugoslavia has 4 CIRIs: 689(1981-1991) 563(1992-1999, 2003-2005) 692(2000-2002) 560(2006-2011)
# Czechia has 2 CIRIs: 242(1981-1992) and 239(1993-2011)
# Russia has 2 CIRIs: 590(1981-1991) and 530(1992-2011)
# CIRI 299 represents German Federal Republic(1981-1990) and Germany(1990-2011)

# combine Yugoslavia
data[data$country_name=="Yugoslavia", "CIRI"] = 689

# combine Czechia
data[data$CIRI==239,"CIRI"] = 242

# combine Russian
data[data$CIRI==530, "CIRI"] = 590

# combine Germany
data = data[data$YEAR!=1990 | data$country_name!="German Federal Republic", ]
data[data$country_name=="German Federal Republic", "country_name"] = "Germany"

country_names = unique(data$country_name)
CIRIs = unique(data$CIRI)
n = length(CIRIs)

CIRI_data = array(array(NA, n*m*horizon), c(n, m, horizon))
C = 3

questions = c("DISAP", "KILL", "POLPRIS", "TORT")
for (h in 1:length(unique(data$YEAR))) {
  year = unique(data$YEAR)[h]
  for(j in 1:m){
    for(i in 1:length(CIRIs)){
      tmp = data[data$YEAR==year & data$CIRI==CIRIs[i], c("YEAR", "country_name","theta_mean","theta_sd", questions[j])]
      if(nrow(tmp)){
        CIRI_data[i,j,h] = tmp[[questions[j]]] + 1
      }
    }
  }
}

CIRI_theta = array(array(NA, n*horizon), c(n, horizon))
for (h in 1:length(unique(data$YEAR))) {
  year = unique(data$YEAR)[h]
  for(j in 1:m){
    for(i in 1:length(CIRIs)){
      tmp = DOIRT_theta[DOIRT_theta$YEAR==year & DOIRT_theta$CIRI==CIRIs[i],
                        c("YEAR", "CTRY","latentmean","latentsd", questions[j])]
      if(nrow(tmp)){
        CIRI_theta[i,h] = tmp[["latentmean"]]
      }
    }
  }
}

# CIRI_theta = (CIRI_theta-mean(CIRI_theta,na.rm = TRUE)) / sd(CIRI_theta,na.rm = TRUE)

if(TYPE=="GP"){
  theta_os = 1
  theta_ls = 6
}else if(TYPE=="CST"){
  theta_os = 1
  theta_ls = 10*horizon
}else{
  theta_os = 1
  theta_ls = 0.1
}


theta_init = array(array(0, n*horizon), c(n, horizon))
theta_init[!is.na(CIRI_theta)] = CIRI_theta[!is.na(CIRI_theta)]
for(h in 1:horizon){
  theta_init[,h] = theta_init[,h] + 0.1*rnorm(n)
}

SAMPLE_ITERS = 100
BURNOUT_ITERS = 100
SEED = 1
THIN = 1
CHAIN = 1
beta_prior_sds =  matrix(1.0, nrow = 2, ncol = ncol(CIRI_data))
beta_proposal_sds =  matrix(0.1, nrow = 2, ncol = ncol(CIRI_data))
samples_all <- gpirtMCMC(CIRI_data, SAMPLE_ITERS,BURNOUT_ITERS,
                     THIN, CHAIN, theta_init = theta_init, 
                     beta_prior_sds = beta_prior_sds, 
                     beta_proposal_sds = beta_proposal_sds,
                     theta_os = theta_os, theta_ls = theta_ls, 
                     vote_codes = NULL, thresholds=NULL,
                     SEED=SEED, constant_IRF = 1)

samples = samples_all[[1]]
SAMPLE_ITERS = SAMPLE_ITERS/THIN
# library(rstan)
# sims <- matrix(rnorm((1+SAMPLE_ITERS)*CHAIN), nrow = 1+SAMPLE_ITERS, ncol = CHAIN)
# theta_rhats = matrix(rnorm(n*horizon), nrow = n, ncol = horizon)
# xs = seq(-5,5,0.01)
# idx = 401:601
# irf_rhats = matrix(rnorm(length(idx)*m), nrow = length(idx), ncol = m)
# for(i in 1:n){
#   for(h in 1:horizon){
#     for(c in 1:CHAIN){
#       pred_theta = colMeans(samples_all[[c]]$theta)
#       mask = (!is.na(CIRI_theta[,h]))
#       for(j in 1:m){
#         mask = mask & (!is.na(CIRI_data[,j,h]))
#       }
#       sims[,c] = sign(cor(CIRI_theta[mask,h],pred_theta[mask,h]))*samples_all[[c]]$theta[,i,h]
#     }
#     theta_rhats[i, h] = Rhat(sims)
#   }
# }

xs = seq(-5,5,0.01)
pred_theta = matrix(0, nrow=n, ncol=horizon)
pred_theta_sd = matrix(0, nrow=n, ncol=horizon)
for(i in 1:n){
  for (h in 1:horizon) {
    # mask = rep(1,SAMPLE_ITERS)
    # for(it in 1:SAMPLE_ITERS){
    #   mask[it] = sign(cor(samples$f[[it]][, 1, h],samples$f[[1]][, 1, h]))
    # }
    tmp = samples$theta[-1,i,h]# *mask
    pred_theta[i,h] = mean(tmp)
    pred_theta_sd[i,h] = sd(tmp)
  }
}

# compare with CIRI theta
cor_theta = rep(0, horizon)
for(h in 1:horizon){
  mask = is.na(CIRI_theta[,h])
  for(j in 1:m){
    mask = mask | is.na(CIRI_data[,j,h])
  }
  cor_theta[h] = cor(CIRI_theta[!mask,h],pred_theta[!mask,h])
}

for(it in 1:SAMPLE_ITERS){
  for(h in 1:horizon){
    for(j in 1:m){
      samples$f[[it]][,j,h] = samples$f[[it]][,j,h] + samples$beta[[it]][1,j,h] + samples$beta[[it]][2,j,h]*samples$theta[it,,h]
      samples$fstar[[it]][,j,h] = samples$fstar[[it]][,j,h] + samples$beta[[it]][1,j,h] + samples$beta[[it]][2,j,h]*xs
    }
  }
}

xs = seq(-5,5,0.01)
idx = 221:751
gpirt_iccs = array(array(0, length(xs[idx])*m*horizon),
                   c(length(xs[idx]),m, horizon))

source("getprob_gpirt.R")
for (h in 1:1) {
  for (j in 1:m) {
    IRFs = matrix(0, nrow=SAMPLE_ITERS, ncol=length(idx))
    for(iter in 1:SAMPLE_ITERS){
      # tmp = sign(cor(samples$fstar[[iter]][, j, h],samples$fstar[[1]][, j, h]))
      IRFs[iter, ] = samples$fstar[[iter]][idx, j, h]# *tmp
    }
    probs = getprobs_gpirt(xs[idx], t(IRFs), samples$threshold)
    tmp = probs %>%
      group_by(xs) %>%
      summarize(icc=sum(order*p))
    gpirt_iccs[,j,h] = tmp$icc
    
    q = ggplot(probs, aes(x=xs, y=p, group=order, color=factor(order))) +
      geom_line(size=2) +ggtitle(paste("GP IRT item ",j, sep="")) +
      theme(plot.title = element_text(hjust = 0.5))
    # print(q)
  }
}

h=1
for (j in 1:m) {
  IRFs = matrix(0, nrow=SAMPLE_ITERS, ncol=length(idx))
  for(iter in 1:SAMPLE_ITERS){
    # tmp = sign(cor(samples$f[[iter]][, j, h],samples$f[[1]][, j, h]))
    IRFs[iter, ] = samples$fstar[[iter]][idx, j, h]# *tmp
  }
  probs = getprobs_gpirt(xs[idx], t(IRFs), samples$threshold)
  tmp = probs %>%
    group_by(xs) %>%
    summarize(icc=sum(order*p))
  gpirt_iccs[,j,h] = tmp$icc
  
  q = ggplot(probs, aes(x=xs, y=p, group=order, color=factor(order))) +
    geom_line(size=2) +ggtitle(paste("GP IRT item ",j, sep="")) +
    theme(plot.title = element_text(hjust = 0.5))
  print(q)
}

plot(samples$theta[10,,10],samples$f[[10]][,1,10])

plot(xs[idx],gpirt_iccs[,1,30], ylim=c(1,3))
mask = is.na(CIRI_data[,1,30])
points(pred_theta[!mask,30], CIRI_data[!mask,1,30])

tmp=rep(0,SAMPLE_ITERS)
for(it in 1:SAMPLE_ITERS){
  tmp[it] = samples$beta[[it]][2,1,20]
}
plot(1:SAMPLE_ITERS,tmp)

# plot(xs[idx],gpirt_iccs[,1,1], ylim=c(1,3))
# mask = is.na(CIRI_data[,1,1])
# points(CIRI_theta[!mask,1], CIRI_data[!mask,1,1])


# plot(1:101,samples$ll/sum(!is.na(CIRI_data)))
# plot(1:31,samples$theta[40,1,])
# plot(samples$theta[20,,10],samples$f[[20]][,1,10], ylim=c(-2,2))
# mask = !is.na(CIRI_data[,1,10])
# points(pred_theta[mask,10], 2*(CIRI_data[mask,1,10]-2))
# plot(1:101,samples$theta[,2,2])

library(countrycode)
continents = countrycode(sourcevar = as.character(country_names),
                         origin = "country.name",
                         destination = "continent")
continents[61] = "Europe"
continents[144] = "Asia"

norm_CIRI_theta = (CIRI_theta - mean(CIRI_theta[!is.na(CIRI_theta)])) / sd(CIRI_theta[!is.na(CIRI_theta)])
norm_pred_theta = pred_theta  / sd(pred_theta[!is.na(CIRI_theta)])

all_CIRI_results = data.frame(matrix(ncol = 5, nrow = 0))
colnames(all_CIRI_results) <- c("session"
                                , "gpirt", "CIRI_theta", "country_name", "continent")
for(h in 1:length(unique_sessions)){
  session_id = unique_sessions[h]
  # nominate scores 
  nominate_scores = CIRI_theta[,h]
  mask = (!is.na(nominate_scores))
  for(j in 1:m){
    mask = mask & (!is.na(CIRI_data[,j,h]))
  }
  nominate_data = data.frame(pred_theta[mask,h], nominate_scores[mask])
  colnames(nominate_data) = c("gpirt", "CIRI_theta")
  nominate_data$session = session_id
  nominate_data$continent = continents[mask]
  nominate_data$country_name = as.character(country_names)[mask]
  all_CIRI_results = rbind(all_CIRI_results, nominate_data)
}

write.csv(all_CIRI_results, file="./results/gpirt_CIRI_results.csv")


folder_path = "./figures/CIRI/"
dir.create(file.path(folder_path), showWarnings = FALSE)

for(h in 1:horizon){
  congress = unique_sessions[h]
  rollcall_ids = 1:4
  senator_ids = CIRIs
  subfolder = as.character(unique_sessions[h])
  dir.create(file.path(folder_path,subfolder), showWarnings = FALSE)
  for(j in 1:length(rollcall_ids)){
    x = pred_theta[,h]
    response = CIRI_data[,j,h]
    irf_plot = data.frame(x,response)
    xs = seq(-5,5,0.01)
    idx = 221:751
    gpirt_plot = data.frame(xs[idx],gpirt_iccs[,j,1])
    colnames(gpirt_plot) = c("xs","icc")
    p = ggplot()+
      geom_point(data = na.omit(irf_plot), aes(x=x,y=response,color=factor(response)),
                 size=4, shape="|") +
      scale_color_manual(name='Level',
                         labels=c("Low", 'Medium','High'),
                         values=c( 'red', "blue",'black'))+
      geom_line(data = gpirt_plot, aes(x=xs,y=icc), size=1, )+
      scale_x_continuous(name=bquote(theta), breaks = seq(-3, 3, by = 1)) + 
      scale_y_discrete(name=NULL, limits=c("low","median","high")) +
      theme(panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=2),
            legend.position = "none",
            axis.text.y = element_text(size=16),
            axis.text.x = element_text(size=16))
    
    ggsave(filename = paste(folder_path, subfolder, "/", as.character(j), ".png",sep = ""),width = 4, height = 3, dpi = 300)
    
    # png(paste(folder_path, subfolder, "/", as.character(j), ".png",sep = ""))
    # print(p)
    # dev.off()
  }
}

# save.image(file='./results/gpirt_CIRI.RData')

# library(gganimate)
# animated_data = data.frame(matrix(ncol = 3, nrow = 0))
# colnames(animated_data) = c("it","theta","f")
# for(it in 1:SAMPLE_ITERS){
#   tmp = sign(cor(samples$fstar[[it]][idx, 1, 1],samples$fstar[[1]][idx, 1, 1]))
#   tmp = data.frame(samples$theta[it,,1], samples$f[[it]][,1,1]*tmp)
#   tmp$it = it
#   colnames(tmp) = c("theta","f","it")
#   animated_data = rbind(animated_data, tmp) 
# }
# p = animated_data %>%
#   ggplot(aes(theta,f,frame=it)) +
#   geom_point() +
#   theme_minimal() +
#   ## gganimate functionality starts here
#   transition_states(it) +
#   ggtitle("iteration: {closest_state}") + 
#   ease_aes()
# 
# animate(p, renderer = gifski_renderer())
# anim_save("f1h1.gif")
# 
# animated_data = data.frame(matrix(ncol = 3, nrow = 0))
# colnames(animated_data) = c("it","t","theta")
# for(it in 1:SAMPLE_ITERS){
#   tmp = data.frame(1:horizon, samples$theta[1+it,1,])
#   tmp$it = it
#   colnames(tmp) = c("t","theta","it")
#   animated_data = rbind(animated_data, tmp)
# }
# p = animated_data %>%
#   ggplot(aes(t,theta,frame=it)) +
#   geom_point() +
#   theme_minimal() +
#   ## gganimate functionality starts here
#   transition_states(it) +
#   ggtitle("iteration: {closest_state}") + 
#   ease_aes()
# 
# animate(p, renderer = gifski_renderer())
# anim_save("theta1.gif")
