gpirt_path = "~/Documents/Github/gpirt"
setwd(gpirt_path)
library(Rcpp)
Rcpp::compileAttributes()
install.packages(gpirt_path, type="source", repos = NULL)# ,lib=R_path, INSTALL_opts = '--no-lock')
setwd("../OrdGPIRT")

library(gpirt)
library(dplyr)
library(ggplot2)
library(stats)
library(haven)
library(ltm)

gpirt_path = "~/Documents/Github/OrdGPIRT"
setwd(gpirt_path)
source("supremecourt_data.R")

theta_os = 1
theta_ls = 7

SEED = 1
SAMPLE_ITERS = 100
BURNOUT_ITERS = 100
THIN = 1
CHAIN = 1
beta_prior_sds =  matrix(3.0, nrow = 3, ncol = ncol(gpirt_data))
theta_prior_sds =  matrix(1.0, nrow = 2, ncol = nrow(gpirt_data))
theta_prior_sds[2,] = 1/horizon
beta_prior_sds[3,] = 0.5
samples_all <- gpirtMCMC(gpirt_data, SAMPLE_ITERS,BURNOUT_ITERS,
                         THIN, CHAIN,
                         beta_prior_sds = beta_prior_sds,
                         theta_prior_sds = theta_prior_sds,
                         theta_os = theta_os, theta_ls = theta_ls, 
                         vote_codes = NULL, thresholds=NULL,
                         SEED=SEED, constant_IRF = 0)

SAMPLE_ITERS = SAMPLE_ITERS/THIN
library(rstan)
sims <- matrix(rnorm((SAMPLE_ITERS)*CHAIN), nrow = SAMPLE_ITERS, ncol = CHAIN)
theta_rhats = matrix(rnorm(n*horizon), nrow = n, ncol = horizon)
for(i in 1:n){
  for(h in 1:horizon){
    for(c in 1:CHAIN){
      sims[,c] = samples_all[[c]]$theta[,i,h]
    }
    theta_rhats[i, h] = Rhat(sims)
  }
}

samples = samples_all[[1]]

# save.image(file='./results/gpirt_SupremeCourt.RData')

xs = seq(-5,5,0.01)
pred_theta = matrix(0, nrow=n, ncol=horizon)
pred_theta_sd = matrix(0, nrow=n, ncol=horizon)
for(i in 1:n){
  for (h in 1:horizon) {
    tmp = samples$theta[-1,i,h]
    pred_theta[i,h] = mean(tmp)
    pred_theta_sd[i,h] = sd(tmp)
  }
}

# for(it in 1:SAMPLE_ITERS){
#   for(h in 1:horizon){
#     for(j in 1:m){
#       samples$f[[it]][,j,h] = samples$f[[it]][,j,h] + samples$beta[[it]][1,j,h] + samples$beta[[it]][2,j,h]*samples$theta[it,,h]
#       samples$fstar[[it]][,j,h] = samples$fstar[[it]][,j,h] + samples$beta[[it]][1,j,h] + samples$beta[[it]][2,j,h]*xs
#     }
#   }
# }

xs = seq(-5,5,0.01)
idx = 301:701
source("getprob_gpirt.R")
gpirt_iccs = array(array(0, length(xs[idx])*m*horizon),
                   c(length(xs[idx]),m, horizon))

for (h in 1:horizon) {
  year = unique(data$term)[h]
  caseIds = as.character(unique(data[data$term==year,"caseId"]))
  for (j in 1:length(caseIds)) {
    IRFs = matrix(0, nrow=SAMPLE_ITERS, ncol=length(idx))
    for(iter in 1:SAMPLE_ITERS){
      IRFs[iter, ] = samples$fstar[[iter]][idx, j, h]
    }
    thresholds = matrix(0,nrow=SAMPLE_ITERS,ncol=C+1)
    for(iter in 1:SAMPLE_ITERS){
      thresholds[iter, ] = samples$threshold[[iter]][j,,h]
    }
    probs = getprobs_gpirt(xs[idx], t(IRFs), thresholds)
    tmp = probs %>%
      group_by(xs) %>%
      summarize(icc=sum(order*p))
    gpirt_iccs[,j,h] = tmp$icc
  }
}

folder_path = "./figures/SupremeCourt/"
dir.create(file.path(folder_path), showWarnings = FALSE)

# 2017 16
HAJUST = c(0.1,1,-0.2,1,0,-0.1,0,0,0)
VAJUST2 = c(0.5,2.2,1.0,2.2,0,-0.8,0,0,0)

for(h in 1:horizon){
  year = unique(data$term)[h]
  caseIds = as.character(unique(data[data$term==year,"caseId"]))
  subfolder = as.character(year)
  dir.create(file.path(folder_path,subfolder), showWarnings = FALSE)
  for(j in 1:length(caseIds)){
    caseName = unique(as.character(data[data$caseId==caseIds[j],"caseName"]))
    x = pred_theta[,h]
    response = gpirt_data[,j,h]
    irf_plot = data.frame(x,response, all_justice_names)
    gpirt_plot = data.frame(xs[idx],gpirt_iccs[,j,h])
    colnames(gpirt_plot) = c("xs","icc")
    VAJUST = 2*(irf_plot$response[!is.na(irf_plot$response)]>=2)-1
    # VAJUST = 2*(irf_plot$response[!is.na(irf_plot$response)]>=3)-1
    ANGLE = -VAJUST*45
    # ANGLE[irf_plot$response[!is.na(irf_plot$response)]==2] = -135
    p = ggplot()+
      geom_point(data = na.omit(irf_plot), aes(x=x,y=response,color=factor(response)),
                 size=4, shape="|") +
      geom_text(data = na.omit(irf_plot), 
                aes(x=x,y=response, label=all_justice_names),
                vjust=1.5*VAJUST+0.5,hjust=0, angle=ANGLE, size=6) +
      scale_color_manual(name='Level',
                         labels=c("Dissent", "Concurrence", "Plurality"),
                         values=c('black', 'blue', 'red'))+
      geom_line(data = gpirt_plot, aes(x=xs,y=icc), size=1)+
      scale_x_continuous(name="x", breaks = seq(-4, 4, by = 1)) +
      scale_y_continuous(name = "Prob of Voting with Majority",
                         breaks = NULL,
                         labels = NULL,
                         limits = c(1,3),
                         sec.axis = sec_axis(~.,
                                             breaks = 1:3,
                                             labels = c("Dissent", "Concurrence", "Plurality")))+
     theme(panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=2),
            legend.position = "none",
            axis.text.y = element_text(size=20,colour = "black"),
            axis.text.x = element_text(size=20,colour = "black"),
            # axis.title.x=element_text(size=20,face="bold",colour = "black"),
            axis.title.x=element_text(size=20),
            axis.title.y=element_text(size=20),
            title = element_text(size=20,colour="black"),
            plot.title = element_text(hjust = 0.5),
           plot.margin=grid::unit(c(1,1,1,1), "mm")) +
      ggtitle(caseName)
    ggsave(filename = paste(folder_path, subfolder, "/",
                            as.character(j), ".png",sep = ""),
           width = 8, height = 5, dpi = 300)
  }
}

all_gpirt_scores = data.frame(matrix(ncol = 7, nrow = 0))
for(i in 1:length(all_justice_ids)){
  for(h in 1:horizon){
    year = unique(data$term)[h]
    if(sum(!is.na(gpirt_data[i,,h]))){
      all_gpirt_scores[nrow(all_gpirt_scores)+1,] = c(year, 
         all_justice_ids[i], pred_theta[i,h], all_justice_names[i], 
         as.integer(0),as.integer(0),pred_theta_sd[i,h])
    }
  }
}

colnames(all_gpirt_scores) <- c("year", "justice_id", "gpirt","name",
                                "x", "y", "sd")

for(i in 1:length(all_justice_ids)){
  mask = !is.na(gpirt_data[i,1,])
  for(j in 1:m){
    mask = mask | !is.na(gpirt_data[i,j,])
  }
  
  loc_x = as.integer(mean(which(mask==1)))
  loc_y = pred_theta[i,loc_x]
  all_gpirt_scores[all_gpirt_scores$justice_id==all_justice_ids[i] &
              all_gpirt_scores$year==unique(data$term)[loc_x], 5]=loc_x
  all_gpirt_scores[all_gpirt_scores$justice_id==all_justice_ids[i] & 
       all_gpirt_scores$year==unique(data$term)[loc_x],6]=loc_y
}

all_gpirt_scores$gpirt = as.numeric(all_gpirt_scores$gpirt)
all_gpirt_scores$x = as.numeric(all_gpirt_scores$x)
all_gpirt_scores$y = as.numeric(all_gpirt_scores$y)
all_gpirt_scores$sd = as.numeric(all_gpirt_scores$sd)

all_gpirt_scores["name"][all_gpirt_scores["name"] == "WHRehnquist"] <- "Rehnquist"
all_gpirt_scores["name"][all_gpirt_scores["name"] == "JPStevens"] <- "Stevens"
all_gpirt_scores["name"][all_gpirt_scores["name"] == "SDOConnor"] <- "OConnor"
all_gpirt_scores["name"][all_gpirt_scores["name"] == "AScalia"] <- "Scalia"
all_gpirt_scores["name"][all_gpirt_scores["name"] == "AMKennedy"] <- "Kennedy"
all_gpirt_scores["name"][all_gpirt_scores["name"] == "DHSouter"] <- "Souter"
all_gpirt_scores["name"][all_gpirt_scores["name"] == "CThomas"] <- "Thomas"
all_gpirt_scores["name"][all_gpirt_scores["name"] == "RBGinsburg"] <- "Ginsburg"

all_gpirt_scores["name"][all_gpirt_scores["name"] == "SGBreyer"] <- "Breyer"
all_gpirt_scores["name"][all_gpirt_scores["name"] == "JGRoberts"] <- "Roberts"
all_gpirt_scores["name"][all_gpirt_scores["name"] == "SAAlito"] <- "Alito"
all_gpirt_scores["name"][all_gpirt_scores["name"] == "SSotomayor"] <- "Sotomayor"
all_gpirt_scores["name"][all_gpirt_scores["name"] == "EKagan"] <- "Kagan"
all_gpirt_scores["name"][all_gpirt_scores["name"] == "NMGorsuch"] <- "Gorsuch"
all_gpirt_scores["name"][all_gpirt_scores["name"] == "BMKavanaugh"] <- "Kavanaugh"
all_gpirt_scores["name"][all_gpirt_scores["name"] == "ACBarrett"] <- "Barrett"

all_gpirt_scores$party = "black"
all_gpirt_scores["party"][all_gpirt_scores["name"] == "Alito"] = "red"
all_gpirt_scores["party"][all_gpirt_scores["name"] == "Gorsuch"] = "red"
all_gpirt_scores["party"][all_gpirt_scores["name"] == "Barrett"] = "red"
all_gpirt_scores["party"][all_gpirt_scores["name"] == "Thomas"] = "red"
all_gpirt_scores["party"][all_gpirt_scores["name"] == "Roberts"] = "red"
all_gpirt_scores["party"][all_gpirt_scores["name"] == "Scalia"] = "red"
all_gpirt_scores["party"][all_gpirt_scores["name"] == "Rehnquist"] = "red"
all_gpirt_scores["party"][all_gpirt_scores["name"] == "Kavanaugh"] = "red"
all_gpirt_scores["party"][all_gpirt_scores["name"] == "Kennedy"] = "red"
all_gpirt_scores["party"][all_gpirt_scores["name"] == "OConnor"] = "red"
all_gpirt_scores["party"][all_gpirt_scores["name"] == "Souter"] = "red"
all_gpirt_scores["party"][all_gpirt_scores["name"] == "Stevens"] = "red"
all_gpirt_scores["party"][all_gpirt_scores["name"] == "Sotomayor"] = "blue"
all_gpirt_scores["party"][all_gpirt_scores["name"] == "Ginsburg"] = "blue"
all_gpirt_scores["party"][all_gpirt_scores["name"] == "Breyer"] = "blue"
all_gpirt_scores["party"][all_gpirt_scores["name"] == "Kagan"] = "blue"

tmp = all_gpirt_scores %>% group_by(name) %>% 
  mutate(meangpirt = mean(gpirt))  %>% 
 arrange(desc(meangpirt))

colfunc1 <- colorRampPalette(c("orange3","darkred"))
colfunc2 <- colorRampPalette(c("darkblue","dodgerblue"))
mypalette = c(colfunc2(6),colfunc1(10))

# all_gpirt_scores[all_gpirt_scores$name!="Barrett",]
tmp = tmp[tmp$name!="Barrett",]
HJUST = c(-2,4.0,5,0.2,0.5,-0.8,1.5,
          -3,-1.7,-2,6,2,-3,0,-1.5)
VJUST = c(-0.5,1.5,-0.3,-0.2,-0.5,1.8,0,
          3.0,1,-0.5,-1.3,1.5,0.5,1.3,0.7)
p = ggplot(tmp, 
           aes(x=year,y=gpirt,group=name,
               color=factor(meangpirt)))+
  scale_y_continuous(name="GD-GPIRT") + 
  scale_x_discrete(breaks=c(2000,2005,2010,2015,2020),labels = c(2000,2005,2010,2015,2020)) + 
  geom_line(size=0.8) +
  geom_text(data=subset(tmp, x!=0),
    aes(x,y,label=name,color=factor(meangpirt)),
    hjust=HJUST, vjust=VJUST, size=6) +
#   scale_color_manual(values = c("aquamarine", "bisque","blue","brown",
#                                "chocolate", "darkgoldenrod", "darkolivegreen","deeppink",
#                                 "firebrick", "gold","grey","khaki",
#                                "orange","sienna","turquoise","yellow")) +
  scale_color_manual(values=mypalette) +
 # scale_linetype_manual(values = c("longdash","solid")) +
  theme(panel.background = element_rect(colour = "white", fill = "white"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "gray"),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size=16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=16))

ggsave(filename = paste(folder_path, "gpirt_supremecourt_2000", ".pdf",sep = ""),
       width = 12, height = 4, dpi = 300)

# for(i in 1:length(all_justice_ids)){
#   id = all_justice_ids[i]
#   name = all_justice_names[i]
#   jpeg(file=paste('./figures/SupremeCourt/', name ,'.jpeg' , sep=''))
#   mask = !is.na(gpirt_data[i,1,])
#   plot((1:horizon)[mask], pred_theta[i,mask])
#   title(name)
#   dev.off()
# }

mqData = read.csv("./data/dynamic_MQ/mqjustices.csv")

ideological_trends = data.frame(matrix(ncol = 4, nrow = 0))
colnames(ideological_trends) <- c("year", "ideology", "avgscore", "sd")

for(year in unique(mqData$term)){
  avgscore = mean(mqData[mqData$term==year & mqData$post_mn<0, "post_mn"])
  ideological_trends[nrow(ideological_trends)+1,] = c(year, "liberal", avgscore)
  avgscore = mean(mqData[mqData$term==year & mqData$post_mn>0, "post_mn"])
  ideological_trends[nrow(ideological_trends)+1,] = c(year, "conservative", avgscore)
}

ideological_trends$avgscore = as.numeric(ideological_trends$avgscore)
ideological_trends$year = as.numeric(ideological_trends$year)

write.csv(ideological_trends, file="./results/mq_dynamic.csv")

mqData = mqData[mqData$term>=2000,]
mqData$x = 0
mqData$y = 0

for(i in 1:length(all_justice_ids)){
  loc_x = all_gpirt_scores[all_gpirt_scores$justice_id==all_justice_ids[i] &
                             all_gpirt_scores$x!=0, "x"]
  loc_y = mqData[mqData$justice==all_justice_ids[i] & 
                   mqData$term==unique(data$term)[loc_x], "post_mn"]
  mqData[mqData$justice==all_justice_ids[i] & 
           mqData$term==unique(data$term)[loc_x],"x"] =loc_x+min(unique(data$term))-1
  mqData[mqData$justice==all_justice_ids[i] & 
           mqData$term==unique(data$term)[loc_x],"y"] =loc_y
}

ggplot(mqData, 
       aes(x=term,y=post_mn,group=factor(justice)))+
  scale_y_continuous(name="Martin-Quinn") + 
  geom_line() +
  geom_text(data=subset(mqData, x!=0),
            aes(x,y,label=justiceName), hjust=1, vjust=0, size=2)

gpirt_supreme_court_scores = data.frame(matrix(ncol = 5, nrow = 0))
colnames(gpirt_supreme_court_scores) <- c("year", "justice_id", "score","name","type")

for(i in 1:length(all_justice_ids)){
  for(h in 1:horizon){
    year = unique(data$term)[h]
    if(sum(!is.na(gpirt_data[i,,h]))){
      gpirt_supreme_court_scores[nrow(gpirt_supreme_court_scores)+1,]  = c(year,
               all_justice_ids[i], pred_theta[i,h], all_justice_names[i], 
               "GPIRT")
      gpirt_supreme_court_scores[nrow(gpirt_supreme_court_scores)+1,] = c(year,
               all_justice_ids[i], mqData$post_mn[mqData$justice==all_justice_ids[i] & 
                                            mqData$term==year], 
              all_justice_names[i], "Martin-Quinn")
    }
  } 
}

gpirt_supreme_court_scores$score = as.numeric(gpirt_supreme_court_scores$score)

write.csv(gpirt_supreme_court_scores, file="./results/gpirt_Supreme_Court_dynamic.csv")

# for(i in 1:length(all_justice_ids)){
#   id = all_justice_ids[i]
#   name = all_justice_names[i]
#   jpeg(file=paste('./figures/SupremeCourt/', name ,'.jpeg' , sep=''))
#   mask = !is.na(gpirt_data[i,2,])
#   tmp = mqData[mqData$justice==all_justice_ids[i],"post_mn"]
#   plot(pred_theta[i,mask],tmp[which(mqData[mqData$justice==all_justice_ids[i],"term"] %in% (1999+(1:horizon)[mask]))])
#   title(name)
#   dev.off()
# }