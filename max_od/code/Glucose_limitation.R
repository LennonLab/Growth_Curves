rm(list = ls())

library("growthcurver")
library('pracma')
library(tidyverse)


setwd("C:/Users/danschw/GitHub/Growth_Curves/max_od/")

d1 <- read.csv("./data/BDM_limit.csv")
d2 <- read.csv("./data/BDM_limit2.csv")
d.wells <-  read.csv("./data/BDM_well_guide.csv")

setwd("C:/Users/danschw/Documents/Github/Growth_Curves/growthcurver/")
# source('find_peaks.r')

d1 <- read.csv("../data/BDM_limit.csv")
d2 <- read.csv("../data/BDM_limit2.csv")
d.wells <-  read.csv("../data/BDM_well_guide.csv")

d.wells$row <- gsub("\\d","",d.wells$Well)
d.wells$col <-gsub("\\D","",d.wells$Well)

#Lets start by plotting the results of each experimet
 #first I'll transform the data to long format to plot with ggplot
d1.long <- 
  gather(d1, key = "well", value = "OD600" , names(d1[,-1]) )
d1.long$exp=1

d2.long <- 
  gather(d2, key = "well", value = "OD600" , names(d2[,-1]) )
d2.long$exp=2

#Now put the 2 experiments together
d.long <-bind_rows(d1.long,d2.long)

d.long$row <- gsub("\\d","",d.long$well)
d.long$col <-gsub("\\D","",d.long$well)
 #add Glucose concentration data
d.long$Gluc_mM <- NA
for (i in 1:nrow(d.long)){
  exp <- d.long$exp[i]
  well <- d.long$well[i]
  tmp <- d.wells[d.wells$exp==exp,]
  tmp <- tmp[tmp$Well == well,]
  d.long$Gluc_mM[i] <- tmp$Gluc_mM
    
}


ggplot(d.long, aes(x=Time, group=interaction(well,exp)))+
  geom_line(aes(y=OD600, color=as.factor(Gluc_mM), linetype=as.factor(exp)),size=1)+
  scale_color_brewer(name="Glucose (mM)",palette =  "Paired")+
  theme_bw()+xlab("Time (hrs)")

ggplot(d.long, aes(x=Time, group=interaction(well,exp)))+
  geom_line(aes(y=OD600, color=as.factor(Gluc_mM)),size=1)+
  scale_color_brewer(name="Glucose (mM)",palette =  "Paired")+
  theme_bw()+xlab("Time (hrs)")+
  facet_grid(.~exp)

ggplot(filter(d.long,exp==2), aes(x=Time, group=interaction(well,exp)))+
  geom_line(aes(y=OD600, color=as.factor(Gluc_mM)),size=1)+
  scale_color_brewer(name="Glucose (mM)",palette =  "Paired")+
  theme_bw()+xlab("Time (hrs)")

# Convert the "time" column from hours to minutes
# d$Time <- 60 * d$Time 

# I want to compare the yields of the different glucose treatments. 
# I already have a working trimming scheme for experiment #1
# The problem is that in experiment #2 after reching the maximum the curves decrease (sporulation?) and then later on rise again (evaporation?) 

# Here is the trimming for exp #1
d.wells$max <- NA
d.wells$maxTime <- NA
d.wells$maxOD <- NA 
for (i in 1:nrow(d.wells)) {
  if (d.wells$exp[i]==1){
    
    d.wells$max[i] <- which.max(d1[,as.character(d.wells$Well[i])])
    d.wells$maxTime[i] <- d1$Time[d.wells$max[i]]
    d.wells$maxOD[i] <- d1[d.wells$max[i],as.character(d.wells$Well[i])]
  }
}

# lets plot that out
ggplot(filter(d.long,exp==1), aes(x=Time))+
  geom_line(aes(y=OD600, color=as.factor(Gluc_mM)),size=1)+
  geom_vline(data = filter(d.wells,exp==1), aes(xintercept=maxTime))+
  geom_hline(data = filter(d.wells,exp==1), aes(yintercept=maxOD))+
  scale_color_brewer(name="Glucose (mM)",palette =  "Paired")+
  theme_bw()+xlab("Time (hrs)")+
  facet_grid(row~col)

# The lines fall on the plot maximum for all but the no-glucose treatment which has no maxima
# applying the same for exp #2:
for (i in 1:nrow(d.wells)) {
  if (d.wells$exp[i]==2){
    
    d.wells$max[i] <- which.max(d2[,as.character(d.wells$Well[i])])
    d.wells$maxTime[i] <- d2$Time[d.wells$max[i]]
    d.wells$maxOD[i] <- d2[d.wells$max[i],as.character(d.wells$Well[i])]
  }
}

# lets plot that out
ggplot(filter(d.long,exp==2), aes(x=Time))+
  geom_line(aes(y=OD600, color=as.factor(Gluc_mM)),size=1)+
  geom_vline(data = filter(d.wells,exp==2), aes(xintercept=maxTime))+
  geom_hline(data = filter(d.wells,exp==2), aes(yintercept=maxOD))+
  scale_color_brewer(name="Glucose (mM)",palette =  "Paired")+
  theme_bw()+xlab("Time (hrs)")+
  facet_grid(row~col)

# It works well for the high glucose wells (15, 10 mM) but not for the ones where the second rise reaches higher OD than the first.
#lets look then just at those low concentrations:
ggplot(filter(d.long,exp==2 & Gluc_mM<10), aes(x=Time, group=interaction(well,exp)))+
  geom_line(aes(y=OD600, color=as.factor(Gluc_mM)),size=1)+
  scale_color_brewer(name="Glucose (mM)",palette =  "Paired")+
  scale_x_continuous(minor_breaks = 0:50)+
  theme_bw()+xlab("Time (hrs)")

# if I trim these at 23 the method above should work (but I can't trim the higher ones too)

for (i in 1:nrow(d.wells)) {
  if (d.wells$exp[i]==2 & d.wells$Gluc_mM[i]<10){
    tmp <- d2[,as.character(d.wells$Well[i])]
    tmp <- tmp[d2$Time<=23]
    d.wells$max[i] <- which.max(tmp)
    d.wells$maxTime[i] <- d2$Time[d.wells$max[i]]
    d.wells$maxOD[i] <- d2[d.wells$max[i],as.character(d.wells$Well[i])]
  }
}

# lets plot that out
ggplot(filter(d.long,exp==2), aes(x=Time))+
  geom_line(aes(y=OD600, color=as.factor(Gluc_mM)),size=1)+
  geom_vline(data = filter(d.wells,exp==2), aes(xintercept=maxTime))+
  geom_hline(data = filter(d.wells,exp==2), aes(yintercept=maxOD))+
  scale_color_brewer(name="Glucose (mM)",palette =  "Paired")+
  theme_bw()+xlab("Time (hrs)")+
  facet_grid(row~col)

#Alright!! now we got the maximum identified for all.
 

#That looks good
#I'll try using also the growthcurver package to estimate yield


# prepare table for results ('SummarizeGrowth' returns 16 values)
gc.fit.data <- as.data.frame(matrix(NA, nrow =nrow(d.wells), ncol=16))
#these are the values returned by 'SummarizeGrowth' in $vals
colnames(gc.fit.data) <- c( "k","k_se","k_p","n0","n0_se","n0_p","r","r_se","r_p","sigma","df","t_mid","t_gen","auc_l","auc_e","note")
# Append results to other data
d.wells <- cbind.data.frame(d.wells,gc.fit.data)
rm(gc.fit.data)
# find min value for each well to mimic the background correction when plotting:
d.wells$min <-NA
for (i in 1:nrow(d.wells)) {
  if (d.wells$exp[i]==1) {
    d.wells$min[i] <- min(d1[,as.character(d.wells$Well[i])])
  }
  if (d.wells$exp[i]==2) {
    d.wells$min[i] <- min(d2[,as.character(d.wells$Well[i])])
  }
}



pdf("./fig/plots.pdf")

for (i in 1:nrow(d.wells)){
  well <- as.character(d.wells$Well[i])
  if (d.wells$exp[i]==1){
    # use crowth curver to fit growth curve on data trimmed by max OD found above
    cur <- SummarizeGrowth(d1$Time, d1[, well], t_trim = d.wells$maxTime[i]*1.01, bg_correct = "min")
    plot(cur, main = paste('exp1  ',well, '  Gluc=',d.wells$Gluc_mM[i],'mM'),
         xlim=c(0,tail(d1$Time,1)),lwd=2)
    points(d1$Time,d1[,well]-d.wells$min[i], type = 'l', col='blue',lwd=2)
    abline(h=cur$vals$k, col='grey',lwd=2)
  } else{
    # use crowth curver to fit growth curve on data trimmed by max OD found above
    cur <- SummarizeGrowth(d2$Time, d2[, well], t_trim = d.wells$maxTime[i]*1.01, bg_correct = "min")
    plot(cur, main = paste('exp1  ',well, '  Gluc=',d.wells$Gluc_mM[i],'mM'),
         xlim=c(0,tail(d2$Time,1)),lwd=2)
    points(d2$Time,d2[,well]-d.wells$min[i], type = 'l', col='blue',lwd=2)
    abline(h=cur$vals$k, col='grey',lwd=2)
  }

  d.wells[i,c( "k","k_se","k_p","n0","n0_se","n0_p","r","r_se","r_p","sigma","df","t_mid","t_gen","auc_l","auc_e","note")] <- unlist(cur$vals)
}
dev.off()


p <- 

ggplot(filter(d.wells, Gluc_mM>0, as.numeric(k)<5 ), aes(x=Gluc_mM, y=as.numeric(k)))+
  # geom_vline(xintercept = 2.73, color="red")+# Dawes chemostat conc. of Glucose
  # geom_errorbar(aes(ymin=k-k_se, ymax=k+k_se, color=rep), width=.1) +
  geom_point(aes(color=as.factor(Gluc_mM)))+ 
  geom_smooth()+
  xlab("Glucose(mM)")+ylab("yield(OD)")+
  scale_color_brewer(name="Glucose (mM)",palette =  "Paired")+
  ggtitle("Glucose limitation in defined media")+
  theme_bw()

ggsave("./fig/growthvurver_yield.pdf",p)



#The growthcurver isn't doing a good job!!!
#I will stick to the maxOD!


#Jay suggested I fit a  Michaelis-Menten curve to this data
#below is an adaptation of his script for this
require("bbmle")
## Run MLE
#starting values 
V = max(d.wells$maxOD) #max rate
K = 5 #conc. of 1/2Vmax
Z = 1 #sd

# Michaelis-Menten
fit <- mle2(d.wells$maxOD ~ dnorm(mean = v * d.wells$Gluc_mM / 
                                           (kmm + d.wells$Gluc_mM), sd = z), start = list(v = V, kmm = K, z = Z), 
            data = d.wells)

# Plot Data

png(filename="./fig/GlucMM.png",
    width = 1200, height = 1200, res = 96*2)

plot.new()
par(mar = c(7, 7, 5, 7))

plot( d.wells$Gluc_mM, d.wells$maxOD, xlim = c(0, 45), 
     ylim = c(0, 1.3), type = "p", 
     pch = 22, bg = "grey", col = "black", cex = 2, ylab = "", xlab = "", 
     cex.lab = 1.5, las = 1, lwd = 2, yaxt = "n", xaxt = "n")
box(lwd=2)

# Add ticks and tick labels
axis(side = 2, lwd.ticks = 2, las = 1, cex.axis = 1.25, 
     labels = c("0.0", "0.25", "0.5", "0.75", "1.0", "1.25"), at = c(0, 0.25, 0.5, 0.75, 1.0, 1.25))

axis(side = 4, labels = F, lwd.ticks = 2, 
     at = c(0, 0.25, 0.5, 0.75, 1.0, 1.25))

axis(side = 1, lwd.ticks = 2, cex.axis = 1.25, las = 1, mgp = c(3, 1, 0),
     labels = seq(0,45,5), at = seq(0,45,5))

axis(side = 3, labels = F, lwd.ticks = 2, las = 1, cex.axis = 1.25, 
     at =seq(0,45,5))

mtext('Glucose (mM)', side = 1, outer = TRUE, cex = 1.5, 
      line = -4, adj = 0.5)

mtext(expression(paste('Max OD')), 
      side = 2, outer = TRUE, cex = 1.5, line = -3, adj = 0.6)

# Plot function
curve((coef(fit)[[1]] * x) / (coef(fit)[[2]] + x), from = min(d.wells$Gluc_mM), to = max(d.wells$Gluc_mM), add = TRUE, lty = 6, lwd = 2.5)
dev.off()
graphics.off()

###########################################
###########################################
###########################################

 # #I'll try and use the lab script
# write_csv(d.trunc,path = "BDM_limit_max.csv")
# 
# # Load Dependencies
# source("../bin/modified_Gomp.R")
# 
# # Run Example (CSV file)
# # CSV file must have the following format: Header followed by whole or decimal numbers 
# # Time  A1  A2  A3  ...
# # 30    0.2 0.1 0.4 ...
# growth.modGomp("BDM_limit_max.csv", "csv_test", synergy=F, temp=F, smooth =T)
# 
# res <- read.csv("../output/csv_test.txt")
# plot(res$L, res$A)

