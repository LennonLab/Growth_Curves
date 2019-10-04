################################################################################
#                                                                              #
#	Lennon Lab Growth Curve Analysis (with Synergy MX Plate Reader)              #
#   Parameter Estimate Code 
# test                                                   #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: M. Muscarella                                                    #
#   Based on growthcurve_code.R Written by: M. Larsen (2013/07/18)             #
#                                                                              #
#	Last update: 2/19/14                                                         #
#                                                                              #
################################################################################

# input = "GrowthCurve_Example.txt"

# This function is currently not working. It is supposed to compare different 
# models but things are complicated. I recommend using the second function for 
# now.

growth.estimate <- function(input=" ", intercept.guess=0.1){
# Input = Raw txt output file from Synergy MX
# Intercept.guess = initial guess of non-grid parameter for y intercept
# Load Code Dependencies
  source("./scripts/read.synergy.R")
  source("./scripts/curve_fit_fxs.R")
  source("./scripts/grid.mle2.R")
  
# Data Input 
  data.in <- read.synergy(input)
  temp.test <- lm(round(Temp, 2) ~ Time, data=data.in)  
  p <- round(anova(temp.test)$'Pr(>F)'[1], 3)
  temp.min <- min(data.in$Temp)
  temp.max <- max(data.in$Temp)
  temp.diff <- temp.max - temp.min
  if (temp.diff < 2) {} else {stop("Stop, check for temperature effects")}
  samples <- colnames(data.in[3:dim(data.in)[2]])
  
# storing umax estimate comparisons
  fx.comp <- matrix(NA,nrow=(dim(data.in)[2])-2,ncol=1+1+2+8)  # Why 1+1+2+8 ????
  colnames(fx.comp)<-c("model","top.mod","fit1","fit2","ci1 2.5 %","ci1 97.5 %",
    "ci2 2.5 %","ci2 97.5 %","ciFI1 2.5 %","ciFI1 97.5 %","ciFI2 2.5 %",
    "ciFI2 97.5 %") 
     
# initialize data storage
  results <- matrix(NA,nrow=(dim(data.in)[2])-2,ncol=13)
  colnames(results)<-c("Curve","best.mod","b0","A","umax","L","dd","topt","z",
    "umax.lw","umax.up","umax.lw.FI","umax.up.FI")
  results<-as.data.frame(results)  
  for(i in 3:dim(data.in)[2]){
  
# Print Operation Status
    print(paste(round(((i-2)/(dim(data.in)[2]-2)*100),2),"% complete", sep = ""), quote=F)
    
# Extract Data
    t <- data.in$Time
    s <- data.in[,i]
    
# Smoothing Function
    s.2 <- as.numeric(filter(s, rep(1/11,11), circular=F, sides=2))
    s.2[1:5] <- s[1:5]
    s.2[(length(s.2)-5):length(s.2)] <- s[(length(s)-5):length(s)]
    s.max <- max(which(s.2 == max(s.2, na.rm=T)))
    t.end <- round(t[s.max],0) + 2
    t.trim <- t[which(t <= t.end)]
    s.trim <- s.2[which(t <= t.end)]
    tmpdata <- data.frame(t.trim, s.trim)
    plot(s.trim ~ t.trim, main=colnames(data.in[i]), ylab="ABS", xlab="Time", pch=19, data=tmpdata)
    
# Set Grid and Start Lists for Each Model
	# Constant Model, dat~dnorm(mean=b0,sd=exp(z))
	grids0<-list()
	start0<-list(b0=intercept.guess,z=-1)
	# Old Modified Gompertz, dat~dnorm(mean=m.gomp(time2,c(b0,A,umax,L)),sd=exp(z))
	grids1<-list(umax=c(0.05,0.1,1),L=c(0.1,5,10,20),z=c(-0.5,-2))
	start1<-list(b0=intercept.guess,A=max(tmpdata[,2]),umax=NA,L=NA,z=NA)
	# New Gompertz, dat~dnorm(mean=new.gomp(time2,c(b0,A,umax,L,dd,topt)),sd=exp(z))
	grids2<-list(umax=c(0.05,0.1,1),L=c(0.1,5,10,20),dd=c(0.01,0.1,1,2),topt=c(0.1,5,10,20,50),z=c(-0.5,-2))
	start2<-list(b0=intercept.guess,A=log(max(tmpdata[,2])),umax=NA,L=NA,dd=NA,topt=NA,z=NA)
	
# Perform grid.mle2 Fits
	fit0<-grid.mle2(minuslogl=s.trim~dnorm(mean=b0,sd=exp(z)),grids=grids0,start=start0,data=tmpdata)
  print("finished fit0")
	fit1<-grid.mle2(minuslogl=s.trim~dnorm(mean=m.gomp(t.trim,c(b0,A,umax,L)),sd=exp(z)),grids=grids1,start=start1,data=tmpdata,method="BFGS")
  print("finished fit1")
	fit2<-grid.mle2(minuslogl=s.trim~dnorm(mean=new.gomp(t.trim,c(b0,A,umax,L,dd,topt)),sd=exp(z)),grids=grids2,start=start2,data=tmpdata, method="BFGS")
  print("finished fit2")	

	# isolate best of each class of models
	best.f0<-fit0$res.mod[[1]]
	best.f1<-fit1$res.mod[[which(fit1$res.mat[,'AIC']==min(fit1$res.mat[,'AIC']))[1]]]
	best.f2<-fit2$res.mod[[which(fit2$res.mat[,'AIC']==min(fit2$res.mat[,'AIC']))[1]]]

	# stash results of model
	# model comparison
	aic.tab<-c(AIC(best.f0),AIC(best.f1),AIC(best.f2))
	top.mod<-which(aic.tab==min(aic.tab))-1
	if(length(which(aic.tab<=min(aic.tab)+2))>1){	# check for ties
		top.mod<-3
	}

	# generate plot of model fits
	pdf(file=paste("./plots/testplot",s,".pdf",sep=""))
	#plot(s.trim~t.trim, ylim=c(0,max(tmpdata[,2])+(0.1*max(tmpdata[,2]))),main=s[1],ylab="ABS",xlab="Time", pch=19,data=tmpdata)
	# update plot with curves
	plot(s.trim ~ t.trim, main=colnames(data.in[i]), ylab="ABS", xlab="Time", pch=19, data=tmpdata)
	curve(coef(best.f0)[[1]]+0*x,0,max(tmpdata$t.trim),col='black',add=T)
	curve(m.gomp(x,coef(best.f1)[1:4]),0,max(tmpdata$t.trim),col='blue',add=T)
	curve(new.gomp(x,coef(best.f2)[1:6]),0,max(tmpdata$t.trim),col='red',add=T)
	dev.off()

	# CI's for the umax term? (for method comparison)
	# best.f0 assumes umax=0
	# attempt profiling

	pf1<-profile(best.f1)
	print("finished pf1")
	pf2<-profile(best.f2)
	print("finished pf2")

  if(class(pf1)=="profile.mle2"){
    ci1<-confint(pf1)['umax',]
  	}else{
      ci1<-c(NA,NA)
		  names(ci1)<-c("2.5 %","97.5 %")
      }
    ciFI1<-confint.FI(best.f1)['umax',]
    if(class(pf2)=="profile.mle2"){
      ci2<-confint(pf2)['umax',]
      }else{
        ci2<-c(NA,NA)
        names(ci2)<-c("2.5 %","97.5 %")
        }
    ciFI2<-confint.FI(best.f2)['umax',]	
    # stash results
    fx.comp[i,1]<-s
    fx.comp[i,2]<-top.mod
    fx.comp[i,3]<-coef(best.f1)['umax']
    fx.comp[i,4]<-coef(best.f2)['umax']
    fx.comp[i,5:ncol(fx.comp)]<-c(ci1,ci2,ciFI1,ciFI2)	
    ### Save coefficients of best model
    ind<-which(aic.tab==min(aic.tab))	#index of model with lowest AIC
    best.mod<-c(best.f0,best.f1,best.f2)[[ind]]
    cfs<-coef(best.mod)
    results$Curve[i]<-s
    results$best.mod[i]<-ind-1	
    results$b0[i]<-cfs['b0']
    results$A[i]<-cfs['A']
    results$umax[i]<-cfs['umax']
    results$L[i]<-cfs['L']
    results$dd[i]<-cfs['dd']
    results$topt[i]<-cfs['topt']
    results$z[i]<-cfs['z']
    if(ind==3){
      results[i,c("umax.lw","umax.up","umax.lw.FI","umax.up.FI")]<-c(ci2,ciFI2)
      }else{
        if(ind==2){
          results[i,c("umax.lw","umax.up","umax.lw.FI","umax.up.FI")]<-c(ci1,ciFI1)				
        }}
#  } else {
#    # generate plot of model fits
#    pdf(file=paste("./newplots2/testplot",s,".pdf",sep=""))
#    plot(s.trim ~ t.trim, main=colnames(data.in[i]), ylab="ABS", xlab="Time", pch=19, data=tmpdata)
#    dev.off()
#    results$Curve[i]<-s
#    results$best.mod[i]<-0
  }
  results1<-results
  write.csv(results,"results.csv")
  }

################################################################################
################################################################################

growth.modGomp <- function(input=" ", intercept.guess=0.1){
# Input = Raw txt output file from Synergy MX
# Intercept.guess = initial guess of non-grid parameter for y intercept
# Load Code Dependencies
  source("read.synergy.r")
  source("curve_fit_fxs.R")
  source("grid.mle2.R")
  
# Data Input 
  data.in <- read.synergy(input)
  temp.test <- lm(round(Temp, 2) ~ Time, data=data.in)  
  p <- round(anova(temp.test)$'Pr(>F)'[1], 3)
  temp.min <- min(data.in$Temp)
  temp.max <- max(data.in$Temp)
  temp.diff <- temp.max - temp.min
  if (temp.diff < 2) {} else {stop("Stop, check for temperature effects")}
  samples <- colnames(data.in[3:dim(data.in)[2]])
     
# Initialize Data Storage
  results <- matrix(NA,nrow=(dim(data.in)[2])-2,ncol=10)
  colnames(results) < -c("Curve","b0","A","umax","L","z",
    "umax.lw","umax.up","umax.lw.FI","umax.up.FI")
  results <- as.data.frame(results)  
  
  for(i in 3:dim(data.in)[2]){
    # Print Operation Status
    print(paste(round(((i-2)/(dim(data.in)[2]-2)*100),2),"% complete", sep = ""), quote=F)
    
    # Extract Data
    t <- data.in$Time
    s <- data.in[,i]
    realdata <- data.frame(t,s)
    
    # Smoothing Function
    s.2 <- as.numeric(filter(s, rep(1/11,11), circular=F, sides=2))
    s.2[1:5] <- s[1:5]
    s.2[(length(s.2)-5):length(s.2)] <- s[(length(s)-5):length(s)]
    s.max <- max(which(s.2 == max(s.2, na.rm=T)))
    t.end <- round(t[s.max],0) + 1
    t.trim <- t[which(t <= t.end)]
    s.trim <- s.2[which(t <= t.end)]
    tmpdata <- data.frame(t.trim, s.trim)
    #plot(s.trim ~ t.trim, main=colnames(data.in[i]), ylab="ABS", xlab="Time", pch=19, data=tmpdata)
    
  # Set Grid and Start Lists for Model
	 # Modified Gompertz, dat~dnorm(mean=m.gomp(time2,c(b0,A,umax,L)),sd=exp(z))
	 grids1<-list(umax=c(0.05,0.1,1),L=c(-5,-0.5,0.1,5,10,20),z=c(-0.5,-2))
	 start1<-list(b0=intercept.guess,A=max(tmpdata[,2]),umax=NA,L=NA,z=NA)

  # Perform grid.mle2 Fits
    fit1<-grid.mle2(minuslogl=s.trim~dnorm(mean=m.gomp(t.trim,c(b0,A,umax,L)),sd=exp(z)),grids=grids1,start=start1,data=tmpdata,method="BFGS")
    print("finished fit")	

	# isolate best of each class of model
    best.f1<-fit1$res.mod[[which(fit1$res.mat[,'AIC']==min(fit1$res.mat[,'AIC']))[1]]]

	# generate plot of model fits
    plot(s ~ t, main=colnames(data.in[i]), ylab="ABS", xlab="Time", pch=19, data=realdata)
    curve(m.gomp(x,coef(best.f1)[1:4]),0,max(realdata$t),col='blue',lwd=2,add=T)
    pdf(file=paste("./plot.out/testplot",colnames(data.in[i]),".pdf",sep=""))
    plot(s ~ t, main=colnames(data.in[i]), ylab="ABS", xlab="Time", pch=19, data=realdata)
    curve(m.gomp(x,coef(best.f1)[1:4]),0,max(tmpdata$t.trim),col='blue',lwd=2,add=T)
    dev.off()

### 2/19/14 Stopping Point: Added modGomp only function. Edited plotting and grids
### The next step is going to be working on stats for curve fits and summary stats
### as well as outputting curve parameters. But good so far.

	# attempt profiling
	pf1<-profile(best.f1)
	print("finished pf1")

  if(class(pf1)=="profile.mle2"){
    ci1<-confint(pf1)['umax',]
  	}else{
      ci1<-c(NA,NA)
		  names(ci1)<-c("2.5 %","97.5 %")
      }
    ciFI1<-confint.FI(best.f1)['umax',]

    # stash results
    fx.comp[i,1]<-s
    fx.comp[i,2]<-coef(best.f1)['umax']
    fx.comp[i,3:ncol(fx.comp)]<-c(ci1,ciFI1)	
    ### Save coefficients of best model
    ind<-which(aic.tab==min(aic.tab))	#index of model with lowest AIC
    best.mod<-c(best.f0,best.f1,best.f2)[[ind]]
    cfs<-coef(best.mod)
    results$Curve[i]<-s
    results$best.mod[i]<-ind-1	
    results$b0[i]<-cfs['b0']
    results$A[i]<-cfs['A']
    results$umax[i]<-cfs['umax']
    results$L[i]<-cfs['L']
    results$dd[i]<-cfs['dd']
    results$topt[i]<-cfs['topt']
    results$z[i]<-cfs['z']
    if(ind==3){
      results[i,c("umax.lw","umax.up","umax.lw.FI","umax.up.FI")]<-c(ci2,ciFI2)
      }else{
        if(ind==2){
          results[i,c("umax.lw","umax.up","umax.lw.FI","umax.up.FI")]<-c(ci1,ciFI1)				
        }}
  }else{
    # generate plot of model fits
    pdf(file=paste("./newplots2/testplot",s,".pdf",sep=""))
    plot(s.trim ~ t.trim, main=colnames(data.in[i]), ylab="ABS", xlab="Time", pch=19, data=tmpdata)
    dev.off()
    results$Curve[i]<-s
    results$best.mod[i]<-0
  }}
results1<-results
write.csv(results,"results.csv")
}
