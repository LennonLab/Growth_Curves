################################################################################
#                                                                              #
#	Lennon Lab Growth Curve Analysis (with Synergy MX Plate Reader)              #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: M. Larsen (2013/07/18)                                           #
# Modified by: M. Muscarella                                                   #
#                                                                              #
#	Last update: 1/16/13                                                         #
#                                                                              #
################################################################################

# Code Dependencies
require(reshape)
require(pracma)
source("read.synergy.r")
source("curve_fit_fxs.R")
source("grid.mle2.R")


dat<-read.csv("growthcurve_130412.txt",header=TRUE,sep="\t")
time<-dat[1]

#remove temperature column and remake dataframe
dat2<-dat[,3:ncol(dat)]
dat3<-data.frame(time,dat2)
m1<-melt(dat3)
colnames(m1)<-c("time","well","abs")
r<-unique(m1$well)
uids<-unique(r)

# initial guesses of non-grid parameters
intercept.guess<-0.1

# storing umax estimate comparisons
fx.comp<-matrix(NA,nrow=length(unique(r)),ncol=1+1+2+8)
colnames(fx.comp)<-c("model","top.mod","fit1","fit2","ci1 2.5 %","ci1 97.5 %","ci2 2.5 %","ci2 97.5 %","ciFI1 2.5 %","ciFI1 97.5 %","ciFI2 2.5 %","ciFI2 97.5 %")
head(fx.comp)

# initialize data storage
results<-matrix(NA,nrow=length(uids),ncol=13)
colnames(results)<-c("Curve","best.mod","b0","A","umax","L","dd","topt","z","umax.lw","umax.up","umax.lw.FI","umax.up.FI")
results<-as.data.frame(results)

#head(results)

# Shit gets real
uids<-unique(r)
i<-1
for(i in 1:length(uids)){ 



#for(i in 5:length(uids)){ 
	print(paste(i/length(uids),"%"))	# status
	# extract data
	s<-uids[i]
	print(s)
	dat4a<-hampel(m1$abs[m1$well==s],5,t0=1)
	dat4<-dat4a$y

	#time2<-m1$time[m1$well==s]
	time2<-as.numeric(seq(1:73))
	tmpdata<-data.frame(dat4,time2,s)	
	plot(dat4~time2, ylim=c(0,max(tmpdata[,1])+0.2),main=s[1],ylab="ABS",xlab="Time", pch=19,data=tmpdata)
		#A=0.8
		#L=20
		#umax=0.05
		#curve(A*exp(-exp(umax*exp(1)/A*(L-x)+1)),1,72,add=T,col="blue")

# Check for flat oxymoronic curves (?)
if(diff(range(tmpdata$dat4))>0.05){
	# set grid and start lists for each model
	# constant model, dat~dnorm(mean=b0,sd=exp(z))
	grids0<-list()
	start0<-list(b0=intercept.guess,z=-1)
	# old modified gompertz, dat~dnorm(mean=m.gomp(time2,c(b0,A,umax,L)),sd=exp(z))
	grids1<-list(umax=c(0.05,0.1,1),L=c(10,20),z=c(-0.5,-2))
	start1<-list(b0=intercept.guess,A=max(tmpdata[,1]),umax=NA,L=NA,z=NA)
	# new gompertz, dat~dnorm(mean=new.gomp(time2,c(b0,A,umax,L,dd,topt)),sd=exp(z))
	grids2<-list(umax=c(0.05,0.1,1),L=c(10,20),dd=c(0.01,0.1,1),topt=c(20,50),z=c(-0.5,-2))
	start2<-list(b0=intercept.guess,A=log(max(tmpdata[,1])),umax=NA,L=NA,dd=NA,topt=NA,z=NA)
	# grid.mle2 fits performed
	fit0<-grid.mle2(minuslogl=dat4~dnorm(mean=b0,sd=exp(z)),grids=grids0,start=start0,data=tmpdata)
  print("finished fit0")
	fit1<-grid.mle2(minuslogl=dat4~dnorm(mean=m.gomp(time2,c(b0,A,umax,L)),sd=exp(z)),grids=grids1,start=start1,data=tmpdata,method="BFGS")
  print("finished fit1")
	fit2<-grid.mle2(minuslogl=dat4~dnorm(mean=new.gomp(time2,c(b0,A,umax,L,dd,topt)),sd=exp(z)),grids=grids2,start=start2,data=tmpdata, method="BFGS")
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
	pdf(file=paste("./newplots2/testplot",s,".pdf",sep=""))
	plot(dat4~time2, ylim=c(0,max(tmpdata[,1])+0.2),main=s[1],ylab="ABS",xlab="Time", pch=19,data=tmpdata)
	# update plot with curves
	curve(coef(best.f0)[[1]]+0*x,0,max(tmpdata$time2),col='black',add=T)
	curve(m.gomp(x,coef(best.f1)[1:4]),0,max(tmpdata$time2),col='blue',add=T)
	curve(new.gomp(x,coef(best.f2)[1:6]),0,max(tmpdata$time2),col='red',add=T)
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
  }else{
    # generate plot of model fits
    pdf(file=paste("./newplots2/testplot",s,".pdf",sep=""))
    plot(dat4~time2, ylim=c(0,max(tmpdata[,1])+0.2),main=s[1],ylab="ABS",xlab="Time", pch=19,data=tmpdata)
    dev.off()
    results$Curve[i]<-s
    results$best.mod[i]<-0
  }}
results1<-results
write.csv(results,"results.csv")



