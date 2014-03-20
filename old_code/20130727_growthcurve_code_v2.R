##################################################
#
#	BioTek Synergy Mx Growth Curve Analysis
#
##################################################
#
#	WRITTEN AND MODIFIED: 2013/07/18 MLL
#
##################################################

#Clear current workspace
rm(list=ls())
getwd()
setwd("/Users/lennonj/Desktop/2013_GrowthCurve")

#convert from matrix form to long form
library(reshape)
library(pracma)
require(bbmle)
require(numDeriv)

source("curve_fit_fxs.R")
source("grid.mle2.R")

#dat<-read.csv("20130727_Rpf_KBS0714.txt",header=TRUE,sep="\t")
dat<-read.table("growthcurve_130412.txt",header=TRUE,sep="\t")

time<-dat[1]

#remove temperature column and remake dataframe
dat2<-dat[,3:ncol(dat)] # just OD data
dat3<-data.frame(time,dat2) # merge time column with OD data

m1<-melt(dat3)

colnames(m1)<-c("time","well","abs") # can change here to be "rep" or just column ID

r<-unique(m1$well) # extract info on number of uniqe elements

uids<-unique(r) # confused, should be (is!) the same as r

# initial guesses of non-grid parameters
intercept.guess<-0.08 # OD at t = 0

# initialize data storage
results<-matrix(NA,nrow=length(uids),ncol=13)
colnames(results)<-c("Curve","best.mod","b0","A","umax","L","dd","topt","z","umax.lw","umax.up","umax.lw.FI","umax.up.FI")
results<-as.data.frame(results)
#head(results)

uids<-unique(r) # why is this in here again?
i<-43 # are we specifying things on a sample by sample basis?  Don't we want to run this for all columns/wells?

#par(mfrow=c(1,2))


		# where does the following for loops start and stop?  Use indentations. Does this go all the way to end?
for(i in i:length(uids)){ 

#for(i in 5:length(uids)){              # I don't understand why were starting with the 5th column of uids?
	print(paste(i/length(uids),"%"))	# status
			
	# extract data
	s<-uids[i]
	print(s)

	dat4a<-hampel(m1$abs[m1$well==s],5,t0=1) # getting rid of outliers
	dat4<-dat4a$y # keeping non-outliers
	#time2<-m1$time[m1$well==s]
	time2<-as.numeric(seq(1:73)) # this needs to be changed for other datasets, right? nrows(dat)?
	tmpdata<-data.frame(dat4,time2,s)
	tmpdata2<-tmpdata[1:which.max(tmpdata[,1]-1),] # don't understand, but gives abs,time (hr), and well

	
	plot(dat4~time2, ylim=c(0,max(tmpdata[,1])+0.2),main=s[1],ylab="ABS",xlab="Time", pch=19,data=tmpdata)
	plot(dat4~time2, ylim=c(0,max(tmpdata2[,1])+0.2),main=s[1],ylab="ABS",xlab="Time", pch=19,data=tmpdata2)
		
	#why are we plotting against tmpdata and tmpdata2? Add comments please	
		
		
		
# Check for flat oxymoronic curves (?)  --> What is "?" for?  If range in abs data is small then what...?
if(diff(range(tmpdata$dat4))>0.05){

	# set grid and start lists for each model
	
	# constant model, dat~dnorm(mean=b0,sd=exp(z))
	grids0<-list()
	start0<-list(b0=intercept.guess,z=-1)
	
	# old modified gompertz, dat~dnorm(mean=m.gomp(time2,c(b0,A,umax,L)),sd=exp(z))
	grids1<-list(umax=c(0.05,0.1,1),L=c(10,20),z=c(-0.5,-2))
	start1<-list(b0=intercept.guess,A=max(tmpdata2[,1]),umax=NA,L=NA,z=NA)

	# grid.mle2 fits performed
	fit0<-grid.mle2(minuslogl=dat4~dnorm(mean=b0,sd=exp(z)),grids=grids0,start=start0,data=tmpdata2)
	
	# Could use some comments in here
	
	
print("finished fit0")
	
	fit1<-grid.mle2(minuslogl=dat4~dnorm(mean=m.gomp(time2,c(b0,A,umax,L)),sd=exp(z)),grids=grids1,start=start1,data=tmpdata2,method="BFGS")
print("finished fit1")

	# isolate best of each class of models
	best.f0<-fit0$res.mod[[1]]
	best.f1<-fit1$res.mod[[which(fit1$res.mat[,'AIC']==min(fit1$res.mat[,'AIC']))[1]]]

# stash results of model

	# model comparison
	aic.tab<-c(AIC(best.f0),AIC(best.f1))
	top.mod<-which(aic.tab==min(aic.tab))-1
	if(length(which(aic.tab<=min(aic.tab)+2))>1){	# check for ties
		top.mod<-3
	}

	# generate plot of model fits
	pdf(file=paste("./newplots2/testplot",s,".pdf",sep=""))

	plot(dat4~time2, ylim=c(0,max(tmpdata2[,1])+0.2),main=s[1],ylab="ABS",xlab="Time", pch=19,data=tmpdata2)

	# update plot with curves
	curve(coef(best.f0)[[1]]+0*x,0,max(tmpdata2$time2),col='black',add=T)
	curve(m.gomp(x,coef(best.f1)[1:4]),0,max(tmpdata2$time2),col='blue',add=T)
	
	dev.off()


	# CI's for the umax term? (for method comparison)
	# best.f0 assumes umax=0
	
	# attempt profiling
	pf1<-profile(best.f1)
	print("finished pf1")

### Save coefficients of best model
	ind<-which(aic.tab==min(aic.tab))	#index of model with lowest AIC
	best.mod<-c(best.f0,best.f1,best.f2)[[ind]] ###I'm getting errors stating that "object 'best.f2' not found"
	cfs<-coef(best.mod) ## also getting errors with this 

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
	plot(dat4~time2, ylim=c(0,max(tmpdata2[,1])+0.2),main=s[1],ylab="ABS",xlab="Time", pch=19,data=tmpdata2)
	dev.off()

	results$Curve[i]<-s
	results$best.mod[i]<-0
}


}


