################################################################################
#                                                                              #
#	Lennon Lab Growth Curve Analysis (with Synergy MX Plate Reader)              #
#   Parameter Estimate Code                                                    #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: M. Muscarella                                                    #
#   Based on growthcurve_code.R Written by: M. Larsen (2013/07/18)             #
#                                                                              #
#	Last update: 3/3/14                                                         #
#                                                                              #
################################################################################

growth.modGomp <- function(input=" ", output=" ", intercept.guess=0.1, synergy=T, skip = ""){
# Input = Raw txt output file from Synergy MX
# Intercept.guess = initial guess of non-grid parameter for y intercept '
# synergy=T --> the data comes from the synergy mx machine.
#   If false data should be deliminated in proper format
# Load Code Dependencies
  source("../bin/read.synergy.R")
  source("../bin/curve_fit_fxs.R")
  source("../bin/grid.mle2.R")

# Data Input
  data.in = read.csv("ControlRpf.csv")
  #temp.test <- lm(round(Temp, 2) ~ Time, data=data.in)
  #p <- round(anova(temp.test)$'Pr(>F)'[1], 3)
  #temp.min <- min(data.in$Temp)
  #temp.max <- max(data.in$Temp)
  #temp.diff <- temp.max - temp.min
  #if (temp.diff < 3) {} else {stop("Stop, check for temperature effects")}
  #samples <- colnames(data.in[3:dim(data.in)[2]])

  # Initialize Data Storage
  results <- matrix(NA,nrow=(dim(data.in)[2])-2,ncol=10)
  colnames(results) <- c("Curve","b0","A","umax","L","z",
    "umax.lw","umax.up","umax.lw.FI","umax.up.FI")
  results <- as.data.frame(results)

  # Creat Output
  outfile <- paste("../output/",output,".txt", sep="")
  titles <- c("Curve","b0","A","umax","L","z", "umax.lw","umax.up","umax.lw.FI","umax.up.FI")
  write.table(as.matrix(t(titles)), file=outfile, append=T, row.names=F, col.names=F, sep=",", quote=FALSE)

  for(i in 3:dim(data.in)[2]){
    # Print Operation Status
    print(paste(round(((i-2)/(dim(data.in)[2]-2)*100),2),"% complete", sep = ""), quote=F)

    # Extract Data
    t <- data.in$Time
    s <- data.in[,i]
    if (max(s) - min(s) < 0.05) next
    else realdata <- data.frame(t,s)

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
	 grids1<-list(umax=c(0.05,0.1,1),L=c(-5,-0.5,0.1,5,10,20),z=c(-2,-0.5))
	 start1<-list(b0=intercept.guess,A=max(tmpdata[,2]),umax=NA,L=NA,z=NA)

  # Perform grid.mle2 Fits
    fit1<-grid.mle2(minuslogl=s.trim~dnorm(mean=m.gomp(t.trim,c(b0,A,umax,L)),sd=exp(z)),grids=grids1,start=start1,data=tmpdata,method="BFGS")
    print("finished fit")

	# isolate best of each class of model
    best.f1<-fit1$res.mod[[which(fit1$res.mat[,'AIC']==min(fit1$res.mat[,'AIC']))[1]]]

	# generate plot of model fits
    plot(s ~ t, main=colnames(data.in[i]), ylab="ABS", xlab="Time", pch=19, data=realdata)
    curve(m.gomp(x,coef(best.f1)[1:4]),0,max(realdata$t),col='blue',lwd=2,add=T)
    pdf(file=paste("../output/testplot",colnames(data.in[i]),".pdf",sep=""))
    plot(s ~ t, main=colnames(data.in[i]), ylab="ABS", xlab="Time", pch=19, data=realdata)
    curve(m.gomp(x,coef(best.f1)[1:4]),0,max(tmpdata$t.trim),col='blue',lwd=2,add=T)
    dev.off()

	# attempt profiling
    pf1<-profile(best.f1)
    print("finished pf1")

    if(class(pf1)=="profile.mle2"){
    ci1<-confint(pf1)['umax',]
  	} else {
      ci1<-c(NA,NA)
		  names(ci1)<-c("2.5%","97.5%")
      }
    ciFI1<-confint.FI(best.f1)['umax',]

    # Save coefficients of model
    cfs<-coef(best.f1)
    results$Curve[i] <- colnames(data.in)[i]
    results$b0[i]<-cfs['b0']
    results$A[i]<-cfs['A']
    results$umax[i]<-cfs['umax']
    results$L[i]<-cfs['L']
    results$z[i]<-cfs['z']
    results[i,c("umax.lw" , "umax.up" , "umax.lw.FI" , "umax.up.FI")] <-c (ci1,ciFI1)
    write.table(results[i,], file=outfile, append=T, row.names=F, col.names=F, sep=",", quote=FALSE)
    }
  results1<-results
  write.csv(results,"results.csv")
}
