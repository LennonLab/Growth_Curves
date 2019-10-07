################################################################################
#                                                                              #
#	Lennon Lab Growth Curve Analysis (using Synergy MX Plate Reader)             #
#   Uses OD data to fit Gomphertz model and estimate parameters                #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: M. Muscarella                                                    #
#   Partially based on growthcurve_code.R Written by: M. Larsen (2013/07/18)   # 
#                                                                              #
# Contributors: M. Muscarella, M. Larsen, J. Lennon V. Kuo, M. Behringer,      #
#               D. Schwartz                                                    #
#                                                                              #
#	Last update: 10/4/2019 by M. Muscarella                                      #
#                                                                              #
################################################################################

growth.modGomp <- function(input = " ", output.name = " ", skip = "", 
                           output.dir = "../output/",
                           intercept.guess = 0.1, delta = 0.05,
                           synergy = T, temp = T, smooth = T, trim = T){
    # Input = Raw txt output file from Synergy MX
    # Intercept.guess = initial guess for y intercept
    # delta = minimum change in OD required for analysis
    # synergy = TRUE if the data comes from the synergy mx machine
    #   If FALSE data should be deliminated in proper format
    # skip = the number of lines to initially skip when importing synergy data
    # temp = TRUE if there temperature data include for QC purposes
    # trim = TRUE if you want to trim the data

  # Load Code Dependencies
  source("../bin/read.synergy.R")
  source("../bin/curve_fit_fxs.R")
  source("../bin/grid.mle2.R")
  
  # Create Directory For Output
  dir.create(output.dir, showWarnings = FALSE)
  
  # Data Input
  if (synergy == T){
    data.in <- read.synergy(input, skip = skip)
    samples <- colnames(data.in[3:dim(data.in)[2]])
  } else {
    if(is.character(input)){
      data.in <- read.csv(input, header=T)
      samples <- colnames(data.in[2:dim(data.in)[2]])
    }else{
      data.in <- input
      samples <- colnames(data.in[2:dim(data.in)[2]])
    }
  }
  
  # Check for Temp Issues
  if (temp == T){
    temp.test <- lm(round(Temp, 2) ~ Time, data = data.in)
    p <- round(anova(temp.test)$'Pr(>F)'[1], 3)
    temp.min <- min(data.in$Temp)
    temp.max <- max(data.in$Temp)
    temp.diff <- temp.max - temp.min
    if (temp.diff < 3) {} else {stop("Stop, check for temperature effects")}
  } else {}

  # Initialize Data Storage
  results <- matrix(NA, nrow=length(samples), ncol=12)
  colnames(results) <- c("Curve", "b0", "A", "umax", "L", "z", "K", "lag", 
                         "umax.lw", "umax.up", "umax.lw.FI", "umax.up.FI")
  results <- as.data.frame(results)

  # Initialize Output
  outfile <- paste(output.dir, paste(output.name, ".txt", sep=""), sep ="/")
  titles <- c("Curve", "b0", "A", "umax", "L", "z", "K", "lag", "umax.lw", "umax.up",
              "umax.lw.FI", "umax.up.FI")
  write.table(as.matrix(t(titles)), file=outfile, append=F, row.names=F,
              col.names=F, sep=",", quote=FALSE)
  
  # Initialize Plotting Device
  outplot <- paste(output.dir, paste(output.name, ".pdf", sep=""), sep ="/")
  pdf(file = outplot, width = 6, height = 4)

  for(i in 1:length(samples)){
    # Print Operation Status
    print(paste("Starting Sample ", samples[i], " (", i, " of ", 
                length(samples), ")", sep = ""), quote=F)

    # Extract Data
    t <- data.in$Time
    s <- data.in[,which(colnames(data.in) == samples[i])]
    if (max(s) - min(s) < delta) {
      plot(s ~ t, main=samples[i], 
           ylab=expression(paste("Absorbance"[600])), 
           xlab=expression(paste("Time (hrs)")),
           pch=19, las = 1, cex.axis = 1, cex.lab = 1.5)
      legend("topleft", legend = bquote(Delta ~ "OD <" ~ .(delta)), bty = "n")
      print(paste("Observed change in OD is not greater than ", delta,
                  " in sample ", samples[i], sep=""), quote=F)
      next
    } else {
      realdata <- data.frame(t,s)
    }

    # Smoothing Function
    if (smooth == T){
      s.1 <- realdata$s
      t.1 <- realdata$t
      s.2 <- as.numeric(stats::filter(s.1, rep(1/11,11), circular=F, sides=2))
      s.2[1:5] <- s.1[1:5]
      s.2[(length(s.2) - 5):length(s.2)] <- s.1[(length(s.1) - 5):length(s.1)]
      # s.max <- max(which(s.2 == max(s.2, na.rm=T)))
      # t.end <- round(t[s.max],0) + 1
      # t.trim <- t[which(t <= t.end)]
      # s.trim <- s.2[which(t <= t.end)]
      tmpdata <- data.frame(t = t.1, s = s.2)
    } else {
      s.2 <- as.numeric(s.1)
      # s.max <- max(which(s.2 == max(s.2, na.rm=T)))
      # t.end <- round(t[s.max],0) + 1
      # t.trim <- t[which(t <= t.end)]
      # s.trim <- s.2[which(t <= t.end)]
      tmpdata <- data.frame(t = t.1, s = s.2)
    }
    
    # Trimming Function
    if (trim == T){
      s.3 <- tmpdata$s
      t.3 <- tmpdata$t
      s.int <- 1
      s.end <- length(s.3)
      s.min <- min(which(s.3 == min(s.3, na.rm=T)))
      s.max <- max(which(s.3 == max(s.3, na.rm=T)))
      
      # define beginning
      if(s.3[s.int] - s.3[s.min] < delta){
        s.grw <- s.int
      } else {
        if(s.min > 0.5 * length(s.3)){
          plot(s.3 ~ t.3, main=samples[i], ylab="ABS", xlab="Time", pch=19)
          legend("topleft", legend = "OD appears to decreases", bty = "n")
          print(paste("Observed OD appears to decrease in sample ", 
                      samples[i], sep=""), quote=F)
          next
        } else {
          s.grw <- min(which(s.3 > (s.3[s.min] + delta)))
        }
      }
      
      # define end
      if((s.3[s.max] - s.3[s.end]) < delta){
        s.die <- s.end
      } else {
        s.die <- max(which(s.3 > s.3[s.max] - delta))
      }
      
      # trim beginning and end
      s.trim <- s.3[s.grw:s.die]
      t.trim <- t.3[s.grw:s.die]
      tmpdata <- data.frame(t = t.trim, s = s.trim)
      
    } else {
      s.3 <- tmpdata$s
      t.3 <- tmpdata$t
      s.trim <- s.3
      t.trim <- t.3
      tmpdata <- data.frame(t = t.trim, s = s.trim)
    }

    # Set Grid and Start Lists for Model
    # Modified Gompertz, dat~dnorm(mean=m.gomp(time2,c(b0,A,umax,L)),sd=exp(z))
    grids1<-list(A=c(0.8 * max(na.omit(tmpdata[,2])), max(na.omit(tmpdata[,2])), 
                     1.2 * max(na.omit(tmpdata[,2]))),
                 umax=c(0.05, 0.1, 1, 2),
                 L=c(1, 25, 50, 100),
                 z=c(-0.5, -0.05, -0.005))
    start1<-list(b0=intercept.guess,A=max(na.omit(tmpdata[,2])),umax=(2 * delta),L=0,z=-0.01)

    # Perform grid.mle2 Fits
    fit1<-grid.mle2(minuslogl=s~dnorm(mean=m.gomp(t,c(b0,A,umax,L)),
                    sd=exp(z)),grids=grids1,start=start1,data=tmpdata,
                    method="BFGS")
    print("finished fit", quote=F)

	  # isolate best of each class of model
    res.mat <- as.data.frame(fit1$res.mat)
    res.mat$sort <- seq(1:dim(fit1$res.mat)[1])
    pos <- res.mat[which(fit1$res.mat[,'A'] > 0), ]
    best <- pos$sort[which(pos$AIC == min(pos$AIC))]
    best.f1<-fit1$res.mod[[best]]
    #best.f1<-fit1$res.mod[[which(fit1$res.mat[,'AIC']==min(fit1$res.mat[,'AIC']))[1]]]

	  # generate plot of model fits
    # plot(s ~ t, main=samples[i], ylab="ABS", xlab="Time", pch=19, data=realdata)
    # curve(m.gomp(x,coef(best.f1)[1:4]),0,max(realdata$t),col='blue',lwd=2,add=T)
    # pdf(file=paste(temp.dir,"testplot",samples[i],".pdf",sep=""))
    par(mar = c(6, 6, 4, 2))
    plot(s ~ t, main=samples[i], 
         ylim = c(0, (coef(best.f1)[2] + coef(best.f1)[1]) * 1.1),
         ylab=expression(paste("Absorbance"[600])), 
         xlab=expression(paste("Time (hrs)")),
         pch=19, las = 1, cex.axis = 1, cex.lab = 1.5, data=realdata)
    curve(m.gomp(x + tmpdata$t[1],coef(best.f1)[1:4]),
          from = 0,
          to = max(tmpdata$t),
          col='blue',lwd=2,add=T)
    L = coef(best.f1)[4]
    A = coef(best.f1)[2]
    b0 = coef(best.f1)[1]
    umax = coef(best.f1)[3]
    #abline(v = L + tmpdata$t[1], lty = 2, col = "gray", lwd = 2)
    abline(v = L, lty = 2, col = "gray", lwd = 2)
    
    abline(h = A, lty = 2, col = "gray", lwd = 2)
    abline(h = A + b0, lty = 2, col = "blue", lwd = 2)
    #abline(-A + b0, umax, lty = 2, col = "red", lwd = 2)
    abline(-umax*(L- (b0/umax)), umax, lty = 2, col = "red", lwd = 2)
    abline(v = L - exp(1)*(b0/umax), lty = 2, col = "blue", lwd = 2)
    abline(h = 0, lty = 3, lwd = 2)
    legend("bottomright", legend = c("umax", "A and L", "Corrected A and L"),
           lty = 2, col = c("red", "gray", "blue"), bty = "n", cex = 0.8)

	# attempt profiling
    pf1<-profile(best.f1)
    print("finished pf1", quote=F)

    if(class(pf1)=="profile.mle2"){
    ci1<-suppressWarnings(confint(pf1)['umax',])
  	} else {
      ci1<-c(NA,NA)
		  names(ci1)<-c("2.5%","97.5%")
      }
    ciFI1<-suppressWarnings(confint.FI(best.f1)['umax',])

    # Save coefficients of model
    cfs<-coef(best.f1)
    results$Curve[i] <- samples[i]
    results$b0[i]<-round(cfs['b0'], 4)
    results$A[i]<-round(cfs['A'], 4)
    results$umax[i]<-round(cfs['umax'], 4)
    results$L[i]<-round(cfs['L'] + tmpdata$t[1], 4)
    results$z[i]<-round(cfs['z'], 4)
    results$K[i] <- round(cfs['A'] + cfs['b0'], 4)
    results$lag[i] <- round(cfs['L'] - exp(1)*(cfs['b0']/cfs['umax']), 4)
    results[i,c("umax.lw" , "umax.up" , "umax.lw.FI" , "umax.up.FI")] <- round(c(ci1,ciFI1), 4)
    write.table(results[i,], file=outfile, append=T, row.names=F, col.names=F, sep=",", quote=FALSE)

    # Print Operation Status
    print(paste("umax for ", samples[i], " = ", results$umax[i], sep = ""), quote = F)
    print(paste(round(((i)/(length(samples))*100),0),"% complete",
                sep = ""), quote=F)
  }
  results <- na.omit(results)
  # results1<-results
  # write.csv(results,"results.csv")
  dev.off()
  graphics.off()
  return(results)
  
}


