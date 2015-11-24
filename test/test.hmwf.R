################################################################################
#                                                                              #
#	Lennon Lab Growth Curve Analysis: HMWF Isolates 1-10                         #
#   Parameter Estimate Code                                                    #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: M. Muscarella                                                    #
#	  Last update: 4/03/14                                                       #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())
setwd('~/GitHub/growth-curves/test/')

# Load Dependencies
source("../bin/modified_Gomp.R")
# Create Directory For Output
dir.create("../output", showWarnings = FALSE)

# Run Example
growth.modGomp("../data/GrowthCurve_18hrs_140402_170936.txt", "HMWF", start=51)







growth.parameters <- read.csv(file="../output/HMWF.txt", header=T)
design <- read.delim(file="../data/hmwf_1-10.txt")

growth <- merge(growth.parameters, design, by.x = "Curve", by.y="Well", all.x=T)

se <- function(x, ...){sd(x, ...)/sqrt(length(x))}
growth.means <- tapply(growth$umax, growth$Isolate, mean)
growth.se <- tapply(growth$umax, growth$Isolate, se)

par(mar=c(8,5,1,1), oma=c(1,1,1,1)+0.1 )
plot1 <- barplot(growth.means, ylim=c(0,1.25), las=2, ylab="Max Growth Rate", cex.lab=2)
arrows(x0 = plot1, y0 = growth.means, y1 = growth.means - growth.se, angle = 90,
    length=0.05, lwd = 2, col="white")
arrows(x0 = plot1, y0 = growth.means, y1 = growth.means + growth.se, angle = 90,
    length=0.05, lwd = 2, col="black")
title(xlab = "HMWF Isolate", cex.lab = 2, line = 6.5)
