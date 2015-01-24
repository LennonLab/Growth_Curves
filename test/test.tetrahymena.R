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
source("../scripts/modified_Gomp.r")
# Create Directory For Output
dir.create("../output", showWarnings = FALSE)

# Run Example
growth.modGomp("../data/GrowthCurve_18hrs_140329_174411.txt", "tetra", start=147, t.cut=10)

growth.parameters <- read.csv(file="../output/tetra.txt", header=T)
design <- read.delim(file="../data/tetra.txt")

growth <- merge(growth.parameters, design, by.x = "Curve", by.y="Well", all.x=T)

se <- function(x, ...){sd(x, ...)/sqrt(length(x))}
growth.means <- tapply(growth$umax, growth$Treatment, mean, na.rm=TRUE)
growth.se <- tapply(growth$umax, growth$Treatment, se, na.rm=TRUE)

par(mar=c(10,5,1,1), oma=c(1,1,1,1)+0.1 )
plot1 <- barplot(growth.means, ylim=c(-0.5,0.2), las=2, ylab="Max Rate", cex.lab=2, yaxt="n", cex.names=1.25, lwd=2)
arrows(x0 = plot1, y0 = growth.means, y1 = growth.means - growth.se, angle = 90,
    length=0.05, lwd = 2, col="black")
arrows(x0 = plot1, y0 = growth.means, y1 = growth.means + growth.se, angle = 90,
    length=0.05, lwd = 2, col="black")
abline(h=0, lwd=2)
axis(side = 2, labels=T, lwd.ticks=2, las=1, lwd=2, cex.axis=1.25)
title(xlab = "Treatment", cex.lab = 2, line = 8.5)
