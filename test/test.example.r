################################################################################
#                                                                              #
#	Lennon Lab Growth Curve Analysis Example                                     #
#   Parameter Estimate Code                                                    #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: M. Muscarella                                                    #
#	  Last update: 1/24/15                                                       #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())
setwd('~/GitHub/Growth_Curves/test/')

# Load Dependencies
source("../bin/modified_Gomp.r")
# Create Directory For Output
dir.create("../output", showWarnings = FALSE)

# Run Example
growth.modGomp("../data/GrowthCurve_Example.txt", "test", skip=31)
growth.modGomp("../data/GrowthCurve_Example2.txt", "test", skip=48)
growth.modGomp("../data/RPF.txt", "RPF", skip = 39)
