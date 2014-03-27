################################################################################
#                                                                              #
#	Lennon Lab Growth Curve Analysis Example                                     #
#   Parameter Estimate Code                                                    #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: M. Muscarella                                                    #
#	  Last update: 3/27/14                                                       #
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
growth.modGomp("../test/GrowthCurve_test.txt", "test", start = 51)