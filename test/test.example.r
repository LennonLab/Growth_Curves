################################################################################
#                                                                              #
#	Lennon Lab Growth Curve Analysis Example                                     #
#   Parameter Estimate Code                                                    #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: M. Muscarella                                                    #
#	  Last update: 11/24/2015 by M. Muscarella & V. Kuo                          #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())
setwd("~/GitHub/Growth_Curves/test")

# Load Dependencies
source("../bin/modified_Gomp.R")
# Create Directory For Output
dir.create("../output", showWarnings = FALSE)

# Run Example (Synergy)
growth.modGomp("../data/GrowthCurve_Example.txt", "test", skip=31)
growth.modGomp("../data/GrowthCurve_Example2.txt", "test", skip=48)
growth.modGomp("../data/RPF.txt", "RPF", skip = 39)

# Run Example (CSV file)
growth.modGomp("../data/ControlRpf.csv", "csv_test", synergy=F, temp=F)
