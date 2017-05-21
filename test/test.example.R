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
setwd("C:/Users/Venus/Github/Growth_Curves/test")

# Load Dependencies
source("../bin/modified_Gomp.R")

# Run Example (CSV file)
# CSV file must have the following format: Header followed by whole or decimal numbers 
# Time  A1  A2  A3  ...
# 30    0.2 0.1 0.4 ...
growth.modGomp("../data/Brentgrowthcurve.csv", "csv_test", synergy=F, temp=F)

growth.modGomp("../data/20170305_GeneralityRpf.csv", "csv_test", synergy=F, temp=F)

# Run Example (Synergy)
growth.modGomp("../data/GrowthCurve_Example.txt", "test", skip=31)
growth.modGomp("../data/GrowthCurve_Example2.txt", "test", skip=48)
growth.modGomp("../data/RPF.txt", "RPF", skip = 39)
growth.modGomp("../test/125_1_6_GrowthCurves.txt", "BehriTest", skip = 38)

# Run Example (CSV file)
growth.modGomp("../data/ControlRpf.csv", "csv_test", synergy=F, temp=F)
growth.modGomp("../data/Pseudo.csv", "csv_test", synergy=F, temp=F, smooth=F)
