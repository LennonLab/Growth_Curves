################################################################################
#                                                                              #
#	Lennon Lab Growth Curve Analysis Example                                     #
#   Parameter Estimate Code                                                    #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: M. Muscarella                                                    #
#	  Last update: 2/19/14                                                       #
#                                                                              #
################################################################################

# Load Dependencies
source("../scripts/modified_Gomp.r")
# Create Directory For Output
dir.create("../output", showWarnings = FALSE)

# Run Example
growth.modGomp("../data/GrowthCurve_Example.txt", "test")