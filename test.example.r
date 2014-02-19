################################################################################
#                                                                              #
#	Lennon Lab Growth Curve Analysis Example                                     #
#   Parameter Estimate Code                                                    #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: M. Muscarella                                                    #
#   Based on growthcurve_code.R Written by: M. Larsen (2013/07/18)             #
#                                                                              #
#	Last update: 2/19/14                                                         #
#                                                                              #
################################################################################

# Load Dependencies
source("growthcurve_models.r")
# Create Directory For Output
dir.create("./plot.out", showWarnings = FALSE)

# Run Example
growth.modGomp("GrowthCurve_Example.txt")