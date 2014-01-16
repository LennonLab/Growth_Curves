################################################################################
# Growth Curve Data Input From Synergy MX                                      #
#                                                                              #
#	Written by M. Muscarella                                                     #
#                                                                              #
#	Last update: 1/16/13                                                         #
#                                                                              #
################################################################################

read.synergy <- function(input = " "){
  data.in <- read.delim(input, skip=32, header=T)
  results.start <- which(data.in == "Results")
  data.out <- data.in[1:(results.start - 1),]
  colnames(data.out)[2] <- "Temp"
  return(data.out)
  }
