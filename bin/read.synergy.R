################################################################################
# Growth Curve Data Input From Synergy MX                                      #
#                                                                              #
#	Written by M. Muscarella                                                     #
#                                                                              #
#	Last update: 1/17/13                                                         #
#                                                                              #
# Features:                                                                    #
#   Reads in the datatable from the symergy mx machine                         #
#   Selects only the matrix of plate data                                      #
#   Converts time to minutes                                                   #
#   Checks the type of all values & changes to numeric if needed               #
#   Outputs data matric {Time, Temp, Well...                                   #
################################################################################

read.synergy <- function(input = " ", skip = ""){
  data.in <- read.delim(input, skip=skip, header=T, as.is=T)
  results.start <- which(data.in == "Results")
  data.out <- data.in[1:(results.start - 2),]
  colnames(data.out)[2] <- "Temp"
  data.out$Time <- as.character(data.out$Time)
  t.h <- as.numeric(lapply(strsplit(data.out$Time, "\\:"), "[", 1))
  t.m <- as.numeric(lapply(strsplit(data.out$Time, "\\:"), "[", 2))
  data.out$Time <- round(t.h + t.m/60, 2)
  for (i in 1:dim(data.out)[2]){
    if (is.numeric(data.out[,i]) == FALSE){data.out[,i] = as.numeric(data.out[,i])}}
  return(data.out)
  }
