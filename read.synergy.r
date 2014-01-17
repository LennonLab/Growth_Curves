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

read.synergy <- function(input = " "){
  data.in <- read.delim(input, skip=32, header=T, as.is=T)
  results.start <- which(data.in == "Results")
  data.out <- data.in[1:(results.start - 1),]
  data.in.2 <- read.delim(input, skip=32, header=T)[1:(results.start - 1),]
  colnames(data.out)[2] <- "Temp"
  data.out$Time <- strftime(strptime(data.out$Time, format="%H:%M:%S"),"%M")
  for (i in 1:dim(data.out)[2]){
    if (is.numeric(data.out[,i]) == FALSE){data.out[,i] = as.numeric(data.out[,i])}}
  return(data.out)
  }                                                            
  
