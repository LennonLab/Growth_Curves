################################################################################
#	Bacterial Growth Curve Models                                              #
#		- Modified Gompertz                                                    #
#		- new 'modified' gompertz                                              #
#         with 2 extra parameters capable of long term decay                   #
#		- confint.FI() function for calculating the Fisher Information 		   #
#  		  confidence  														   #
#      intervals; useful when profiling fails to converge/flat likelihood      #
#      surfaces                                                                #
#                                                                              #
#	Code by Colin T. Kremer                                                    #
#                                                                              #
#	Last update: 4/17/12                                                       #
#                                                                              #
################################################################################
# Tools for fitting modified modified gompertz                                 #
# Reference: ___________                                                       #
################################################################################

# Modified Gompertz Equation
m.gomp <- function(t, pars){
	b0 <-   pars[1]
	A <-    pars[2]
	umax <- pars[3]
	L <-    pars[4]
	b0+A*exp(-exp(umax*exp(1)*(L-t)/A+1))
  }

# New 'Modified' Gompertz Equation
# dd = controls the rate of decay in abundance experienced
# topt =  after this time point
# topt is constrained to be positive, and equal to or exceeding a positive offset of L
# A is also constrained to be positive
new.gomp <- function(t,pars){
	b0 <- pars[1]	     # intercept abundance (for curves that aren't standardized)
	A <- exp(pars[2])	 # formerly the carrying capacity (this interpretation fails in new model)
	umax <- pars[3]	     # maximum growth rate
	L <- pars[4]		 # lag before onset of exponential growth
	d <- pars[5]		 # controls rate of decreas after maximum population abundance
	topt <- pars[6]	     # controls onset of population decrease
	# constrain decay to set in after linear umax portion would reach A if unconstrained
	sig <- A/umax+L
    #denom <-1/(1 + ((t-L)^2/(exp(topt)+sig)^2)^d)
	denom <- 1/(1 + ((t)^2/(exp(topt)+sig)^2)^d)
	b0 + denom*(A*exp(-exp((umax*exp(1)*(L-t))/A+1)))
  }

# function to generate confidence intervals based on Fisher Information criteria
confint.FI <- function(model){
	cfs <- coef(model)
	ses <- sqrt(diag(vcov(model)))	# standard errors
	lw <- cfs-1.96*ses
	up <- cfs+1.96*ses
	res <- cbind(lw,up)
	dimnames(res) <- list(names(cfs),c("2.5 %","97.5 %"))
	res
  }

################################################################################
# examples                                                                     #
# confint.FI(fit.ng2)                                                          #
# coef(fit.ng2)                                                                #
################################################################################
