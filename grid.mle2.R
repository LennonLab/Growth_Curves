
####################################
#
#	Maximum likelihood model fitting using grid starts
#		- to calculate mle2 fits for a grid of initial parameter guesses
#		- initially designed for batch fitting bacterial culture growth curves
#		- designed as a wrapper for mle2() from Ben Bolker's bbmle package
#
#	Code by Colin T. Kremer
#
#	Last update: 4/21/12
#
#	To implement:
#		- improved passing of control parameters, etc to mle2 and optim
#		- clean up structure of function output
#
####################################

############# To fix:
#
# - output format is still a bit messy - learn more about nested lists.
# - passing control parameters to optim
#		- see functions match.call() and use of call in mle2 code.
#		- for now, handled in a crude way by forcing all models to use control 
#		  pars for parscale, hardwired into the code.
#
################

# Function requires:
#
# minuslogl = negative log likelihood function; passed to mle2
# grids = list of multiple starting guesses for parameters found in minuslogl 
# start = list of single starting guesses for parameters in minuslogl
# data = data frame containing the data used in the minuslogl function

# Function returns:

# A list containing 2 things:
#	1) res.mat = matrix of coefficient estimates and model AIC value
#	2) res.model = indexed list of mle2 objects, one for each fit, ordered as in res.mat

# for useage, see examples in grid.mle2_wexamples.R

grid.mle2<-function(minuslogl,grids,start,data,...){
	require(bbmle)

	if(length(grids)>=1){	# if grids have been supplied,

		# all combinations of grid variables and values
		grid.starts<-as.matrix(expand.grid(grids))
		ncombos<-dim(grid.starts)[[1]]
		
		# cycle through each combo
		res.mat<-matrix(NA,nrow=ncombos,ncol=I(length(start)+1))
		res.mod<-list()
		for(i in 1:dim(grid.starts)[[1]]){
		
			# some how need to match grid parameters to start lists.
			mod.start<-as.list(grid.starts[i,])	
			new.start<-start
			new.start[names(start) %in% names(mod.start)]<-mod.start

#			res.fit<-mle2(minuslogl=minuslogl,start=new.start,data=data,...)	

			pscale<-as.numeric(new.start)
			names(pscale)<-names(new.start)
			res.fit<-mle2(minuslogl=minuslogl,start=new.start,control=list(parscale=pscale),data=data,...)	

			res.mat[i,]<-c(coef(res.fit),AIC(res.fit))		
			res.mod[[i]]<-res.fit
		}
		colnames(res.mat)<-c(names(coef(res.fit)),"AIC")
	}else{	# otherwise, no grids; perform mle2 fit as usual
		res.mod<-list()
		pscale<-as.numeric(start)
		names(pscale)<-names(start)
		res.fit<-mle2(minuslogl=minuslogl,start=start,control=list(parscale=pscale),data=data,...)	

		res.mat<-c(coef(res.fit),AIC(res.fit))		
		res.mod[[1]]<-res.fit
		names(res.mat)<-c(names(coef(res.fit)),"AIC")
	}
	res<-list(res.mat=res.mat,res.mod=as.list(res.mod))
	res
}





