#################################################################################
##
##   R package twinkle by Alexios Ghalanos Copyright (C) 2014.
##   This file is part of the R package twinkle.
##
##   The R package twinkle is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package twinkle is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################


.simregressors = function(model, xregsim, vregsim, ssim, N, n, m.sim)
{
	modelinc = model$modelinc
	mxn = modelinc[3]
	vxn = modelinc[39]
	sxn = modelinc[20]
	mar = max(c(modelinc[c(2,6,10,14,4,8,12,16,17)], ifelse(modelinc[49]==1, max(model$modeldata$ylags), 0)))
	m = max(c(modelinc[32:33], mar))
	# use coredata to extract pre-values since we do not work with xts in simulation
	if(mxn>0){
		if(is.null(xregsim)){
			xregsim = vector(mode = "list", length = m.sim)
			for(i in 1:m.sim) xregsim[[i]] = matrix(0, ncol = mxn, nrow = n)
		}
		if(!is.null(xregsim))
		{
			if(!is.list(xregsim)) stop("\nstarsim-->error: xregsim should be a list of length m.sim")
			if(length(xregsim) != m.sim){
				msd = vector(mode = "list", length = m.sim)
				for(i in 1:m.sim) msd[[i]] = as.matrix(xregsim[[1]])
				xregsim = msd
				warning("\nstarsim-->warning: length of xregsim list not equal to m.sim...\nreplicating first list element m.sim times.\n")
			}
			for(i in 1:m.sim){
				if(dim(as.matrix(xregsim[[i]]))[2] != mxn ) 
					stop(paste("\nstarsim-->error: xregsim ", i," has wrong no. of column", sep=""))
				if(dim(as.matrix(xregsim[[i]]))[1] != n )
					stop(paste("\nstarsim-->error: xregsim ", i," has wrong no. of rows", sep=""))
			}		
		}
		mexsimlist = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim){
			if(mar>0){
				premexdata = coredata(model$modeldata$mexdata[min(N,N-mar+1):N,,drop=FALSE])
				mexsimlist[[i]] = matrix(rbind(premexdata, as.matrix(xregsim[[i]])), ncol = mxn)
			} else{
				mexsimlist[[i]] = as.matrix(xregsim[[i]])
			}
		}
	} else{
		mexsimlist = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim) mexsimlist[[i]]=0
	}
	if(vxn>0){
		if(is.null(vregsim)){
			vregsim = vector(mode = "list", length = m.sim)
			for(i in 1:m.sim) vregsim[[i]] = matrix(0, ncol = vxn, nrow = n)
		}
		if(!is.null(vregsim))
		{
			if(!is.list(vregsim)) 
				stop("\nstarsim-->error: vregsim should be a list of length m.sim")
			if(length(vregsim) != m.sim){
				vsd = vector(mode = "list", length = m.sim)
				for(i in 1:m.sim) vsd[[i]] = as.matrix(vregsim[[1]])
				vregsim = vsd
				warning("\nstarsim-->warning: length of vregsim list not equal to m.sim...\nreplicating first list element m.sim times.\n")
			}
			for(i in 1:m.sim){
				if(dim(as.matrix(vregsim[[i]]))[2] != vxn ) 
					stop(paste("\nstarsim-->error: vregsim ", i," has wrong no. of column", sep=""))
				if(dim(as.matrix(vregsim[[i]]))[1] != n )
					stop(paste("\nstarsim-->error: vregsim ", i," has wrong no. of rows", sep=""))
			}
		}
		vexsimlist = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim){
			if(m>0){
				prevexdata = coredata(model$modeldata$vexdata[min(N,N-m+1):N,,drop=FALSE])
				vexsimlist[[i]] = matrix(rbind(prevexdata, as.matrix(vregsim[[i]])), ncol = vxn)
			} else{
				vexsimlist[[i]] = as.matrix(vregsim[[i]])
			}
		}
	} else{
		vexsimlist = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim) vexsimlist[[i]]=0
	}
	
	if(modelinc[49]==2 && sxn>0){
		if(is.null(ssim)){
			stop("\nstarsim-->error: ssim cannot be NULL in a STAR model which uses 's' in the probability dynamics!\n")
		}
		if(!is.null(ssim))
		{
			if(!is.list(ssim)) 
				stop("\nstarsim-->error: ssim should be a list of length m.sim")
			if(length(ssim) != m.sim){
				vsd = vector(mode = "list", length = m.sim)
				for(i in 1:m.sim) vsd[[i]] = as.matrix(ssim[[1]])
				ssim = vsd
				warning("\nstarsim-->warning: length of ssim list not equal to m.sim...\nreplicating first list element m.sim times.\n")
			}
			for(i in 1:m.sim){
				if(dim(as.matrix(ssim[[i]]))[2] != sxn ) 
					stop(paste("\nstarsim-->error: ssim ", i," has wrong no. of column", sep=""))
				if(dim(as.matrix(ssim[[i]]))[1] != n )
					stop(paste("\nstarsim-->error: ssim ", i," has wrong no. of rows", sep=""))
			}
		}
		ssimlist = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim){
			if(mar>0){
				presdata = coredata(model$modeldata$s[min(N,N-mar+1):N,,drop=FALSE])
				ssimlist[[i]] = matrix(rbind(presdata, as.matrix(ssim[[i]])), ncol = sxn)
			} else{
				ssimlist[[i]] = as.matrix(ssim[[i]])
			}
		}
	} else{
		ssimlist = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim) ssimlist[[i]]=0
	}
	return(list(mexsimlist = mexsimlist, vexsimlist = vexsimlist, ssimlist = ssimlist))
}

# use similar setup to aim
.simregressorspath = function(model, xregsim, vregsim, ssim, n, m.sim)
{
	modelinc = model$modelinc
	mxn = modelinc[3]
	vxn = modelinc[39]
	sxn = modelinc[20]
	# do we have external data (their indices match so we can loop and check)?
	N = NA
	if(modelinc[49]==2){
		N = nrow(model$modeldata$s)
	}
	if(modelinc[3]>0){
		N = nrow(model$modeldata$mexdata)
	}
	if(modelinc[39]>0){
		N = nrow(model$modeldata$vexdata)
	}
	# AR, MA (state), MA (linear) terms + ylags
	mar = max(c(modelinc[c(2,6,10,14,4,8,12,16,17)], ifelse(modelinc[49]==1, max(model$modeldata$ylags), 0)))
	m = max(c(modelinc[32:33], mar))
	# use coredata to extract pre-values since we do not work with xts in simulation
	if(mxn>0){
		if(is.null(xregsim)){
			xregsim = vector(mode = "list", length = m.sim)
			for(i in 1:m.sim){
				xregsim[[i]] = matrix(0, ncol = mxn, nrow = n)
			}
		}
		if(!is.null(xregsim))
		{
			if(!is.list(xregsim)) stop("\nstarpath-->error: xregsim should be a list of length m.sim")
			if(length(xregsim) != m.sim){
				msd = vector(mode = "list", length = m.sim)
				for(i in 1:m.sim) msd[[i]] = as.matrix(xregsim[[1]])
				xregsim = msd
				warning("\nstarpath-->warning: length of xregsim list not equal to m.sim...\nreplicating first list element m.sim times.\n")
			}
			for(i in 1:m.sim){
				if(dim(as.matrix(xregsim[[i]]))[2] != mxn ) 
					stop(paste("\nstarpath-->error: xregsim ", i," has wrong no. of column", sep=""))
				if(dim(as.matrix(xregsim[[i]]))[1] != n )
					stop(paste("\nstarpath-->error: xregsim ", i," has wrong no. of rows", sep=""))
			}		
		}
		mexsimlist = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim){
			if(mar>0){
				premexdata = coredata(model$modeldata$mexdata[min(N,N-mar+1):N,,drop=FALSE])
				mexsimlist[[i]] = matrix(rbind(premexdata, as.matrix(xregsim[[i]])), ncol = mxn)
			} else{
				mexsimlist[[i]] = as.matrix(xregsim[[i]])
			}
		}
	} else{
		mexsimlist = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim) mexsimlist[[i]]=0
	}
	if(vxn>0){
		if(is.null(vregsim)){
			vregsim = vector(mode = "list", length = m.sim)
			for(i in 1:m.sim) vregsim[[i]] = matrix(0, ncol = vxn, nrow = n)
		}
		if(!is.null(vregsim))
		{
			if(!is.list(vregsim)) 
				stop("\nstarpath-->error: vregsim should be a list of length m.sim")
			if(length(vregsim) != m.sim){
				vsd = vector(mode = "list", length = m.sim)
				for(i in 1:m.sim) vsd[[i]] = as.matrix(vregsim[[1]])
				vregsim = vsd
				warning("\nstarpath-->warning: length of vregsim list not equal to m.sim...\nreplicating first list element m.sim times.\n")
			}
			for(i in 1:m.sim){
				if(dim(as.matrix(vregsim[[i]]))[2] != vxn ) 
					stop(paste("\nstarpath-->error: vregsim ", i," has wrong no. of column", sep=""))
				if(dim(as.matrix(vregsim[[i]]))[1] != n )
					stop(paste("\nstarpath-->error: vregsim ", i," has wrong no. of rows", sep=""))
			}
		}
		vexsimlist = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim){
			if(m>0){
				prevexdata = coredata(model$modeldata$vexdata[min(N,N-m+1):N,,drop=FALSE])
				vexsimlist[[i]] = matrix(rbind(prevexdata, as.matrix(vregsim[[i]])), ncol = vxn)
			} else{
				vexsimlist[[i]] = as.matrix(vregsim[[i]])
			}
		}
	} else{
		vexsimlist = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim) vexsimlist[[i]]=0
	}
	
	if(modelinc[49]==2 && sxn>0){
		if(is.null(ssim)){
			stop("\nstarpath-->error: ssim cannot be NULL in a STAR model which uses 's' in the probability dynamics!\n")
		}
		if(!is.null(ssim))
		{
			if(!is.list(ssim)) 
				stop("\nstarpath-->error: ssim should be a list of length m.sim")
			if(length(ssim) != m.sim){
				vsd = vector(mode = "list", length = m.sim)
				for(i in 1:m.sim) vsd[[i]] = as.matrix(ssim[[1]])
				ssim = vsd
				warning("\nstarpath-->warning: length of ssim list not equal to m.sim...\nreplicating first list element m.sim times.\n")
			}
			for(i in 1:m.sim){
				if(dim(as.matrix(ssim[[i]]))[2] != sxn ) 
					stop(paste("\nstarpath-->error: ssim ", i," has wrong no. of column", sep=""))
				if(dim(as.matrix(ssim[[i]]))[1] != n )
					stop(paste("\nstarpath-->error: ssim ", i," has wrong no. of rows", sep=""))
			}
		}
		ssimlist = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim){
			if(mar>0){
				presdata = coredata(model$modeldata$s[min(N,N-mar+1):N,,drop=FALSE])
				ssimlist[[i]] = matrix(rbind(presdata, as.matrix(ssim[[i]])), ncol = sxn)
			} else{
				ssimlist[[i]] = as.matrix(ssim[[i]])
			}
		}
	} else{
		ssimlist = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim) ssimlist[[i]]=0
	}
	return(list(mexsimlist = mexsimlist, vexsimlist = vexsimlist, ssimlist = ssimlist))
}


.custzdist = function(custom.dist, distribution, mu = 0, sigma=1, skew = 0, shape=5, lambda=-0.5, m.sim, n)
{
	# ToDo: clean up (make use of rdist which is now vectorized)
	if(is.na(custom.dist$name) | is.na(custom.dist$distfit)[1]){
		z = matrix(rdist(distribution, n = m.sim*n, mu = mu, sigma = sigma, skew = skew, shape = shape, lambda = lambda), nrow = n, ncol = m.sim)
	}
	if(!is.na(custom.dist$name) && !is.na(custom.dist$distfit)[1]){
		
		if(is.matrix(custom.dist$distfit))
		{
			if(dim(custom.dist$distfit)[2]!=m.sim) stop("\ncolumn dimension of custom innovations matrix must be equal to m.sim\n")
			if(dim(custom.dist$distfit)[1]!=n) stop("\nrow dimension of custom innovations matrix must be equal to n.sim+n.start\n")
			z = custom.dist$distfit
		} else{
			if(!is.character(custom.dist$name)) stop("\ncustom distribution must be a character string\n")
			temp = paste("r", custom.dist$name, sep="")
			if(is.null(custom.dist$distfit)) stop("\ncustom distribution missing a distfit object\n")
			# if this is often used we might consider using apply to
			# use a new seed for each m.sim
			.rdist = eval(parse(text=paste(temp)))
			tmp = .rdist(n*m.sim, custom.dist$distfit)
			z = matrix(as.numeric(tmp), ncol = m.sim, nrow = n, byrow=TRUE)
		}
	}
	return(z)
}

dstar1sim = function(arglist)
{
	# unless the user passes his own fixed probabilities (i.e. weights in the 1-state case)
	# we set the value equal to the previous value.
	probs = rep(tail(arglist$probs,1), arglist$N)
	return(list(probs = probs))
}


dstar2sim = function(arglist)
{
	ipars = arglist$ipars
	model = arglist$model
	idx = model$pidx
	XL = arglist$XL
	modelinc = model$modelinc
	if(modelinc[21]>0){
		beta = ipars[idx["s1.beta",1],1]
	} else{
		beta = 0
	}
	gamma = ipars[idx["s1.gamma",1],1]
	alpha = ipars[idx["s1.alpha",1]:idx["s1.alpha",2],1]
	cnst = ipars[idx["s1.c",1],1]
	pmu = gamma*(-cnst + as.numeric(XL%*%alpha))
	if(modelinc[21]>0) pmu = .recfilter(as.double(pmu), as.double(beta), init = as.double(tail(arglist$pmu, 1)))
	probs = matrix(NA, ncol = 2, nrow = NROW(pmu))
	probs[,1] = 1/(1+exp(-pmu))
	probs[,2] = 1 - probs[,1]
	return(list(probs = probs, pmu = matrix(pmu, ncol=1)))
}


dstar3sim = function(arglist)
{
	ipars = arglist$ipars
	probs = arglist$probs
	model = arglist$model
	idx = model$pidx
	XL = arglist$XL
	modelinc = model$modelinc
	if(modelinc[21]>0){
		beta1 = ipars[idx["s1.beta",1],1]
	} else{
		beta1 = 0
	}
	if(modelinc[25]>0){
		beta2 = ipars[idx["s2.beta",1],1]
	} else{
		beta2 = 0
	}
	alpha1 = ipars[idx["s1.alpha",1]:idx["s1.alpha",2],1]
	alpha2 = ipars[idx["s2.alpha",1]:idx["s2.alpha",2],1]
	cnst1 = ipars[idx["s1.c",1],1]
	cnst2 = ipars[idx["s2.c",1],1]
	
	gamma1 = ipars[idx["s1.gamma",1],1]
	gamma2 = ipars[idx["s2.gamma",1],1]
	
	pmu1 = gamma1*(-cnst1 + as.numeric(XL%*%alpha1))
	pmu2 = gamma2*(-cnst2 + as.numeric(XL%*%alpha2))
	if(modelinc[21]>0) pmu1 = .recfilter(as.double(pmu1), as.double(beta1), init = as.double(tail(arglist$pmu[,1],1)))
	if(modelinc[25]>0) pmu2 = .recfilter(as.double(pmu2), as.double(beta2), init = as.double(tail(arglist$pmu[,2],1)))
	p1 = 1/(1+exp(pmu1))
	p2 = 1/(1+exp(pmu2))
	p12 = p1+p2
	p3 = 1/(1+p12)
	p1 = p1/(1+p12)
	p2 = p2/(1+p12)
	probs = matrix(NA, ncol = 3, nrow = NROW(pmu1))
	probs[,1] = p1
	probs[,2] = p2
	probs[,3] = p3
	return(list(probs = probs, pmu = cbind(pmu1, pmu2)))
}

dstar4sim = function(arglist)
{
	ipars = arglist$ipars
	probs = arglist$probs
	model = arglist$model
	idx = model$pidx
	XL = arglist$XL
	modelinc = model$modelinc
	if(modelinc[21]>0){
		beta1 = ipars[idx["s1.beta",1],1]
	} else{
		beta1 = 0
	}
	if(modelinc[25]>0){
		beta2 = ipars[idx["s2.beta",1],1]
	} else{
		beta2 = 0
	}
	if(modelinc[29]>0){
		beta3 = ipars[idx["s3.beta",1],1]
	} else{
		beta3 = 0
	}
	alpha1 = ipars[idx["s1.alpha",1]:idx["s1.alpha",2],1]
	alpha2 = ipars[idx["s2.alpha",1]:idx["s2.alpha",2],1]
	alpha3 = ipars[idx["s3.alpha",1]:idx["s3.alpha",2],1]
	cnst1 = ipars[idx["s1.c",1],1]
	cnst2 = ipars[idx["s2.c",1],1]
	cnst3 = ipars[idx["s3.c",1],1]
	
	gamma1 = ipars[idx["s1.gamma",1],1]
	gamma2 = ipars[idx["s2.gamma",1],1]
	gamma3 = ipars[idx["s3.gamma",1],1]
	
	pmu1 = gamma1*(-cnst1 + as.numeric(XL%*%alpha1))
	pmu2 = gamma2*(-cnst2 + as.numeric(XL%*%alpha2))
	pmu3 = gamma3*(-cnst3 + as.numeric(XL%*%alpha3))
	
	if(modelinc[21]>0) pmu1 = .recfilter(as.double(pmu1), as.double(beta1), init = as.double(tail(arglist$pmu[,1],1)))
	if(modelinc[25]>0) pmu2 = .recfilter(as.double(pmu2), as.double(beta2), init = as.double(tail(arglist$pmu[,2],1)))
	if(modelinc[29]>0) pmu3 = .recfilter(as.double(pmu3), as.double(beta3), init = as.double(tail(arglist$pmu[,3],1)))
	p1 = 1/(1+exp(pmu1))
	p2 = 1/(1+exp(pmu2))
	p3 = 1/(1+exp(pmu3))
	p123 = p1+p2+p3
	p4 = 1/(1+p123)
	p1 = p1/(1+p123)
	p2 = p2/(1+p123)
	p3 = p3/(1+p123)
	probs = matrix(NA, ncol = 4, nrow = NROW(pmu1))
	probs[,1] = p1
	probs[,2] = p2
	probs[,3] = p3
	probs[,4] = p4
	return(list(probs = probs, pmu = cbind(pmu1, pmu2, pmu3)))
}

dstar1path = function(arglist)
{
	# unless the user passes his own fixed probabilities (i.e. weights in the 1-state case)
	# we set the value equal to the previous value.
	probs = rep(tail(arglist$probs,1), arglist$N)
	return(list(probs = probs))
}


dstar2path = function(arglist)
{
	ipars = arglist$ipars
	model = arglist$model
	idx = model$pidx
	XL = arglist$XL
	modelinc = model$modelinc
	if(modelinc[20]>0){
		beta = ipars[idx["s1.beta",1],1]
	} else{
		beta = 0
	}
	alpha = ipars[idx["s1.alpha",1]:idx["s1.alpha",2],1]
	cnst = ipars[idx["s1.c",1],1]
	gamma = ipars[idx["s1.gamma",1],1]
	
	pmu = gamma*(-cnst + as.numeric(XL%*%alpha))
	if(modelinc[21]>0){
		init1 = ( gamma*(-cnst+sum(alpha * colMeans(XL))) ) /(1-beta)
		pmu = .recfilter(as.double(pmu), as.double(beta), init = as.double(init1))
	}
	probs = matrix(NA, ncol = 2, nrow = NROW(pmu))
	probs[,1] = 1/(1+exp(-pmu))
	probs[,2] = 1 - probs[,1]
	return(list(probs = probs, pmu = matrix(pmu, ncol=1)))
}


dstar3path = function(arglist)
{
	ipars = arglist$ipars
	probs = arglist$probs
	model = arglist$model
	idx = model$pidx
	XL = arglist$XL
	modelinc = model$modelinc
	if(modelinc[21]>0){
		beta1 = ipars[idx["s1.beta",1],1]
	} else{
		beta1 = 0
	}
	if(modelinc[25]>0){
		beta2 = ipars[idx["s2.beta",1],1]
	} else{
		beta2 = 0
	}
	alpha1 = ipars[idx["s1.alpha",1]:idx["s1.alpha",2],1]
	alpha2 = ipars[idx["s2.alpha",1]:idx["s2.alpha",2],1]
	cnst1 = ipars[idx["s1.c",1],1]
	cnst2 = ipars[idx["s2.c",1],1]
	gamma1 = ipars[idx["s1.gamma",1],1]
	gamma2 = ipars[idx["s2.gamma",1],1]
	
	pmu1 = gamma1*(-cnst1 + as.numeric(XL%*%alpha1))
	pmu2 = gamma2*(-cnst2 + as.numeric(XL%*%alpha2))
	if(modelinc[21]>0){
		init1 = (gamma1*(-cnst1+sum(alpha1 * colMeans(XL)))) /(1-beta1)
		pmu1 = .recfilter(as.double(pmu1), as.double(beta1), init = as.double(init1))
	}
	if(modelinc[25]>0){
		init2 = (gamma2*(-cnst2+sum(alpha2 * colMeans(XL)))) /(1-beta2)
		pmu2 = .recfilter(as.double(pmu2), as.double(beta2), init = as.double(init2))
	}
	p1 = 1/(1+exp(pmu1))
	p2 = 1/(1+exp(pmu2))
	p12 = p1+p2
	p3 = 1/(1+p12)
	p1 = p1/(1+p12)
	p2 = p2/(1+p12)
	probs = matrix(NA, ncol = 3, nrow = NROW(pmu1))
	probs[,1] = p1
	probs[,2] = p2
	probs[,3] = p3
	return(list(probs = probs, pmu = cbind(pmu1, pmu2)))
}

dstar4path = function(arglist)
{
	ipars = arglist$ipars
	probs = arglist$probs
	model = arglist$model
	idx = model$pidx
	XL = arglist$XL
	modelinc = model$modelinc
	if(modelinc[21]>0){
		beta1 = ipars[idx["s1.beta",1],1]
	} else{
		beta1 = 0
	}
	if(modelinc[25]>0){
		beta2 = ipars[idx["s2.beta",1],1]
	} else{
		beta2 = 0
	}
	if(modelinc[29]>0){
		beta3 = ipars[idx["s3.beta",1],1]
	} else{
		beta3 = 0
	}
	alpha1 = ipars[idx["s1.alpha",1]:idx["s1.alpha",2],1]
	alpha2 = ipars[idx["s2.alpha",1]:idx["s2.alpha",2],1]
	alpha3 = ipars[idx["s3.alpha",1]:idx["s3.alpha",2],1]
	cnst1 = ipars[idx["s1.c",1],1]
	cnst2 = ipars[idx["s2.c",1],1]
	cnst3 = ipars[idx["s3.c",1],1]
	gamma1 = ipars[idx["s1.gamma",1],1]
	gamma2 = ipars[idx["s2.gamma",1],1]
	gamma3 = ipars[idx["s3.gamma",1],1]
	pmu1 = gamma1*(-cnst1 + as.numeric(XL%*%alpha1))
	pmu2 = gamma2*(-cnst2 + as.numeric(XL%*%alpha2))
	pmu3 = gamma3*(-cnst3 + as.numeric(XL%*%alpha3))
	
	if(modelinc[21]>0){
		init1 = (gamma1*(-cnst1+sum(alpha1 * colMeans(XL)))) /(1-beta1)
		pmu1 = .recfilter(as.double(pmu1), as.double(beta1), init = as.double(init1))
	}
	if(modelinc[25]>0){
		init2 = (gamma2*(-cnst2+sum(alpha2 * colMeans(XL)))) /(1-beta2)
		pmu2 = .recfilter(as.double(pmu2), as.double(beta2), init = as.double(init2))
	}
	if(modelinc[29]>0){
		init3 = (gamma3*(-cnst3+sum(alpha3 * colMeans(XL)))) /(1-beta3)
		pmu3 = .recfilter(as.double(pmu3), as.double(beta3), init = as.double(init3))
	}
	p1 = 1/(1+exp(pmu1))
	p2 = 1/(1+exp(pmu2))
	p3 = 1/(1+exp(pmu3))
	p123 = p1+p2+p3
	p4 = 1/(1+p123)
	p1 = p1/(1+p123)
	p2 = p2/(1+p123)
	p3 = p3/(1+p123)
	probs = matrix(NA, ncol = 4, nrow = NROW(pmu1))
	probs[,1] = p1
	probs[,2] = p2
	probs[,3] = p3
	probs[,4] = p4
	return(list(probs = probs, pmu = cbind(pmu1, pmu2, pmu3)))
}

# simulate 1-ahead under a range of values for probability-dynamics and external regressors
# .simulate1ahead = function()
# allow VAR or ARMA to be used for each series to create a forecast with error bands.