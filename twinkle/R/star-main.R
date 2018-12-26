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

.starfit = function(spec, data, out.sample = 0, solver = "optim", solver.control = list(), 
		fit.control = list(stationarity = 0, fixed.se = 0, rec.init = 'all'), cluster = NULL, 
		n = 25, ...)
{
	if(solver=="strategy"){
		if(is.null(solver.control$trace)) trace=0 else trace = solver.control$trace
		fit = .starfit.strategy(spec = spec, data = data, out.sample = out.sample, fit.control = fit.control, n = n, trace = trace)
		
	} else{
		fit = .starfit.main(spec = spec, data = data, out.sample = out.sample, solver = solver, 
				solver.control = solver.control, fit.control = fit.control, cluster = cluster)
	}
	return(fit)
}

.starfit.main = function(spec, data, out.sample = 0, solver = "optim", solver.control = list(), 
		fit.control = list(stationarity = 0, fixed.se = 0, rec.init = 'all'), cluster = NULL, 
		...)
{
	tic = Sys.time()
	model = spec@model
	modelinc = model$modelinc
	if(is.null(fit.control$stationarity)) fit.control$stationarity = FALSE
	if(is.null(solver.control$trace)) trace = 0 else trace = solver.control$trace
	if(is.null(fit.control$fixed.se)) fit.control$fixed.se = FALSE
	if(is.null(fit.control$rec.init)) fit.control$rec.init = 'all'
	mm = match(names(fit.control), c("stationarity", "fixed.se", "rec.init"))
	if(any(is.na(mm))){
		idx = which(is.na(mm))
		enx = NULL
		for(i in 1:length(idx)) enx = c(enx, names(fit.control)[idx[i]])
		warning(paste(c("unidentified option(s) in fit.control:\n", enx), sep="", collapse=" "), call. = FALSE, domain = NULL)
	}
	n.start = round(out.sample,0)
	# extract data and index
	xdata = .extractdata(data)
	if(!is.numeric(out.sample)) stop("\nstarfit-->error: out.sample must be numeric\n")
	if(as.numeric(out.sample)<0) stop("\nstarfit-->error: out.sample must be positive\n")
	n = length(xdata$data)
	#---------------------------------------------------------
	# check x (state regressors)
	#---------------------------------------------------------
	if(modelinc[47]==0){
		if(modelinc[49]==2){
			chk = all(index(data)==index(model$modeldata$s))
			if(!is.logical(chk) | chk == FALSE){
				print(paste("\n",chk,sep=""))
				stop("\nstarfit-->error: data and s indices do not match\n")
			}
			XL = model$modeldata$s[1:(n-n.start), , drop = FALSE]
		} else{
			# we apply yfun here if it is needed for efficiency:
			if(modelinc[48]==1){
				ytmp = xts(model$modeldata$fun(data[1:(n-n.start),1]), index(data[1:(n-n.start),1]))
			} else{
				ytmp = data[1:(n-n.start),1]
			}
			XL = NULL
			for(i in 1:length(model$modeldata$ylags)){
				if(i==1) XL = lag(ytmp, model$modeldata$ylags[i]) else XL = cbind(XL, lag(ytmp, model$modeldata$ylags[i]))
			}
			XL[is.na(XL)]=0
		}
		# convert to matrix
		XL = coredata(XL)
	} else{
		XL = NULL
	}
	# XL is the lagged regressor dataset for the state dynamics
	#---------------------------------------------------------
	
	#---------------------------------------------------------
	# check fixed.probs (in case a set of fixed state probabilities
	# have been provided)
	#---------------------------------------------------------
	if(modelinc[47]==1){
		chk = all(index(data)==index(model$fixed.prob))
		if(!is.logical(chk) | chk == FALSE){
			print(paste("\n",chk,sep=""))
			stop("\nstarfit-->error: data and fixed.probs indices do not match\n")
		}
		probs = coredata(model$fixed.prob[1:(n-n.start), , drop = FALSE])
	} else{
		probs = matrix(0, ncol = modelinc[46], nrow = (n-n.start))
	}
	#---------------------------------------------------------
	
	#---------------------------------------------------------
	# check ARX external regressors (X)
	#---------------------------------------------------------
	if(modelinc[3] > 0){
		chk = all(index(data)==index(model$modeldata$mexdata))
		if(!is.logical(chk) | chk == FALSE){
			print(paste("\n",chk,sep=""))
			stop("\nstarfit-->error: data and xreg indices do not match\n")
		}
		# lag mexdata according to xlags
		mexdata = coredata(model$modeldata$mexdata[1:(n-n.start), , drop = FALSE])
	} else{
		mexdata = NULL
	}
	
	if(modelinc[39] > 0){
		chk = all(index(data)==index(model$modeldata$vexdata))
		if(!is.logical(chk) | chk == FALSE){
			print(paste("\n",chk,sep=""))
			stop("\nstarfit-->error: data and vreg indices do not match\n")
		}
		# lag mexdata according to xlags
		vexdata = coredata(model$modeldata$vexdata[1:(n-n.start), , drop = FALSE])
	} else{
		vexdata = NULL
	}
	#---------------------------------------------------------
	# index the data to start at point n-out.sample
	data = xdata$data[1:(n-n.start)]
	index = xdata$index[1:(n-n.start)]
	origdata = xdata$data
	origindex = xdata$index
	period = xdata$period
	#---------------------------------------------------------
	# create a temporary environment to store values (deleted at end of function)
	# and the list of parameters/data passed to the likelihood estimation functions
	starenv = new.env(hash = TRUE)
	
	arglist = list()
	arglist$starenv <- starenv
	# parallel mode flag (related to use of starenv environment)
	arglist$pmode = 0
	pidx = model$pidx
	arglist$index = index
	arglist$trace = trace
	arglist$fit.control = fit.control
	model$modeldata$T = T = length(as.numeric(data))
	recinit = .checkrec(fit.control$rec.init, T)
	arglist$recinit = recinit
	dist = model$modeldesc$distribution
	arglist$data = data
	arglist$XL = XL
	# in the ugarchfilter, N acts as n.old to avoid lookahead bias in filtering
	# for initialization values
	arglist$N = T
	arglist$mexdata = mexdata
	arglist$vexdata = vexdata
	arglist$probs = probs
	arglist$model = model
	ipars = model$pars
	# Optimization Starting Parameters Vector & Bounds
	ipars = .starstart(ipars, arglist = arglist)
	arglist$ipars = ipars
	# we now split out any fixed parameters
	estidx = as.logical( ipars[,4] )
	arglist$estidx = estidx	
	npars = sum(estidx)
	if(any(ipars[,2]==1)){
		if(npars == 0){
			if(fit.control$fixed.se==0) {
				# if all parameters are fixed an no standard erros are to
				# be calculated then we return a ugarchfilter object
				warning("\nstarfit-->warning: all parameters fixed...returning starfilter object instead\n")
				return(starfilter(data = xts(origdata, origindex), spec = spec, out.sample = out.sample))
			} else{
				# if all parameters are fixed but we require standard errors, we
				# skip the solver
				use.solver = 0
				ipars[ipars[,2]==1, 4] = 1
				ipars[ipars[,2]==1, 2] = 0
				arglist$ipars = ipars
				estidx = as.logical( ipars[,4] )
				arglist$estidx = estidx	
			}
		} else{
			# with some parameters fixed we extract them (to be rejoined at end)
			# so that they do not enter the solver
			use.solver = 1
		}
	} else{
		use.solver = 1
	}
	# start counter
	assign("star_llh", 1, envir = starenv)
	arglist$fit.control = fit.control
	# need to control for use of parallel environment which cannot handle the
	# calls between the locally stored environment
	if(modelinc[50]>0){
		if(modelinc[50]==4){
			# There is no 1-state mixture so just use the 2-state for the switching
			llfun = switch(modelinc[46],
					.stars2dLLH2mix,
					.stars2dLLH2mix,
					.stars3dLLH3mix,
					.stars4dLLH4mix)
		} else{
			llfun = switch(modelinc[46],
					.stars1dLLH,
					.stars2dLLH,
					.stars3dLLH,
					.stars4dLLH)
		}
	} else{
		llfun = switch(modelinc[46],
				.stars1sLLH,
				.stars2sLLH,
				.stars3sLLH,
				.stars4sLLH)
	}
	if(use.solver){
		parscale = rep(1, length = npars)
		names(parscale) = rownames(ipars[estidx,])
		arglist$returnType = "llh"
		# switch depending on states:
		solution = .starsolver(solver, pars = ipars[estidx, 1], fun = llfun, 
				Ifn = NULL, ILB = NULL, IUB = NULL, gr = NULL, hessian = NULL, 
				parscale = parscale, control = solver.control, arglist = arglist,
				cluster = cluster)
		sol = solution$sol
		hess = solution$hess
		timer = Sys.time()-tic		
		if(!is.null(sol$pars)){
			ipars[estidx, 1] = sol$pars
			if(modelinc[50]>0 && modelinc[31]==0){
				# call it once more to get omega
				ipars[pidx["omega",1], 1] = get("omega", starenv)
			}
		} else{
			ipars[estidx, 1] = NA
		}
		arglist$ipars = ipars
		convergence = sol$convergence
	} else{
		hess = NULL
		timer = Sys.time()-tic
		convergence = 0
		sol = list()
		sol$message = "all parameters fixed"
	}
	fit = list()
	# check convergence else write message/return
	ipars2 = ipars
	
	if(convergence == 0){
		if(sum(ipars[,2]) > 0 && fit.control$fixed.se == 1){
			ipars[ipars[,2]==1, 4] = 1
			ipars[ipars[,2]==1, 2] = 0
			arglist$ipars = ipars
			estidx = as.logical( ipars[,4] )
			arglist$estidx = estidx	
		}
		fit = .starmodelpost(f = llfun, T = T, timer = timer, convergence = convergence, message = sol$message, hess, 
				arglist = arglist)
		model$modeldata$data = origdata
		model$modeldata$index = origindex
		model$modeldata$period = period
		model$pars = ipars
		model$pars[, 1] = fit$ipars[,1]
		fit$ipars[, 4] = ipars2[, 4]
		fit$ipars[, 2] = ipars2[, 2]
	} else{
		fit$message = sol$message
		fit$convergence = 1
		model$modeldata$data = origdata
		model$modeldata$index = origindex
		model$modeldata$period = period
	}
	# make model list to return some usefule information which
	# will be called by other functions (show, plot, sim etc)
	model$n.start = n.start
	fit$solver = solution
	ans = new("STARfit",
			fit = fit,
			model = model)
	return(ans)
}



.starfitx = function(spec, data, out.sample = 0, solver.control = list(),
		fit.control = list(stationarity = 0, fixed.se = 0, rec.init = 'all'), n = 10, 
		only.start = FALSE, usepars=FALSE, transform.result = TRUE, ...)
{
	tic = Sys.time()
	model = spec@model
	modelinc = model$modelinc
	if(is.null(fit.control$stationarity)) fit.control$stationarity = FALSE
	if(is.null(solver.control$trace)) trace = 0 else trace = solver.control$trace
	if(is.null(fit.control$fixed.se)) fit.control$fixed.se = FALSE
	if(is.null(fit.control$rec.init)) fit.control$rec.init = 'all'
	mm = match(names(fit.control), c("stationarity", "fixed.se", "rec.init"))
	if(any(is.na(mm))){
		idx = which(is.na(mm))
		enx = NULL
		for(i in 1:length(idx)) enx = c(enx, names(fit.control)[idx[i]])
		warning(paste(c("unidentified option(s) in fit.control:\n", enx), sep="", collapse=" "), call. = FALSE, domain = NULL)
	}
	n.start = round(out.sample,0)
	# extract data and index
	xdata = .extractdata(data)
	if(!is.numeric(out.sample)) stop("\nstarfit-->error: out.sample must be numeric\n")
	if(as.numeric(out.sample)<0) stop("\nstarfit-->error: out.sample must be positive\n")
	n = length(xdata$data)
	#---------------------------------------------------------
	# check x (state regressors)
	#---------------------------------------------------------
	if(modelinc[47]==0){
		if(modelinc[49]==2){
			chk = all(index(data)==index(model$modeldata$s))
			if(!is.logical(chk) | chk == FALSE){
				print(paste("\n",chk,sep=""))
				stop("\nstarfit-->error: data and x indices do not match\n")
			}
			XL = model$modeldata$s[1:(n-n.start), , drop = FALSE]
		} else{
			# we apply yfun here if it is needed for efficiency:
			if(modelinc[48]==1){
				ytmp = xts(model$modeldata$fun(data[1:(n-n.start),1]), index(data[1:(n-n.start),1]))
			} else{
				ytmp = data[1:(n-n.start),1]
			}
			XL = NULL
			for(i in 1:length(model$modeldata$ylags)){
				if(i==1) XL = lag(ytmp, model$modeldata$ylags[i]) else XL = cbind(XL, lag(ytmp, model$modeldata$ylags[i]))
			}
			XL[is.na(XL)]=0
		}
		XL = coredata(XL)
	} else{
		XL = NULL
	}
	# XL is the lagged regressor dataset for the state dynamics
	#---------------------------------------------------------
	
	#---------------------------------------------------------
	# check fixed.probs (in case a set of fixed state probabilities
	# have been provided)
	#---------------------------------------------------------
	if(modelinc[47]==1){
		chk = all(index(data)==index(model$fixed.prob))
		if(!is.logical(chk) | chk == FALSE){
			print(paste("\n",chk,sep=""))
			stop("\nstarfit-->error: data and fixed.probs indices do not match\n")
		}
		probs = coredata(model$fixed.prob[1:(n-n.start), , drop = FALSE])
	} else{
		probs = matrix(0, ncol = modelinc[46], nrow = (n-n.start))
	}
	#---------------------------------------------------------
	
	#---------------------------------------------------------
	# check ARX external regressors (X)
	#---------------------------------------------------------
	if(modelinc[3] > 0){
		chk = all(index(data)==index(model$modeldata$mexdata))
		if(!is.logical(chk) | chk == FALSE){
			print(paste("\n",chk,sep=""))
			stop("\nstarfit-->error: data and xreg indices do not match\n")
		}
		# lag mexdata according to xlags
		mexdata = coredata(model$modeldata$mexdata[1:(n-n.start), , drop = FALSE])
	} else{
		mexdata = NULL
	}
	
	if(modelinc[39] > 0){
		chk = all(index(data)==index(model$modeldata$vexdata))
		if(!is.logical(chk) | chk == FALSE){
			print(paste("\n",chk,sep=""))
			stop("\nstarfit-->error: data and vreg indices do not match\n")
		}
		# lag mexdata according to xlags
		vexdata = coredata(model$modeldata$vexdata[1:(n-n.start), , drop = FALSE])
	} else{
		vexdata = NULL
	}
	#---------------------------------------------------------
	# index the data to start at point n-out.sample
	data = xdata$data[1:(n-n.start)]
	index = xdata$index[1:(n-n.start)]
	origdata = xdata$data
	origindex = xdata$index
	period = xdata$period
	#---------------------------------------------------------
	# create a temporary environment to store values (deleted at end of function)
	# and the list of parameters/data passed to the likelihood estimation functions
	starenv = new.env(hash = TRUE)
	
	arglist = list()
	arglist$starenv <- starenv
	# parallel mode flag (related to use of starenv environment)
	arglist$pmode = 0
	pidx = model$pidx
	arglist$index = index
	arglist$trace = trace
	arglist$fit.control = fit.control
	model$modeldata$T = T = length(as.numeric(data))
	recinit = .checkrec(fit.control$rec.init, T)
	arglist$recinit = recinit
	dist = model$modeldesc$distribution
	arglist$data = data
	arglist$XL = XL
	# in the ugarchfilter, N acts as n.old to avoid lookahead bias in filtering
	# for initialization values
	arglist$N = T
	arglist$mexdata = mexdata
	arglist$vexdata = vexdata
	arglist$probs = probs
	arglist$model = model
	ipars = model$pars
	# Optimization Starting Parameters Vector & Bounds
	ipars = .starstart(ipars, arglist = arglist)
	arglist$ipars = ipars
	# we now split out any fixed parameters
	estidx = as.logical( ipars[,4] )
	arglist$estidx = estidx	
	npars = sum(estidx)
	if(any(ipars[,2]==1)){
		if(npars == 0){
			if(fit.control$fixed.se==0){
				# if all parameters are fixed an no standard erros are to
				# be calculated then we return a ugarchfilter object
				warning("\nstarfit-->warning: all parameters fixed...returning starfilter object instead\n")
				return(starfilter(data = xts(origdata, origindex), spec = spec, out.sample = out.sample))
			} else{
				# if all parameters are fixed but we require standard errors, we
				# skip the solver
				use.solver = 0
				ipars[ipars[,2]==1, 4] = 1
				ipars[ipars[,2]==1, 2] = 0
				arglist$ipars = ipars
				estidx = as.logical( ipars[,4] )
				arglist$estidx = estidx	
			}
		} else{
			# with some parameters fixed we extract them (to be rejoined at end)
			# so that they do not enter the solver
			use.solver = 1
		}
	} else{
		use.solver = 1
	}
	# start counter
	assign("star_llh", 1, envir = starenv)
	arglist$fit.control = fit.control
	arglist$transform=TRUE
	# need to control for use of parallel environment which cannot handle the
	# calls between the locally stored environment
	if(modelinc[50]>0){
		if(modelinc[50]==4){
			llfun = switch(modelinc[46],
					.stars2dLLH2mix,
					.stars2dLLH2mix,
					.stars3dLLH3mix,
					.stars4dLLH4mix)
		} else{
			llfun = switch(modelinc[46],
					.stars1dLLH,
					.stars2dLLH,
					.stars3dLLH,
					.stars4dLLH)
		}
	} else{
		llfun = switch(modelinc[46],
				.stars1sLLH,
				.stars2sLLH,
				.stars3sLLH,
				.stars4sLLH)
	}
	parscale = rep(1, length = npars)
	names(parscale) = rownames(ipars[estidx,])
	arglist$returnType = "llh"
	if(only.start){
		solver.control$searchtype="likelihood"
		# solvertype is default  (unconstrained)
		solution = .starparsearch(llfun, arglist = arglist, solver = "optim", solver.control = solver.control, cluster = NULL,
				solvertype="constrained")
		pars = solution[1,]
		names(pars) = rownames(ipars[estidx,])
		arglist$transform=FALSE
		return(list(coef = pars, lik = llfun(pars,arglist)))
	} else{
		# switch depending on states:
		solution = .bfgssolver(pars = ipars[estidx, 1], fun = llfun, gr = NULL, control = solver.control,  arglist = arglist, cluster = NULL, 
				usepars = usepars, transform.result = transform.result)
		if(solution$sol$convergence==0 && !is.null(solution$sol$pars)){
			pars = solution$sol$pars
			names(pars) = rownames(ipars[estidx,])
		} else{
			pars = rep(NA, length(estidx))
			solution$sol$convergence=1
		}
		return(list(coef = pars, lik = solution$sol$value, convergence = solution$sol$convergence))
	}
}



.starfit.strategy = function(spec, data, out.sample = 0, fit.control = list(stationarity = 0, fixed.se = 0, rec.init = 'all'), n = 25, trace = 1, ...)
{
	tic = Sys.time()
	model = spec@model
	modelinc = model$modelinc
	solver.control = list(trace=0, method="BFGS", reltol=1e-6, maxit=1000)
	idx = spec@model$pos.matrix
	pest = spec@model$pars[,4]
	ptyp = spec@model$pars[,8]
	plin = 1:length(which(pest==1 & ptyp==1))
	nn = max(plin, 0)
	if(length(which(pest==1 & ptyp==2))>0){
		pnlin = (nn+1):(nn+length(which(pest==1 & ptyp==2)))
		nn = max(pnlin, 0)
	}
	if(length(which(pest==1 & ptyp==3))>0) pnoth = (nn+1):(nn+length(which(pest==1 & ptyp==3)))
	if(length(plin)>0) idx_lin = plin else idx_lin = NULL
	if(length(pnlin)==0) stop("\nstarfit-->error: strategy cannot be used when (most probably) state probabilities are fixed...")
	idx_state = pnlin
	idx_other = pnoth
	sol = .starfitx(spec, data, out.sample = out.sample, solver.control = solver.control,
			fit.control = fit.control, only.start=TRUE)
	# this returns a set of (transformed) parameters
	if(trace) cat("\niter: 1  (initial): ", sol$lik)
	solver.control = list(trace=0, method="BFGS", reltol=1e-12, maxit=1000)
	minset = list(cf = sol$coef, lik = sol$lik)
	for(i in 1:n)
	{
		specx = spec
		setfixed(specx)<-as.list(minset$cf[c(idx_lin, idx_other)])
		solver.control = list(trace=0, method="BFGS", reltol=1e-12, maxit=10000,parsearch=TRUE)
		sol1 = .starfitx(specx, data, out.sample = out.sample, solver.control = solver.control,
				fit.control = fit.control, transform.result=TRUE, usepars=FALSE)
		
		if(sol1$convergence==0 && sol1$lik<minset$lik){
			if(trace) cat("\niter:",i+1,"   (state): ", sol1$lik)
			minset$lik = sol1$lik
			minset$cf[c(idx_state)]  = sol1$coef
		} else{
			if(trace) cat("\niter:",i+1,"   (state): no change")
		}
		specx = spec
		setfixed(specx)<-as.list(minset$cf[idx_state])
		solver.control = list(trace=0, method="BFGS", reltol=1e-12, maxit=4000, parsearch=TRUE)
		sol2 = .starfitx(specx, data, out.sample = out.sample, solver.control = solver.control,
				fit.control = fit.control, transform.result=TRUE)
		if(sol2$convergence==0 && sol2$lik<minset$lik){
			if(trace) cat("\niter:",i+1,"  (linear): ", sol2$lik)
			minset$lik = sol2$lik
			minset$cf[c(idx_lin, idx_other)]  = sol2$coef
		} else{
			if(trace) cat("\niter:",i+1,"  (linear): no change")	
		}
		specx = spec
		solver.control = list(trace=0, method="BFGS", reltol=1e-8, maxit=4000, parsearch=TRUE)
		sol2 = .starfitx(specx, data, out.sample = out.sample, solver.control = solver.control,
				fit.control = fit.control, usepars=FALSE, transform.result=TRUE)
		if(sol2$convergence==0 && sol2$lik<minset$lik){
			if(trace) cat("\niter:",i+1," (restart): ", sol2$lik)
			minset$lik = sol2$lik
			minset$cf  = sol2$coef
		} else{
			if(trace) cat("\niter:",i+1," (restart): no change")	
		}
	}
	if(trace) cat("\n")
	solver.control = list(trace=trace, method="BFGS", reltol=1e-12, maxit=10000, parsearch=FALSE, usepars=FALSE)
	specx = spec
	setstart(specx)<-as.list(minset$cf)
	sol = starfit(specx, data, out.sample = out.sample, solver = "optim", solver.control = solver.control,
			fit.control = fit.control, cluster = NULL)
	return(sol)
}



.starfilter = function(spec, data, out.sample = 0, n.old = NULL, rec.init = 'all')
{
	tic = Sys.time()
	model = spec@model
	modelinc = model$modelinc
	fit.control=list()
	fit.control$stationarity = FALSE
	fit.control$rec.init = rec.init
	n.start = round(out.sample,0)
	# extract data and index
	xdata = .extractdata(data)
	if(!is.numeric(out.sample)) stop("\nstarfilter-->error: out.sample must be numeric\n")
	if(as.numeric(out.sample)<0) stop("\nstarfilter-->error: out.sample must be positive\n")
	n = length(xdata$data)
	
	#---------------------------------------------------------
	# check x (state regressors)
	#---------------------------------------------------------
	if(modelinc[47]==0){
		if(modelinc[49]==2){
			chk = all.equal(index(data), index(model$modeldata$s))
			if(!is.logical(chk) | chk == FALSE){
				print(paste("\n",chk,sep=""))
				stop("\nstarfilter-->error: data and x indices do not match\n")
			}
			XL = model$modeldata$s[1:(n-n.start), , drop = FALSE]
		} else{
			if(modelinc[48]==1){
				ytmp = xts(model$modeldata$fun(data[1:(n-n.start),1]), index(data[1:(n-n.start),1]))
			} else{
				ytmp = data[1:(n-n.start),1]
			}
			XL = NULL
			for(i in 1:length(model$modeldata$ylags)){
				if(i==1) XL = lag(ytmp, model$modeldata$ylags[i]) else XL = cbind(XL, lag(ytmp, model$modeldata$ylags[i]))
			}
			XL[is.na(XL)]=0
		}
		XL = coredata(XL)
	} else{
		XL = NULL
	}
	# XL is the lagged regressor dataset for the state dynamics
	#---------------------------------------------------------
	
	#---------------------------------------------------------
	# check fixed.probs (in case a set of fixed state probabilities
	# have been provided)
	#---------------------------------------------------------
	if(modelinc[47]==1){
		chk = all.equal(index(data), index(model$fixed.prob))
		if(!is.logical(chk) | chk == FALSE){
			print(paste("\n",chk,sep=""))
			stop("\nstarfilter-->error: data and fixed.probs indices do not match\n")
		}
		probs = coredata(model$fixed.prob[1:(n-n.start), , drop = FALSE])
	} else{
		probs = matrix(0, ncol = modelinc[46], nrow = (n-n.start))
	}
	#---------------------------------------------------------
	
	#---------------------------------------------------------
	# check ARX external regressors (X)
	#---------------------------------------------------------
	if(modelinc[3] > 0){
		chk = all(index(data)==index(model$modeldata$mexdata))
		if(!is.logical(chk) | chk == FALSE){
			print(paste("\n",chk,sep=""))
			stop("\nstarfit-->error: data and xreg indices do not match\n")
		}
		# lag mexdata according to xlags
		mexdata = coredata(model$modeldata$mexdata[1:(n-n.start), , drop = FALSE])
	} else{
		mexdata = NULL
	}
	
	if(modelinc[39] > 0){
		chk = all(index(data)==index(model$modeldata$vexdata))
		if(!is.logical(chk) | chk == FALSE){
			print(paste("\n",chk,sep=""))
			stop("\nstarfit-->error: data and vreg indices do not match\n")
		}
		# lag mexdata according to xlags
		vexdata = coredata(model$modeldata$vexdata[1:(n-n.start), , drop = FALSE])
	} else{
		vexdata = NULL
	}
	#---------------------------------------------------------
	# index the data to start at point n-out.sample
	data = xdata$data[1:(n-n.start)]
	index = xdata$index[1:(n-n.start)]
	origdata = xdata$data
	origindex = xdata$index
	period = xdata$period
	if(!is.null(n.old)) Nx = n.old else Nx = length(data)
	#---------------------------------------------------------
	# create a temporary environment to store values (deleted at end of function)
	# and the list of parameters/data passed to the likelihood estimation functions
	starenv = new.env(hash = TRUE)
	arglist = list()
	arglist$starenv <- starenv
	# parallel mode flag (related to use of starenv environment)
	arglist$pmode = 0
	pidx = model$pidx
	arglist$index = index
	arglist$trace = 0
	arglist$N = Nx
	arglist$fit.control = list()
	model$modeldata$T = T = length(as.numeric(data))
	recinit = .checkrec(rec.init, Nx)
	arglist$recinit = recinit
	dist = model$modeldesc$distribution
	arglist$data = data
	arglist$XL = XL
	# in the ugarchfilter, N acts as n.old to avoid lookahead bias in filtering
	# for initialization values
	arglist$mexdata = mexdata
	arglist$vexdata = vexdata
	arglist$probs = probs
	arglist$model = model
	arglist$ipars = model$pars
	# we now split out any fixed parameters
	estidx = as.logical( model$pars[,3] )
	arglist$estidx = estidx	
	assign("star_llh", 1, envir = starenv)
	arglist$fit.control = fit.control
	arglist$transform = 0
	arglist$returnType ="all"
	# need to control for use of parallel environment which cannot handle the
	# calls between the locally stored environment
	if(modelinc[50]>0){
		if(modelinc[50]==4){
			llfun = switch(modelinc[46],
					.stars2dLLH2mix,
					.stars2dLLH2mix,
					.stars3dLLH3mix,
					.stars4dLLH4mix)
		} else{
			llfun = switch(modelinc[46],
					.stars1dLLH,
					.stars2dLLH,
					.stars3dLLH,
					.stars4dLLH)
		}
	} else{
		llfun = switch(modelinc[46],
				.stars1sLLH,
				.stars2sLLH,
				.stars3sLLH,
				.stars4sLLH)
	}
	pars = model$pars[estidx,1]
	temp = llfun(pars, arglist)
	
	filter = list()
	filter$z = temp$z
	filter$LLH = -temp$llh
	filter$log.likelihoods = temp$LHT
	filter$residuals = temp$res
	filter$condm = temp$condm
	filter$constm = temp$constm
	filter$probability = temp$probs
	filter$pmu = temp$pmu
	if(sum(modelinc[31:40])>0){
		if(modelinc[50]==4) filter$sigma = temp$h else filter$sigma = sqrt(temp$h)
	}
	filter$ipars = model$pars
	model$modeldata$data = origdata
	model$modeldata$index = origindex
	model$modeldata$period = period
	model$n.start = out.sample	
	sol = new("STARfilter",
			filter = filter,
			model = model)
	return(sol)
}


# if mc.sims= NULL use 12*N
# CONTINUE HERE
.starforecast.fit = function(fitORspec, data = NULL, n.ahead = 1, n.roll = 0, out.sample = 0, 
		external.forecasts = list(xregfor = NULL, vregfor = NULL, sfor = NULL, probfor = NULL), 
		method = c("an.parametric", "an.kernel", "mc.empirical", "mc.parametric", "mc.kernel"), 
		mc.sims = NULL, ...)
{
	fit    = fitORspec
	data   = fit@model$modeldata$data
	Nor    = length(as.numeric(data))
	index  = fit@model$modeldata$index
	period = fit@model$modeldata$period
	ns = fit@model$n.start
	N = Nor - ns
	model = fit@model
	ipars = fit@fit$ipars
	modelinc = model$modelinc
	if(n.ahead==1) method="n.ahead-1" else method = method[1]
	idx = model$pidx
	if( n.roll > ns ) stop("\nstarforecast-->error: n.roll must not be greater than out.sample!")
	pars = fit@fit$coef
	ipars = fit@fit$ipars
	# check if necessary the external regressor forecasts provided first
	xreg = .forcregressors(model, external.forecasts$xregfor, 
			external.forecasts$vregfor, external.forecasts$sfor, 
			external.forecasts$probfor, n.ahead, Nor, out.sample = ns, n.roll)
	if(modelinc[50]>0) dynamic = TRUE else dynamic = FALSE
	# need 1 extra for case that n.roll==out.sample+1 (extend the index)
	fcreq = ifelse(ns >= (n.ahead+n.roll), n.ahead+n.roll, ns+1)
	pindex = as.POSIXct(index)
	newindex = c(pindex[1:(N + fcreq-1)], generatefwd(pindex[N+fcreq-1], length.out = 1, by = period))
	
	# This is only for the 1-ahead filter
	if(modelinc[3]>0) mxfx = xts(xreg$mxf[1:(N + fcreq),,drop=FALSE], newindex[1:(N + fcreq)]) else mxfx = NULL
	if(modelinc[39]>0) vxfx = xts(xreg$vxf[1:(N + fcreq),,drop=FALSE], newindex[1:(N + fcreq)]) else vxfx = NULL
	if(modelinc[20]>0 && modelinc[49]==2) sxfx = xts(xreg$sxf[1:(N + fcreq),,drop=FALSE], newindex[1:(N + fcreq)]) else sxfx = NULL
	if(modelinc[47]>0) pxfx = xts(xreg$pxf[1:(N + fcreq),,drop=FALSE], newindex[1:(N + fcreq)]) else pxfx = NULL
	if(modelinc[50]==2){
		if(modelinc[32]>0) kappa = pneg(ipars[idx["ghlambda",1],1], ipars[idx["shape",1],1], ipars[idx["skew",1],1], model$modeldesc$distribution)
	} else if(modelinc[50]==3){
		kappa = egarchKappa(ipars[idx["ghlambda",1],1], ipars[idx["shape",1],1], ipars[idx["skew",1],1], model$modeldesc$distribution)
	} else{
		kappa = 1
	}
	if(dynamic) vmodel = c("sGARCH","gjrGARCH","eGARCH","mixture")[modelinc[50]] else vmodel = NULL
	if(dynamic && modelinc[50]<4) usegarch = TRUE else usegarch = FALSE 
	
	fspec = starspec(mean.model = list(
					states = modelinc[46],
					include.intercept = modelinc[c(1,5,9,13)], 
					arOrder = model$modeldata$arOrder,
					maOrder = model$modeldata$maOrder,
					matype = ifelse(modelinc[17]>0, "linear","state"),
					statevar = c("y","s")[modelinc[49]],
					s = sxfx, statear = model$modeldesc$statear,
					ylags = model$modeldata$ylags, 
					xreg = mxfx, yfun = model$modeldata$fun),
			variance.model = list(dynamic = dynamic, 
					model = vmodel, garchOrder = c(modelinc[32], modelinc[33]), 
					submodel = NULL, vreg = vxfx),  
			distribution.model = model$modeldesc$distribution, fixed.pars = as.list(pars),
			fixed.probs = pxfx)
	
	tmp =  xts(c(data[1:(N + fcreq-1)],0), newindex)
	flt = starfilter(data = tmp, spec = fspec, n.old = N)
	if(dynamic) sigmafilter = flt@filter$sigma else sigmafilter = NULL
	resfilter = flt@filter$residuals
	zfilter = flt@filter$z
	pmufilter = flt@filter$pmu
	probfilter = flt@filter$probability
	condmfilter = flt@filter$condm
	constmfilter = flt@filter$constm
	
	sigmaFor = seriesFor =  matrix(NA, ncol = n.roll+1, nrow = n.ahead)
	condmFor = probFor = array(NA, dim=c(n.ahead, n.roll+1, modelinc[46]), 
			dimnames = list(paste("T+", 1:n.ahead, sep=""), as.character(index[N:(N+n.roll)]), paste("State_",1:modelinc[46],sep="")))
	colnames(seriesFor) = colnames(sigmaFor) = as.character(index[N:(N+n.roll)])
	rownames(seriesFor) = rownames(sigmaFor) = paste("T+", 1:n.ahead, sep="")
	# The 1-ahead
	if(dynamic) sigmaFor[1,] = tail(sigmafilter, n.roll+1) else sigmaFor = NULL
	seriesFor[1,] = tail(rowSums(condmfilter*probfilter), n.roll+1)
	for(i in 1:modelinc[46]){
		probFor[1,,i] = tail(probfilter[,i], n.roll+1)
		condmFor[1,,i] = tail(condmfilter[,i], n.roll+1)
	}
	esim = ydist = NULL
	
	if(n.ahead>1){
		arglist = list()
		arglist$ipars = ipars
		arglist$modelinc = modelinc
		arglist$ylags = model$modeldata$ylags
		arglist$distribution = model$modeldesc$distribution
		arglist$idx = idx
		arglist$n.ahead = n.ahead
		arglist$kappa = kappa
		arglist$method = method
		arglist$mc.sims = mc.sims
		if(modelinc[48]>0) arglist$yfun = model$modeldata$fun
		arglist$gfun = switch(modelinc[46],"starxsim1", "starxsim2","starxsim3","starxsim4")
		if(method!="analytic"){
			# the distribution of the returns using the mc methods
			# each slot in the list will contain a matrix of dim=c(sims, n.ahead)
			esim = ydist = vector(mode="list", length=n.roll+1)
		} else{
			esim = ydist = NULL
		}
		pxfi = mxfi = vxfi = xxfi = NULL
		for(i in 1:(n.roll+1)){
			arglist$np = N + i - 1
			if(usegarch) arglist$omega = rep(ipars[idx["omega",1]:idx["omega",2], 1], N+i+n.ahead-1)
			if(dynamic) arglist$h = c(sigmafilter[1:(N+i-1)], rep(0, n.ahead))
			arglist$epsx = c(resfilter[1:(N+i-1)], rep(0, n.ahead))
			# include the filtered n.ahead=1
			arglist$y = c(data[1:(N+i-1)], seriesFor[1,i], rep(0, n.ahead-1))
			arglist$z = c(zfilter[1:(N+i-1)], rep(0, n.ahead))
			arglist$constm = constmfilter[1:(N+i-1),,drop=FALSE]
			# need to include in pmu the 1-ahead value
			arglist$pmu = pmufilter[1:(N+i),,drop=FALSE]
			arglist$probs = probfilter[1:(N+i),,drop=FALSE]
			arglist$residuals = resfilter
			arglist$zresiduals = zfilter
			# forecast of externals is provided outside the system
			if(modelinc[3]>0)  arglist$mxfi = xreg$mxf[1:(N+i-1+n.ahead), , drop = FALSE]
			if(modelinc[39]>0) arglist$vxfi = xreg$vxf[1:(N+i-1+n.ahead), , drop = FALSE]
			if(modelinc[47]>0) arglist$pxfi = xreg$pxf[1:(N+i-1+n.ahead), , drop = FALSE]
			if(modelinc[49]==2) arglist$sxfi = xreg$sxf[1:(N+i-1+n.ahead), , drop = FALSE]
			if(usegarch){
				ans = switch(vmodel, 
						"sGARCH"   = .nsgarchforecast(arglist), 
						"gjrGARCH" = .ngjrgarchforecast(arglist), 
						"eGARCH"   = .negarchforecast(arglist))
				sigmaFor[,i] = tail(ans$h, n.ahead)
			} else{
				if(modelinc[50]==4){
					ans = starfmix(arglist)
					sigmaFor[2:n.ahead,i] = tail(ans$h, n.ahead-1)
				} else{
					ans = starf(arglist)
				}
			}
			seriesFor[2:n.ahead,i] = tail(ans$y, n.ahead-1)
			for(k in 1:modelinc[46]) probFor[2:n.ahead,i,k] = tail(ans$probs[,k], n.ahead-1)
			ydist[[i]] = ans$ydist
			esim[[i]] = ans$esim
		}
		
	}
	fcst = list()
	fcst$n.ahead = n.ahead
	fcst$N = N+ns
	fcst$n.start = ns
	fcst$n.roll = n.roll
	fcst$yDist = ydist
	fcst$probFor = probFor
	fcst$eSim = esim
	if(dynamic) fcst$sigmaFor = sigmaFor
	fcst$seriesFor = seriesFor
	if(dynamic) model$modeldata$sigma = flt@filter$sigma
	model$modeldata$residuals = flt@filter$residuals
	fcst$method = method
	ans = new("STARforecast",
			forecast = fcst,
			model = model)
	return(ans)
}

.starforecast.spec = function(fitORspec, data = NULL, n.ahead = 1, n.roll = 0, out.sample = 0, 
		external.forecasts = list(xregfor = NULL, vregfor = NULL, sfor = NULL, probfor = NULL), 
		method = c("an.parametric", "an.kernel", "mc.empirical", "mc.parametric", "mc.kernel"), 
		mc.sims = NULL, ...)
{
	spec    = fitORspec
	if(is.null(data)) stop("\nstarforecast-->error: data must not be NULL when using a specification!")
	if(n.ahead==1) method="n.ahead-1" else method = method[1]
	xdata = .extractdata(data)
	Nor = length(as.numeric(xdata$data))
	data = xdata$data
	N = length(as.numeric(data))
	index = xdata$index
	period = xdata$period
	ns = out.sample
	if( n.roll > ns ) stop("\nstarforecast-->error: n.roll must not be greater than out.sample!")
	N = Nor - ns
	model = spec@model
	ipars = model$pars
	pars = unlist(model$fixed.pars)
	parnames = names(pars)
	modelnames = .checkallfixed(spec)
	if(is.na(all(match(modelnames, parnames), 1:length(modelnames)))) {
		cat("\nstarforecast-->error: parameters names do not match specification\n")
		cat("Expected Parameters are: ")
		cat(paste(modelnames))
		cat("\n")
		stop("Exiting", call. = FALSE)
	}
	# once more into the spec
	setfixed(spec)<-as.list(pars)
	model = spec@model
	idx = model$pidx
	ipars = model$pars
	modelinc = model$modelinc
	model$modeldata$data = data
	model$modeldata$index = index
	model$modeldata$period = period
	if(modelinc[50]>0) dynamic = TRUE else dynamic = FALSE
	
	# check if necessary the external regressor forecasts provided first
	xreg = .forcregressors(model, external.forecasts$xregfor, 
			external.forecasts$vregfor, external.forecasts$sfor, 
			external.forecasts$probfor, n.ahead, Nor, out.sample = ns, n.roll)
	# need 1 extra for case that n.roll==out.sample+1 (extend the index)
	fcreq = ifelse(ns >= (n.ahead+n.roll), n.ahead+n.roll, ns+1)
	pindex = as.POSIXct(index)
	newindex = c(pindex[1:(N + fcreq-1)], generatefwd(pindex[N+fcreq-1], length.out = 1, by = period))
	
	# This is only for the 1-ahead filter
	if(modelinc[3]>0) mxfx = xts(xreg$mxf[1:(N + fcreq),,drop=FALSE], newindex[1:(N + fcreq)]) else mxfx = NULL
	if(modelinc[39]>0) vxfx = xts(xreg$vxf[1:(N + fcreq),,drop=FALSE], newindex[1:(N + fcreq)]) else vxfx = NULL
	if(modelinc[20]>0 && modelinc[49]==2) sxfx = xts(xreg$sxf[1:(N + fcreq),,drop=FALSE], newindex[1:(N + fcreq)]) else sxfx = NULL
	if(modelinc[47]>0) pxfx = xts(xreg$pxf[1:(N + fcreq),,drop=FALSE], newindex[1:(N + fcreq)]) else pxfx = NULL
	if(modelinc[50]==2){
		if(modelinc[32]>0) kappa = pneg(ipars[idx["ghlambda",1],1], ipars[idx["shape",1],1], ipars[idx["skew",1],1], model$modeldesc$distribution)
	} else if(modelinc[50]==3){
		kappa = egarchKappa(ipars[idx["ghlambda",1],1], ipars[idx["shape",1],1], ipars[idx["skew",1],1], model$modeldesc$distribution)
	} else{
		kappa = 1
	}
	if(dynamic) vmodel = c("sGARCH","gjrGARCH","eGARCH","mixture")[modelinc[50]] else vmodel = NULL
	if(dynamic && modelinc[50]<4) usegarch = TRUE else usegarch = FALSE 
	
	fspec = starspec(mean.model = list(
					states = modelinc[46],
					include.intercept = modelinc[c(1,5,9,13)], 
					arOrder = model$modeldata$arOrder,
					maOrder = model$modeldata$maOrder,
					matype = ifelse(modelinc[17]>0, "linear","state"),
					statevar = c("y","s")[modelinc[49]],
					s = sxfx, statear = model$modeldesc$statear,
					ylags = model$modeldata$ylags, 
					xreg = mxfx, yfun = model$modeldata$fun),
			variance.model = list(dynamic = dynamic, 
					model = vmodel, garchOrder = c(modelinc[32], modelinc[33]), 
					submodel = NULL, vreg = vxfx),  
			distribution.model = model$modeldesc$distribution, fixed.pars = as.list(pars),
			fixed.probs = pxfx)
	
	tmp =  xts(c(data[1:(N + fcreq-1)],0), newindex)
	flt = starfilter(data = tmp, spec = fspec, n.old = N)
	if(dynamic) sigmafilter = flt@filter$sigma else sigmafilter = NULL
	resfilter = flt@filter$residuals
	zfilter = flt@filter$z
	pmufilter = flt@filter$pmu
	probfilter = flt@filter$probability
	condmfilter = flt@filter$condm
	constmfilter = flt@filter$constm
	
	sigmaFor = seriesFor =  matrix(NA, ncol = n.roll+1, nrow = n.ahead)
	condmFor = probFor = array(NA, dim=c(n.ahead, n.roll+1, modelinc[46]), 
			dimnames = list(paste("T+", 1:n.ahead, sep=""), as.character(index[N:(N+n.roll)]), paste("State_",1:modelinc[46],sep="")))
	colnames(seriesFor) = colnames(sigmaFor) = as.character(index[N:(N+n.roll)])
	rownames(seriesFor) = rownames(sigmaFor) = paste("T+", 1:n.ahead, sep="")
	# The 1-ahead
	if(dynamic) sigmaFor[1,] = tail(sigmafilter, n.roll+1) else sigmaFor = NULL
	seriesFor[1,] = tail(rowSums(condmfilter*probfilter), n.roll+1)
	for(i in 1:modelinc[46]){
		probFor[1,,i] = tail(probfilter[,i], n.roll+1)
		condmFor[1,,i] = tail(condmfilter[,i], n.roll+1)
	}
	esim = ydist = NULL
	if(n.ahead>1){
		arglist = list()
		arglist$ipars = ipars
		arglist$modelinc = modelinc
		arglist$ylags = model$modeldata$ylags
		arglist$distribution = model$modeldesc$distribution
		arglist$idx = idx
		arglist$n.ahead = n.ahead
		arglist$kappa = kappa
		arglist$method = method
		arglist$mc.sims = mc.sims
		if(modelinc[48]>0) arglist$yfun = model$modeldata$fun
		arglist$gfun = switch(modelinc[46],"starxsim1", "starxsim2","starxsim3","starxsim4")
		if(method!="analytic"){
			# the distribution of the returns using the mc methods
			# each slot in the list will contain a matrix of dim=c(sims, n.ahead)
			esim = ydist = vector(mode="list", length=n.roll+1)
		} else{
			esim = ydist = NULL
		}
		pxfi = mxfi = vxfi = xxfi = NULL
		for(i in 1:(n.roll+1)){
			arglist$np = N + i - 1
			if(usegarch) arglist$omega = rep(ipars[idx["omega",1]:idx["omega",2], 1], N+i+n.ahead-1)
			if(usegarch) arglist$h = c(sigmafilter[1:(N+i-1)], rep(0, n.ahead))
			arglist$epsx = c(resfilter[1:(N+i-1)], rep(0, n.ahead))
			# include the filtered n.ahead=1
			arglist$y = c(data[1:(N+i-1)], seriesFor[1,i], rep(0, n.ahead-1))
			arglist$z = c(zfilter[1:(N+i-1)], rep(0, n.ahead))
			arglist$constm = constmfilter[1:(N+i-1),,drop=FALSE]
			# need to include in pmu the 1-ahead value
			arglist$pmu = pmufilter[1:(N+i),,drop=FALSE]
			arglist$probs = probfilter[1:(N+i),,drop=FALSE]
			arglist$residuals = resfilter
			# forecast of externals is provided outside the system
			if(modelinc[3]>0)  arglist$mxfi = xreg$mxf[1:(N+i-1+n.ahead), , drop = FALSE]
			if(modelinc[39]>0) arglist$vxfi = xreg$vxf[1:(N+i-1+n.ahead), , drop = FALSE]
			if(modelinc[47]>0) arglist$pxfi = xreg$pxf[1:(N+i-1+n.ahead), , drop = FALSE]
			if(modelinc[49]==2) arglist$sxfi = xreg$sxf[1:(N+i-1+n.ahead), , drop = FALSE]
			if(usegarch){
				ans = switch(vmodel, 
						"sGARCH"   = .nsgarchforecast(arglist), 
						"gjrGARCH" = .ngjrgarchforecast(arglist), 
						"eGARCH"   = .negarchforecast(arglist))
				sigmaFor[,i] = tail(ans$h, n.ahead)
			} else{
				if(modelinc[50]==4){
					ans = starfmix(arglist)
					sigmaFor[2:n.ahead,i] = tail(ans$h, n.ahead-1)
				} else{
					ans = starf(arglist)
				}
			}
			seriesFor[2:n.ahead,i] = tail(ans$y, n.ahead-1)
			for(k in 1:modelinc[46]) probFor[2:n.ahead,i,k] = tail(ans$probs[,k], n.ahead-1)
			ydist[[i]] = ans$ydist
			esim[[i]] = ans$esim
		}
	}
	fcst = list()
	fcst$n.ahead = n.ahead
	fcst$N = N+ns
	fcst$n.start = ns
	fcst$n.roll = n.roll
	fcst$yDist = ydist
	fcst$probFor = probFor
	fcst$eSim = esim
	if(dynamic) fcst$sigmaFor = sigmaFor
	fcst$seriesFor = seriesFor
	if(dynamic) model$modeldata$sigma = flt@filter$sigma
	model$modeldata$residuals = flt@filter$residuals
	fcst$method = method
	ans = new("STARforecast",
			forecast = fcst,
			model = model)
	return(ans)
}


.starsim.dynamic = function(fit, n.sim = 1000, n.start = 0, m.sim = 1, presigma = NA, 
		prereturns = NA, preresiduals = NA, rseed = NA, 
		custom.dist = list(name = NA, distfit = NA), xregsim = NULL, 
		vregsim = NULL, ssim = NULL, probsim = NULL)
{
	# some checks
	if(is.na(rseed[1])){
		sseed = as.integer(runif(1,0,as.integer(Sys.time())))
	} else{
		if(length(rseed) != m.sim) sseed = as.integer(rseed[1]) else sseed = rseed[1:m.sim]
	}
	# Enlarge Series:
	# need to allow for arfima case:
	n = n.sim + n.start
	data = fit@model$modeldata$data
	N = length(as.numeric(data))
	data = data[1:(N - fit@model$n.start)]
	N = length(as.numeric(data))
	resids = fit@fit$residuals	
	sig = fit@fit$sigma
	pmu = fit@fit$pmu
	probs = fit@fit$probability
	condm = fit@fit$condm
	model = fit@model
	modelinc = model$modelinc
	yfun = model$modeldata$fun
	# max of AR & MA terms and STAR dynamics in the case that y is used
	mar = max(c(modelinc[c(2,6,10,14,4,8,12,16,17)], ifelse(modelinc[49]==1, max(model$modeldata$ylags), 0)))
	m = max(c(modelinc[32:33], mar))
	
	idx = model$pidx
	ipars = fit@fit$ipars
	ylags = model$modeldata$ylags
	# check if necessary the external regressor forecasts provided first
	xreg = .simregressors(model, xregsim, vregsim, ssim, N, n, m.sim)
	mexsim = xreg$mexsimlist
	vexsim = xreg$vexsimlist
	ssim   = xreg$ssimlist
	# Random Samples from the Distribution
	if(length(sseed) == 1){
		set.seed(sseed)
		z = .custzdist(custom.dist = custom.dist, distribution = model$modeldesc$distribution, n = n, m.sim = m.sim, 
				mu = 0, sigma = 1, lambda = ipars[idx["ghlambda",1], 1], 
				skew = ipars[idx["skew",1], 1], shape = ipars[idx["shape",1],1])
	} else{
		z = matrix(0, ncol = m.sim, nrow = n)
		for(i in 1:m.sim){
			set.seed(sseed[i])
			z[,i] = .custzdist(custom.dist = custom.dist, distribution = model$modeldesc$distribution, n = n, m.sim = 1, 
					mu = 0, sigma = 1, lambda = ipars[idx["ghlambda",1], 1], skew = ipars[idx["skew",1], 1], 
					shape = ipars[idx["shape",1],1])
		}
	}
	z = rbind(matrix(tail(fit@fit$z, m), nrow = m, ncol = m.sim), z)
	
	# create the presample information
	if(!is.na(presigma[1])){
		presigma = as.vector(presigma)
		if(length(presigma)<m) stop(paste("\nstarsim-->error: presigma must be of length ", m, sep=""))
	}
	if(!is.na(prereturns[1])){
		prereturns = as.vector(prereturns)
		if(length(prereturns)<mar) stop(paste("\nstarsim-->error: prereturns must be of length ", mar, sep=""))
	}
	if(!is.na(preresiduals[1])){
		preresiduals = as.vector(preresiduals)
		if(length(preresiduals)<mar) stop(paste("\nstarsim-->error: preresiduals must be of length ", mar, sep=""))
		preres = preresiduals
	} else{
		preres = as.numeric(tail(residuals(fit), m))
	}
	if(is.na(presigma[1])){
		presigma  = tail(sig, m)
	}
	if(is.na(prereturns[1])){
		prereturns = tail(data, mar)
	}
	
	# input vectors/matrices
	# outpus matrices
	sigmaSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	seriesSim = matrix(0, ncol = m.sim, nrow = n.sim)
	residSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	condmSim = probSim  =  vector(mode="list", length = m.sim)
	z[is.na(z) | is.nan(z) | !is.finite(z)] = 0
	constm = matrix(0, ncol = modelinc[46], nrow = n)
	# 2 cases for the star model:
	# case 1: probability determined by external regressors in which case
	# the simulation is quite fast
	# case 2: probability determined by lag(fun(y)) in which case much slower...
	# ...need to consider re-writing probability dynamics in C and also calling an
	# user defined functions from C.
	arglist = list()
	arglist$ipars = ipars
	arglist$idx = idx
	arglist$pmu = pmu
	arglist$model = model
	arglist$probs = probs
	
	for(i in 1:m.sim){
		h = c(presigma^2, rep(0, n))
		x = c(prereturns, rep(0, n))
		
		if(modelinc[50]==2){
			ngrd = which(preres<0)
			tmpr = rep(0, length(preres))
			tmpr[ngrd] = preres[ngrd] * preres[ngrd]
			nres = c(tmpr,   rep(0, n))
		} else{
			nres = 0
		}
		if(modelinc[50]==3){
			kappa = egarchKappa(ipars[idx["ghlambda",1],1], ipars[idx["shape",1],1], ipars[idx["skew",1],1], model$modeldesc$distribution)	
		} else{
			kappa = 1
		}
		res = c(preres, rep(0, n))
		
		ans1 = try(.C("garchsimC", model = as.integer(modelinc[1:60]), pars = as.double(ipars[,1]), 
						idx = as.integer(idx[,1]-1), h = as.double(h), z = as.double(z[,i]), 
						res = as.double(res), e = as.double(res*res), nres = as.double(nres),
						meanz = as.double(kappa), vexdata = as.double(vexsim[[i]]), 
						T = as.integer(n+m), m = as.integer(m), PACKAGE = "twinkle"), silent = TRUE)
		if(inherits(ans1, "try-error")) stop("\nstarsim-->error: error in calling C function....\n")
		
		sigmaSim[,i] = ans1$h[(n.start + m + 1):(n+m)]^(1/2)
		residSim[,i] = ans1$res[(n.start + m + 1):(n+m)]
		# since m = max(garchOrder, mar)
		simres = tail(ans1$res,n+mar)
		for(j in 1:modelinc[46]){
			if(modelinc[c(1,5,9,13)[j]]>0) constm[,j] = as.numeric(ipars[idx[paste("s",j,".phi0",sep=""),1],1])
		}
		if(modelinc[3]>0){
			for(j in 1:modelinc[46]) constm[,j] = constm[,j] + matrix( ipars[idx[paste("s",j,".xi",sep=""),1]:idx[paste("s",j,".xi",sep=""),2], 1], ncol = modelinc[3] ) %*% t( matrix( mexsim[[i]], ncol = modelinc[3] ) )
		}
		constm = rbind(matrix(0, ncol=modelinc[46], nrow=mar), constm)
		if(modelinc[49]==2){
			arglist$XL = ssim[[i]]
			# sim starts at (m=mar) and end at n+mar
			# get probabilities
			ptmp = switch(modelinc[46], dstar1sim(arglist), dstar2sim(arglist), dstar3sim(arglist), dstar4sim(arglist))
			rsim = switch(modelinc[46],
					try(.C("starxsim1", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), x = as.double(x), 
									s = double(modelinc[46]*(n+mar)),res = as.double(simres), prob = as.double(ptmp$probs), constm = as.double(constm), 
									m = as.integer(mar), T = as.integer(mar+n), PACKAGE="twinkle"), silent = TRUE),
					try(.C("starxsim2", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), x = as.double(x), 
									s = double(modelinc[46]*(n+mar)),res = as.double(simres), prob = as.double(ptmp$probs), constm = as.double(constm), 
									m = as.integer(mar), T = as.integer(mar+n), PACKAGE="twinkle"), silent = TRUE),
					try(.C("starxsim3", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), x = as.double(x), 
									s = double(modelinc[46]*(n+mar)),res = as.double(simres), prob = as.double(ptmp$probs), constm = as.double(constm), 
									m = as.integer(mar), T = as.integer(mar+n), PACKAGE="twinkle"), silent = TRUE),
					try(.C("starxsim4", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), x = as.double(x), 
									s = double(modelinc[46]*(n+mar)),res = as.double(simres), prob = as.double(ptmp$probs), constm = as.double(constm), 
									m = as.integer(mar), T = as.integer(mar+n), PACKAGE="twinkle"), silent = TRUE)
					)
			seriesSim[,i] = tail(rsim$x, n.sim)
			probSim[[i]]  = tail(ptmp$probs, n.sim)
			condmSim[[i]]  = tail(matrix(rsim$s, ncol = modelinc[46]), n.sim)
			colnames(condmSim[[i]]) = colnames(probSim[[i]]) = paste("state[",1:modelinc[46],"]",sep="")
			rownames(condmSim[[i]]) = rownames(probSim[[i]]) = paste("T+",1:n.sim,sep="")			
		} else{
			# NOTE: would probably benefit from using .Call and performing the whole process in C++ with a call to R for "yfun"...
			ndx = length(data)
			y = c(data, rep(0, n))
			y[(N-mar+1):N] = prereturns
			
			psim = matrix(0, ncol = modelinc[46], nrow = mar+n)
			xcondm = matrix(0, ncol = modelinc[46], nrow = mar+n)
			for(j in 1:n){
				if(modelinc[48]==1) ytmp = c(yfun(as.numeric(y[1:(ndx+j-1)])), NA) else ytmp = c(as.numeric(y[1:(ndx+j-1)]),NA)		
				XL = matrix(NA, ncol = length(ylags), nrow = ndx+j)
				nXL = nrow(XL)
				for(k in 1:ncol(XL)){
					XL[(ylags[k]+1):nXL,k] = ytmp[1:(nXL-ylags[k])]
				}
				XL[is.na(XL)]=0
				arglist$XL = XL
				ptmp = switch(modelinc[46], dstar1sim(arglist), dstar2sim(arglist), dstar3sim(arglist), dstar4sim(arglist))
				psim[mar+j,] = tail(ptmp$probs, 1)
				rsim = switch(modelinc[46],
						try(.C("starxsim1", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
										x = as.double(x), s = double(modelinc[46]*(mar+j)), res = as.double(simres[1:(mar+j)]),
										prob = as.double(psim[1:(mar+j),,drop=FALSE]), constm = as.double(constm[1:(mar+j),,drop=FALSE]),
										m = as.integer(mar+j-1), T = as.integer(mar+j), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim2", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
										x = as.double(x), s = double(modelinc[46]*(mar+j)), res = as.double(simres[1:(mar+j)]),
										prob = as.double(psim[1:(mar+j),,drop=FALSE]), constm = as.double(constm[1:(mar+j),,drop=FALSE]),
										m = as.integer(mar+j-1), T = as.integer(mar+j), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim3", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
										x = as.double(x), s = double(modelinc[46]*(mar+j)), res = as.double(simres[1:(mar+j)]),
										prob = as.double(psim[1:(mar+j),,drop=FALSE]), constm = as.double(constm[1:(mar+j),,drop=FALSE]),
										m = as.integer(mar+j-1), T = as.integer(mar+j), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim4", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
										x = as.double(x), s = double(modelinc[46]*(mar+j)), res = as.double(simres[1:(mar+j)]),
										prob = as.double(psim[1:(mar+j),,drop=FALSE]), constm = as.double(constm[1:(mar+j),,drop=FALSE]),
										m = as.integer(mar+j-1), T = as.integer(mar+j), PACKAGE="twinkle"), silent = TRUE)
				)
				x[mar+j] = rsim$x[mar+j]
				y[ndx+j] = rsim$x[mar+j]
				xcondm[mar+j,] = tail(matrix(rsim$s, ncol = modelinc[46]),1)
			}
			seriesSim[,i] = tail(x, n.sim)
			probSim[[i]]  = tail(ptmp$probs, n.sim)
			condmSim[[i]] = tail(xcondm, n.sim)
			colnames(probSim[[i]]) = paste("state[",1:modelinc[46],"]",sep="")
			rownames(probSim[[i]]) = paste("T+",1:n.sim,sep="")
		}
	}
	residSim = tail(residSim, n.sim)
	rownames(residSim) = rownames(seriesSim) = rownames(sigmaSim) = paste("T+",1:n.sim,sep="")
	# check:
	# u=3
	# xtmp = (ipars["s1.phi0",1] + ipars["s1.phi1",1]*x[u-1]+ipars["s1.phi2",1]*x[u-2])*psim[u,1]+ (ipars["s2.phi0",1] + ipars["s2.phi1",1]*x[u-1]+ipars["s2.phi2",1]*x[u-2])*(psim[u,2]) + simres[u]
	sim = list(sigmaSim = sigmaSim, seriesSim = seriesSim, residSim = residSim, probSim = probSim)
	model$modeldata$sigma = sigma
	sol = new("STARsim",
			simulation = sim,
			model = model,
			seed = as.integer(sseed))
	return(sol)
}

# CONTINUE HERE: Change code to include 3 and 4 mixtures
.starsim.mixture = function(fit, n.sim = 1000, n.start = 0, m.sim = 1, prereturns = NA, preresiduals = NA, rseed = NA, 
		custom.dist = list(name = NA, distfit = NA), xregsim = NULL, ssim = NULL, probsim = NULL)
{
	# some checks
	if(is.na(rseed[1])){
		sseed = as.integer(runif(1,0,as.integer(Sys.time())))
	} else{
		if(length(rseed) != m.sim) sseed = as.integer(rseed[1]) else sseed = rseed[1:m.sim]
	}
	# Enlarge Series:
	# need to allow for arfima case:
	n = n.sim + n.start
	data = fit@model$modeldata$data
	N = length(as.numeric(data))
	data = data[1:(N - fit@model$n.start)]
	N = length(as.numeric(data))
	resids = fit@fit$residuals	
	modelinc = fit@model$modelinc
	
	sig1 = coef(fit)["s1.sigma"]
	sig2 = coef(fit)["s2.sigma"]
	if(modelinc[46]==3)	sig3 = coef(fit)["s3.sigma"] else sig3 = NULL
	if(modelinc[46]==4)	sig4 = coef(fit)["s4.sigma"] else sig4 = NULL
	
	pmu = fit@fit$pmu
	probs = fit@fit$probability
	condm = fit@fit$condm
	model = fit@model
	yfun = model$modeldata$fun
	# max of AR terms and STAR dynamics in the case that y is used
	mar = max(c(modelinc[c(2,6,10,14,4,8,12,16,17)], ifelse(modelinc[49]==1, max(model$modeldata$ylags), 0)))
	idx = model$pidx
	ipars = fit@fit$ipars
	ylags = model$modeldata$ylags
	gfun = switch(modelinc[46],"starxsim1", "starxsim2","starxsim3","starxsim4")
	# check if necessary the external regressor forecasts provided first
	xreg = .simregressors(model, xregsim, NULL, ssim, N, n, m.sim)	
	mexsim = xreg$mexsimlist
	ssim   = xreg$ssimlist
	# Random Samples from the Distribution
	# will postmultiply by mixture sigma once we have probs
	if(!is.na(custom.dist$name)){
		zresidSim = custom.dist$distfit
	} else{
		if(length(sseed) == 1){
			set.seed(sseed)
			zresidSim = matrix(rdist(model$modeldesc$distribution, n*m.sim, mu = 0, sigma = 1, lambda = ipars[idx["ghlambda",1], 1], 
							skew = ipars[idx["skew",1], 1], shape = ipars[idx["shape",1],1]), ncol = m.sim, nrow = n)
		} else{
			zresidSim = matrix(0, ncol = m.sim, nrow = n)
			for(i in 1:m.sim){
				set.seed(sseed[i])
				zresidSim[,i] = rdist(model$modeldesc$distribution, n, mu = 0, sigma = 1, lambda = ipars[idx["ghlambda",1], 1], 
						skew = ipars[idx["skew",1], 1], shape = ipars[idx["shape",1],1])
			}
		}
	}

	# create the presample information
	if(!is.na(prereturns[1])){
		prereturns = as.vector(prereturns)
		if(length(prereturns)<mar) stop(paste("\nstarsim-->error: prereturns must be of length ", mar, sep=""))
	}
	if(!is.na(preresiduals[1])){
		preresiduals = as.vector(preresiduals)
		if(length(preresiduals)<mar) stop(paste("\nstarsim-->error: preresiduals must be of length ", mar, sep=""))
		preres = preresiduals
	} else{
		preres = as.numeric(tail(residuals(fit), mar))
	}
	if(is.na(prereturns[1])){
		prereturns = tail(data, mar)
	}
	# outpus matrices
	residSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	sigmaSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	seriesSim = matrix(0, ncol = m.sim, nrow = n.sim)
	condmSim = probSim  =  vector(mode="list", length = m.sim)
	constm = matrix(0, ncol = modelinc[46], nrow = n)
	# 2 cases for the star model:
	# case 1: probability determined by external regressors in which case
	# the simulation is quite fast
	# case 2: probability determined by lag(fun(y)) in which case much slower...
	# ...need to consider re-writing probability dynamics in C and also calling an
	# user defined functions from C.
	arglist = list()
	arglist$ipars = ipars
	arglist$idx = idx
	arglist$pmu = pmu
	arglist$model = model
	arglist$probs = probs
	for(i in 1:m.sim){
		x = c(prereturns, rep(0, n))	
		
		for(j in 1:modelinc[46]){
			if(modelinc[c(1,5,9,13)[j]]>0) constm[,j] = as.numeric(ipars[idx[paste("s",j,".phi0",sep=""),1],1])
		}
		if(modelinc[3]>0){
			for(j in 1:modelinc[46]) constm[,j] = constm[,j] + matrix( ipars[idx[paste("s",j,".xi",sep=""),1]:idx[paste("s",j,".xi",sep=""),2], 1], ncol = modelinc[3] ) %*% t( matrix( mexsim[[i]], ncol = modelinc[3] ) )			
		}
		constm = rbind(matrix(0, ncol=modelinc[46], nrow=mar), constm)
		if(modelinc[49]==2){
			arglist$XL = ssim[[i]]
			# sim starts at (m=mar) and end at n+mar
			# get probabilities
			ptmp = switch(modelinc[46], dstar1sim(arglist), dstar2sim(arglist), dstar3sim(arglist), dstar4sim(arglist))
			sigsim = switch(modelinc[46],
					tail(ptmp$probs[,1]*sig1 + ptmp$probs[,2]*sig2, n),
					tail(ptmp$probs[,1]*sig1 + ptmp$probs[,2]*sig2, n),
					tail(ptmp$probs[,1]*sig1 + ptmp$probs[,2]*sig2 + ptmp$probs[,3]*sig3, n),
					tail(ptmp$probs[,1]*sig1 + ptmp$probs[,2]*sig2 + ptmp$probs[,3]*sig3 + ptmp$probs[,4]*sig4, n))
			
			simres = c(preres, zresidSim[,i]*sigsim)
			
			rsim = switch(modelinc[46],
					try(.C("starxsim1", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
									x = as.double(x), s = double(modelinc[46]*(n+mar)), res = as.double(simres), prob = as.double(ptmp$probs), 
									constm = as.double(constm), m = as.integer(mar), T = as.integer(mar+n), PACKAGE="twinkle"), silent = TRUE),
					try(.C("starxsim2", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
									x = as.double(x), s = double(modelinc[46]*(n+mar)), res = as.double(simres), prob = as.double(ptmp$probs), 
									constm = as.double(constm), m = as.integer(mar), T = as.integer(mar+n), PACKAGE="twinkle"), silent = TRUE),
					try(.C("starxsim3", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
									x = as.double(x), s = double(modelinc[46]*(n+mar)), res = as.double(simres), prob = as.double(ptmp$probs), 
									constm = as.double(constm), m = as.integer(mar), T = as.integer(mar+n), PACKAGE="twinkle"), silent = TRUE),
					try(.C("starxsim4", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
									x = as.double(x), s = double(modelinc[46]*(n+mar)), res = as.double(simres), prob = as.double(ptmp$probs), 
									constm = as.double(constm), m = as.integer(mar), T = as.integer(mar+n), PACKAGE="twinkle"), silent = TRUE)
			)			
			seriesSim[,i] = tail(rsim$x, n.sim)
			residSim[,i] = tail(zresidSim[,i]*sigsim,n.sim)
			sigmaSim[,i] = tail(sigsim,n.sim)
			probSim[[i]]  = tail(ptmp$probs, n.sim)
			condmSim[[i]]  = tail(matrix(rsim$s, ncol = modelinc[46]), n.sim)
			colnames(condmSim[[i]]) = colnames(probSim[[i]]) = paste("state[",1:modelinc[46],"]",sep="")
			rownames(condmSim[[i]]) = rownames(probSim[[i]]) = paste("T+",1:n.sim,sep="")			
		} else{
			# NOTE: would probably benefit from using .Call and performing the whole process in C++ with a call to R for "yfun"...
			sigsim = simres = c(rep(0, mar+n))
			simres[1:mar] = preres
			ndx = length(data)
			y = c(data, rep(0, n))
			psim = matrix(0, ncol = modelinc[46], nrow = mar+n)
			xcondm = matrix(0, ncol = modelinc[46], nrow = mar+n)
			for(j in 1:n){
				if(modelinc[48]==1) ytmp = c(yfun(as.numeric(y[1:(ndx+j-1)])), NA) else ytmp = c(as.numeric(y[1:(ndx+j-1)]),NA)		
				XL = matrix(NA, ncol = length(ylags), nrow = ndx+j)
				nXL = nrow(XL)
				for(k in 1:ncol(XL)){
					XL[(ylags[k]+1):nXL,k] = ytmp[1:(nXL-ylags[k])]
				}
				XL[is.na(XL)]=0
				arglist$XL = XL
				ptmp = switch(modelinc[46], dstar1sim(arglist), dstar2sim(arglist), dstar3sim(arglist), dstar4sim(arglist))
				psim[mar+j,] = tail(ptmp$probs, 1)
				
				msimx = switch(modelinc[46],
						psim[mar+j,1]*sig1+psim[mar+j,2]*sig2,
						psim[mar+j,1]*sig1+psim[mar+j,2]*sig2,
						psim[mar+j,1]*sig1+psim[mar+j,2]*sig2+psim[mar+j,3]*sig3,
						psim[mar+j,1]*sig1+psim[mar+j,2]*sig2+psim[mar+j,3]*sig3+psim[mar+j,4]*sig4)
				
				sigsim[mar+j] = msimx
				simres[mar+j] = zresidSim[j,i]*sigsim[mar+j]
				rsim = switch(modelinc[46],
						try(.C("starxsim1", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
										x = as.double(x), s = double(modelinc[46]*(mar+j)), res = as.double(simres[1:(mar+j)]),
										prob = as.double(psim[1:(mar+j),,drop=FALSE]), constm = as.double(constm[1:(mar+j),,drop=FALSE]),
										m = as.integer(mar+j-1), T = as.integer(mar+j), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim2",  model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
										x = as.double(x), s = double(modelinc[46]*(mar+j)), res = as.double(simres[1:(mar+j)]),
										prob = as.double(psim[1:(mar+j),,drop=FALSE]), constm = as.double(constm[1:(mar+j),,drop=FALSE]),
										m = as.integer(mar+j-1), T = as.integer(mar+j), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim3", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
										x = as.double(x), s = double(modelinc[46]*(mar+j)), res = as.double(simres[1:(mar+j)]),
										prob = as.double(psim[1:(mar+j),,drop=FALSE]), constm = as.double(constm[1:(mar+j),,drop=FALSE]),
										m = as.integer(mar+j-1), T = as.integer(mar+j), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim4",  model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
										x = as.double(x), s = double(modelinc[46]*(mar+j)), res = as.double(simres[1:(mar+j)]),
										prob = as.double(psim[1:(mar+j),,drop=FALSE]), constm = as.double(constm[1:(mar+j),,drop=FALSE]),
										m = as.integer(mar+j-1), T = as.integer(mar+j), PACKAGE="twinkle"), silent = TRUE)
				)
				x[mar+j] = rsim$x[mar+j]
				y[ndx+j] = rsim$x[mar+j]
				xcondm[mar+j,] = tail(matrix(rsim$s, ncol = modelinc[46]),1)
			}
			seriesSim[,i] = tail(x, n,sim)
			residSim[,i] = tail(simres, n.sim)
			sigmaSim[,i] = tail(sigsim, n.sim)
			probSim[[i]]  = tail(ptmp$probs, n.sim)
			condmSim[[i]] = tail(xcondm, n.sim)
			colnames(probSim[[i]]) = paste("state[",1:modelinc[46],"]",sep="")
			rownames(probSim[[i]]) = paste("T+",1:n.sim,sep="")
		}
	}
	rownames(residSim) = rownames(seriesSim) = paste("T+",1:n.sim,sep="")
	# check:
	# u=3
	# xtmp = (ipars["s1.phi0",1] + ipars["s1.phi1",1]*x[u-1]+ipars["s1.phi2",1]*x[u-2])*psim[u,1]+ (ipars["s2.phi0",1] + ipars["s2.phi1",1]*x[u-1]+ipars["s2.phi2",1]*x[u-2])*(psim[u,2]) + simres[u]
	sim = list(seriesSim = seriesSim, residSim = residSim, probSim = probSim, sigmaSim = sigmaSim)
	sol = new("STARsim",
			simulation = sim,
			model = model,
			seed = as.integer(sseed))
	return(sol)
}

.starsim.static = function(fit, n.sim = 1000, n.start = 0, m.sim = 1, prereturns = NA, preresiduals = NA, rseed = NA, 
		custom.dist = list(name = NA, distfit = NA), xregsim = NULL, ssim = NULL, probsim = NULL)
{
	# some checks
	if(is.na(rseed[1])){
		sseed = as.integer(runif(1,0,as.integer(Sys.time())))
	} else{
		if(length(rseed) != m.sim) sseed = as.integer(rseed[1]) else sseed = rseed[1:m.sim]
	}
	# Enlarge Series:
	# need to allow for arfima case:
	n = n.sim + n.start
	data = fit@model$modeldata$data
	N = length(as.numeric(data))
	data = data[1:(N - fit@model$n.start)]
	N = length(as.numeric(data))
	resids = fit@fit$residuals	
	sig = coef(fit)["sigma"]
	pmu = fit@fit$pmu
	probs = fit@fit$probability
	condm = fit@fit$condm
	model = fit@model
	modelinc = model$modelinc
	yfun = model$modeldata$fun
	# max of AR terms and STAR dynamics in the case that y is used
	mar = max(c(modelinc[c(2,6,10,14,4,8,12,16,17)], ifelse(modelinc[49]==1, max(model$modeldata$ylags), 0)))
	idx = model$pidx
	ipars = fit@fit$ipars
	ylags = model$modeldata$ylags
	gfun = switch(modelinc[46],"starxsim1", "starxsim2","starxsim3","starxsim4")
	# check if necessary the external regressor forecasts provided first
	xreg = .simregressors(model, xregsim, NULL, ssim, N, n, m.sim)	
	mexsim = xreg$mexsimlist
	ssim = xreg$ssimlist
	# Random Samples from the Distribution
	if(length(sseed) == 1){
		set.seed(sseed)
		residSim = .custzdist(custom.dist = custom.dist, distribution = model$modeldesc$distribution, n = n, m.sim = m.sim, 
				mu = 0, sigma = sig, lambda = ipars[idx["ghlambda",1], 1], 
				skew = ipars[idx["skew",1], 1], shape = ipars[idx["shape",1],1])
	} else{
		residSim = matrix(0, ncol = m.sim, nrow = n)
		for(i in 1:m.sim){
			set.seed(sseed[i])
			residSim[,i] = .custzdist(custom.dist = custom.dist, distribution = model$modeldesc$distribution, n = n, m.sim = 1, 
					mu = 0, sigma = sig, lambda = ipars[idx["ghlambda",1], 1], skew = ipars[idx["skew",1], 1], 
					shape = ipars[idx["shape",1],1])
		}
	}
	# create the presample information
	if(!is.na(prereturns[1])){
		prereturns = as.vector(prereturns)
		if(length(prereturns)<mar) stop(paste("\nstarsim-->error: prereturns must be of length ", mar, sep=""))
	}
	if(!is.na(preresiduals[1])){
		preresiduals = as.vector(preresiduals)
		if(length(preresiduals)<mar) stop(paste("\nstarsim-->error: preresiduals must be of length ", mar, sep=""))
		preres = preresiduals
	} else{
		preres = as.numeric(tail(residuals(fit), mar))
	}
	if(is.na(prereturns[1])){
		prereturns = tail(data, mar)
	}
	# input vectors/matrices
	# outpus matrices
	seriesSim = matrix(0, ncol = m.sim, nrow = n.sim)
	condmSim = probSim  =  vector(mode="list", length = m.sim)
	constm = matrix(0, ncol = modelinc[46], nrow = n)
	# 2 cases for the star model:
	# case 1: probability determined by external regressors in which case
	# the simulation is quite fast
	# case 2: probability determined by lag(fun(y)) in which case much slower...
	# ...need to consider re-writing probability dynamics in C and also calling an
	# user defined functions from C.
	arglist = list()
	arglist$ipars = ipars
	arglist$idx = idx
	arglist$pmu = pmu
	arglist$model = model
	arglist$probs = probs
	for(i in 1:m.sim){
		x = c(prereturns, rep(0, n))	
		
		simres = c(preres, residSim[,i])
		# ToDo: change to accomodate modelinc[20]
		for(j in 1:modelinc[46]){
			if(modelinc[c(1,5,9,13)[j]]>0) constm[,j] = as.numeric(ipars[idx[paste("s",j,".phi0",sep=""),1],1])
		}
		if(modelinc[3]>0){
			for(j in 1:modelinc[46]) constm[,j] = constm[,j] + matrix( ipars[idx[paste("s",j,".xi",sep=""),1]:idx[paste("s",j,".xi",sep=""),2], 1], ncol = modelinc[3] ) %*% t( matrix( mexsim[[i]], ncol = modelinc[3] ) )			
		}
		constm = rbind(matrix(0, ncol=modelinc[46], nrow=mar), constm)
		if(modelinc[49]==2){
			arglist$XL = ssim[[i]]
			# sim starts at (m=mar) and end at n+mar
			# get probabilities
			ptmp = switch(modelinc[46], dstar1sim(arglist), dstar2sim(arglist), dstar3sim(arglist), dstar4sim(arglist))
			
			rsim = switch(modelinc[46],
					try(.C("starxsim1", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
									x = as.double(x), s = double(modelinc[46]*(n+mar)), res = as.double(simres),
									prob = as.double(ptmp$probs), constm = as.double(constm), m = as.integer(mar),
									T = as.integer(mar+n), PACKAGE="twinkle"), silent = TRUE),
					try(.C("starxsim2", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
									x = as.double(x), s = double(modelinc[46]*(n+mar)), res = as.double(simres),
									prob = as.double(ptmp$probs), constm = as.double(constm), m = as.integer(mar),
									T = as.integer(mar+n), PACKAGE="twinkle"), silent = TRUE),
					try(.C("starxsim3",model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
									x = as.double(x), s = double(modelinc[46]*(n+mar)), res = as.double(simres),
									prob = as.double(ptmp$probs), constm = as.double(constm), m = as.integer(mar),
									T = as.integer(mar+n), PACKAGE="twinkle"), silent = TRUE),
					try(.C("starxsim4", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
									x = as.double(x), s = double(modelinc[46]*(n+mar)), res = as.double(simres),
									prob = as.double(ptmp$probs), constm = as.double(constm), m = as.integer(mar),
									T = as.integer(mar+n), PACKAGE="twinkle"), silent = TRUE)
			)
			seriesSim[,i] = tail(rsim$x, n.sim)
			probSim[[i]]  = tail(ptmp$probs, n.sim)
			condmSim[[i]]  = tail(matrix(rsim$s, ncol = modelinc[46]), n.sim)
			colnames(condmSim[[i]]) = colnames(probSim[[i]]) = paste("state[",1:modelinc[46],"]",sep="")
			rownames(condmSim[[i]]) = rownames(probSim[[i]]) = paste("T+",1:n.sim,sep="")			
		} else{
			# NOTE: would probably benefit from using .Call and performing the whole process in C++ with a call to R for "yfun"...
			ndx = length(data)
			# y is padded with pre-returns
			y = c(data, rep(0, n))
			y[(N-mar+1):N] = prereturns
			psim = matrix(0, ncol = modelinc[46], nrow = mar+n)
			xcondm = matrix(0, ncol = modelinc[46], nrow = mar+n)
			for(j in 1:n){
				if(modelinc[48]==1) ytmp = c(yfun(as.numeric(y[1:(ndx+j-1)])), NA) else ytmp = c(as.numeric(y[1:(ndx+j-1)]),NA)		
				XL = matrix(NA, ncol = length(ylags), nrow = ndx+j)
				nXL = nrow(XL)
				for(k in 1:ncol(XL)){
					XL[(ylags[k]+1):nXL,k] = ytmp[1:(nXL-ylags[k])]
				}
				XL[is.na(XL)]=0
				arglist$XL = XL
				ptmp = switch(modelinc[46], dstar1sim(arglist), dstar2sim(arglist), dstar3sim(arglist), dstar4sim(arglist))
				psim[mar+j,] = tail(ptmp$probs, 1)
				rsim = switch(modelinc[46],
						try(.C("starxsim1", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
										x = as.double(x), s = double(modelinc[46]*(mar+j)),res = as.double(simres[1:(mar+j)]),
										prob = as.double(psim[1:(mar+j),,drop=FALSE]), constm = as.double(constm[1:(mar+j),,drop=FALSE]),
										m = as.integer(mar+j-1), T = as.integer(mar+j), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim2", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
										x = as.double(x), s = double(modelinc[46]*(mar+j)),res = as.double(simres[1:(mar+j)]),
										prob = as.double(psim[1:(mar+j),,drop=FALSE]), constm = as.double(constm[1:(mar+j),,drop=FALSE]),
										m = as.integer(mar+j-1), T = as.integer(mar+j), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim3", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
										x = as.double(x), s = double(modelinc[46]*(mar+j)),res = as.double(simres[1:(mar+j)]),
										prob = as.double(psim[1:(mar+j),,drop=FALSE]), constm = as.double(constm[1:(mar+j),,drop=FALSE]),
										m = as.integer(mar+j-1), T = as.integer(mar+j), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim4", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
										x = as.double(x), s = double(modelinc[46]*(mar+j)),res = as.double(simres[1:(mar+j)]),
										prob = as.double(psim[1:(mar+j),,drop=FALSE]), constm = as.double(constm[1:(mar+j),,drop=FALSE]),
										m = as.integer(mar+j-1), T = as.integer(mar+j), PACKAGE="twinkle"), silent = TRUE)
				)
				x[mar+j] = rsim$x[mar+j]
				y[ndx+j] = rsim$x[mar+j]
				xcondm[mar+j,] = tail(matrix(rsim$s, ncol = modelinc[46]),1)
			}
			seriesSim[,i] = tail(x, n.sim)
			probSim[[i]]  = tail(ptmp$probs, n.sim)
			condmSim[[i]] = tail(xcondm, n.sim)
			colnames(probSim[[i]]) = paste("state[",1:modelinc[46],"]",sep="")
			rownames(probSim[[i]]) = paste("T+",1:n.sim,sep="")
		}
	}
	residSim = tail(residSim, n.sim)
	rownames(residSim) = rownames(seriesSim) = paste("T+",1:n.sim,sep="")
	# check:
	# u=3
	# xtmp = (ipars["s1.phi0",1] + ipars["s1.phi1",1]*x[u-1]+ipars["s1.phi2",1]*x[u-2])*psim[u,1]+ (ipars["s2.phi0",1] + ipars["s2.phi1",1]*x[u-1]+ipars["s2.phi2",1]*x[u-2])*(psim[u,2]) + simres[u]
	sim = list(seriesSim = seriesSim, residSim = residSim, probSim = probSim, condmSim = condmSim)
	sol = new("STARsim",
			simulation = sim,
			model = model,
			seed = as.integer(sseed))
	return(sol)
}


.starpath.static = function(spec, n.sim = 1000, n.start = 0, m.sim = 1, prereturns = NA,
		preresiduals = NA, rseed = NA, custom.dist = list(name = NA, distfit = NA), 
		xregsim = NULL, ssim = NULL, probsim = NULL)
{
	# some checks
	if(is.na(rseed[1])){
		sseed = as.integer(runif(1,0,as.integer(Sys.time())))
	} else{
		if(length(rseed) != m.sim) sseed = as.integer(rseed[1]) else sseed = rseed[1:m.sim]
	}
	# Enlarge Series:
	# need to allow for arfima case:
	model = spec@model
	ipars = model$pars
	pars = ipars[ipars[,"Fixed"]==1,1]
	parnames = names(pars)
	modelnames = .checkallfixed(spec)
	if(is.na(all(match(modelnames, parnames), 1:length(modelnames)))) {
		cat("\nstarpath-->error: parameters names do not match specification\n")
		cat("Expected Parameters are: ")
		cat(paste(modelnames))
		cat("\n")
		stop("Exiting", call. = FALSE)
	}
	# once more into the spec
	setfixed(spec)<-as.list(pars)
	model = spec@model
	ipars = model$pars
	idx = model$pidx
	modelinc = model$modelinc
	# Enlarge Series:
	n = n.sim + n.start
	m = spec@model$maxOrder
	N = 0
	xreg = .simregressorspath(model, xregsim = xregsim, vregsim = NULL, ssim = ssim, n, m.sim)
	mexsim = xreg$mexsimlist
	ssim = xreg$ssimlist
	distribution = model$modeldesc$distribution	
	sig = ipars["sigma",1]
	yfun = model$modeldata$fun
	# max of AR terms and STAR dynamics in the case that y is used
	mar = max(c(modelinc[c(2,6,10,14,4,8,12,16,17)], ifelse(modelinc[49]==1, max(model$modeldata$ylags), 0)))
	ylags = model$modeldata$ylags
	gfun = switch(modelinc[46],"starxsim1", "starxsim2","starxsim3","starxsim4")
	# check if necessary the external regressor forecasts provided first
	# Random Samples from the Distribution
	if(length(sseed) == 1){
		set.seed(sseed)
		residSim = .custzdist(custom.dist = custom.dist, distribution = model$modeldesc$distribution, n = n, m.sim = m.sim, 
				mu = 0, sigma = 1, lambda = ipars[idx["ghlambda",1], 1], 
				skew = ipars[idx["skew",1], 1], shape = ipars[idx["shape",1],1])*sig
	} else{
		residSim = matrix(0, ncol = m.sim, nrow = n)
		for(i in 1:m.sim){
			set.seed(sseed[i])
			residSim[,i] = .custzdist(custom.dist = custom.dist, distribution = model$modeldesc$distribution, n = n, m.sim = 1, 
					mu = 0, sigma = 1, lambda = ipars[idx["ghlambda",1], 1], skew = ipars[idx["skew",1], 1], 
					shape = ipars[idx["shape",1],1])*sig
		}
	}
	# create the presample information
	if(!is.na(prereturns[1])){
		prereturns = as.vector(prereturns)
		if(length(prereturns)<mar) stop(paste("\nstarpath-->error: prereturns must be of length ", mar, sep=""))
		prereturnsx = tail(prereturns, mar)
		prereturnsy = prereturns
		yfn = length(prereturnsy)
	} else{
		stop("\nstarpath-->error: prereturns cannot be NA when using starpath method.")
	}
	if(!is.na(preresiduals[1])){
		preresiduals = as.vector(preresiduals)
		if(length(preresiduals)<mar) stop(paste("\nstarpath-->error: preresiduals must be of length ", mar, sep=""))
		preres = tail(preresiduals, mar)
	} else{
		preres = rep(0, mar)
	}
	# input vectors/matrices
	# outpus matrices
	seriesSim = matrix(0, ncol = m.sim, nrow = n.sim)
	condmSim = probSim  =  vector(mode="list", length = m.sim)
	constm = matrix(0, ncol = modelinc[46], nrow = n)
	# 2 cases for the star model:
	# case 1: probability determined by external regressors in which case
	# the simulation is quite fast
	# case 2: probability determined by lag(fun(y)) in which case much slower...
	# ...need to consider re-writing probability dynamics in C and also calling an
	# user defined functions from C.
	arglist = list()
	arglist$ipars = ipars
	arglist$idx = idx
	arglist$model = model
	arglist$probs = matrix(0, ncol=modelinc[46], nrow=n)
	
	for(i in 1:m.sim){
		x = c(prereturnsx, rep(0, n))	
		
		simres = c(preres, residSim[,i])
		for(j in 1:modelinc[46]){
			if(modelinc[c(1,5,9,13)[j]]>0) constm[,j] = as.numeric(ipars[idx[paste("s",j,".phi0",sep=""),1],1])
		}
		if(modelinc[3]>0){
			for(j in 1:modelinc[46]) constm[,j] = constm[,j] + matrix( ipars[idx[paste("s",j,".xi",sep=""),1]:idx[paste("s",j,".xi",sep=""),2], 1], ncol = modelinc[3] ) %*% t( matrix( mexsim[[i]], ncol = modelinc[3] ) )			
		}
		constm = rbind(matrix(0, ncol=modelinc[46], nrow=mar), constm)
		if(modelinc[49]==2){
			arglist$XL = ssim[[i]]
			# sim starts at (m=mar) and end at n+mar
			# get probabilities
			ptmp = switch(modelinc[46], dstar1path(arglist), dstar2path(arglist), dstar3path(arglist), dstar4path(arglist))
			rsim = switch(modelinc[46],
					try(.C("starxsim1", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
									x = as.double(x), s = double(modelinc[46]*(n+mar)), res = as.double(simres), prob = as.double(ptmp$probs), 
									constm = as.double(constm), m = as.integer(mar), T = as.integer(mar+n), PACKAGE="twinkle"), silent = TRUE),
					try(.C("starxsim2", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
									x = as.double(x), s = double(modelinc[46]*(n+mar)), res = as.double(simres), prob = as.double(ptmp$probs), 
									constm = as.double(constm), m = as.integer(mar), T = as.integer(mar+n), PACKAGE="twinkle"), silent = TRUE),
					try(.C("starxsim3", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
									x = as.double(x), s = double(modelinc[46]*(n+mar)), res = as.double(simres), prob = as.double(ptmp$probs), 
									constm = as.double(constm), m = as.integer(mar), T = as.integer(mar+n), PACKAGE="twinkle"), silent = TRUE),
					try(.C("starxsim4", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
									x = as.double(x), s = double(modelinc[46]*(n+mar)), res = as.double(simres), prob = as.double(ptmp$probs), 
									constm = as.double(constm), m = as.integer(mar), T = as.integer(mar+n), PACKAGE="twinkle"), silent = TRUE)
			)
			seriesSim[,i] = tail(rsim$x, n.sim)
			probSim[[i]]  = tail(ptmp$probs, n.sim)
			condmSim[[i]]  = tail(matrix(rsim$s, ncol = modelinc[46]), n.sim)
			colnames(condmSim[[i]]) = colnames(probSim[[i]]) = paste("state[",1:modelinc[46],"]",sep="")
			rownames(condmSim[[i]]) = rownames(probSim[[i]]) = paste("T+",1:n.sim)	
		} else{
			# NOTE: would probably benefit from using .Call and performing the whole process in C++ with a call to R for "yfun"...
			# y is padded with pre-returns
			x = c(prereturnsx, rep(0, n))
			y = c(prereturnsx, rep(0, n))
			# NOTE: sometimes the user may need to pass a vector which is longer since yfun may have
			# a lag structure (see base example with vDF dataset).
			yf = c(prereturnsy, rep(0, n))
			psim = matrix(0, ncol = modelinc[46], nrow = mar+n)
			xcondm = matrix(0, ncol = modelinc[46], nrow = mar+n)
			for(j in 1:n){
				if(modelinc[48]==1){
					ytmp = c(yfun(as.numeric(yf[1:(yfn+j-1)])), NA)
					ytmp = tail(ytmp, mar+j)
				} else{
					ytmp = c(as.numeric(y[1:(mar+j-1)]),NA)	
				}
				XL = matrix(NA, ncol = length(ylags), nrow = mar+j)
				nXL = nrow(XL)
				for(k in 1:ncol(XL)){
					XL[(ylags[k]+1):nXL,k] = ytmp[1:(nXL-ylags[k])]
				}
				XL[is.na(XL)]=0
				arglist$XL = XL
				ptmp = switch(modelinc[46], dstar1path(arglist), dstar2path(arglist), dstar3path(arglist), dstar4path(arglist))
				psim[mar+j,] = tail(ptmp$probs, 1)
				rsim = switch(modelinc[46],
						try(.C("starxsim1", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
										x = as.double(x), s = double(modelinc[46]*(mar+j)),res = as.double(simres[1:(mar+j)]),
										prob = as.double(psim[1:(mar+j),,drop=FALSE]), constm = as.double(constm[1:(mar+j),,drop=FALSE]),
										m = as.integer(mar+j-1), T = as.integer(mar+j), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim2", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
										x = as.double(x), s = double(modelinc[46]*(mar+j)),res = as.double(simres[1:(mar+j)]),
										prob = as.double(psim[1:(mar+j),,drop=FALSE]), constm = as.double(constm[1:(mar+j),,drop=FALSE]),
										m = as.integer(mar+j-1), T = as.integer(mar+j), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim3", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
										x = as.double(x), s = double(modelinc[46]*(mar+j)),res = as.double(simres[1:(mar+j)]),
										prob = as.double(psim[1:(mar+j),,drop=FALSE]), constm = as.double(constm[1:(mar+j),,drop=FALSE]),
										m = as.integer(mar+j-1), T = as.integer(mar+j), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim4", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
										x = as.double(x), s = double(modelinc[46]*(mar+j)),res = as.double(simres[1:(mar+j)]),
										prob = as.double(psim[1:(mar+j),,drop=FALSE]), constm = as.double(constm[1:(mar+j),,drop=FALSE]),
										m = as.integer(mar+j-1), T = as.integer(mar+j), PACKAGE="twinkle"), silent = TRUE)
				)
				x[mar+j] = rsim$x[mar+j]
				y[mar+j] = rsim$x[mar+j]
				yf[yfn+j] = y[mar+j]
				xcondm[mar+j,] = tail(matrix(rsim$s, ncol = modelinc[46]),1)
			}
			seriesSim[,i] = tail(x, n.sim)
			probSim[[i]]  = tail(ptmp$probs, n.sim)
			condmSim[[i]] = tail(xcondm, n.sim)
			colnames(probSim[[i]]) = paste("state[",1:modelinc[46],"]",sep="")
			rownames(probSim[[i]]) = paste("T+",1:n.sim)
		}
	}
	residSim = tail(residSim, n.sim)
	rownames(residSim) = rownames(seriesSim) = paste("T+",1:n.sim)
	# check:
	# u=3
	# xtmp = (ipars["s1.phi0",1] + ipars["s1.phi1",1]*x[u-1]+ipars["s1.phi2",1]*x[u-2])*psim[u,1]+ (ipars["s2.phi0",1] + ipars["s2.phi1",1]*x[u-1]+ipars["s2.phi2",1]*x[u-2])*(psim[u,2]) + simres[u]
	sim = list(seriesSim = seriesSim, residSim = residSim, probSim = probSim, condmSim = condmSim)
	sol = new("STARpath",
			path = sim,
			model = model,
			seed = as.integer(sseed))
	return(sol)
}


.starpath.dynamic = function(spec, n.sim = 1000, n.start = 0, m.sim = 1, presigma = NA, 
		prereturns = NA, preresiduals = NA, rseed = NA, 
		custom.dist = list(name = NA, distfit = NA), xregsim = NULL, 
		vregsim = NULL, ssim = NULL, probsim = NULL)
{
	# some checks
	# some checks
	if(is.na(rseed[1])){
		sseed = as.integer(runif(1,0,as.integer(Sys.time())))
	} else{
		if(length(rseed) != m.sim) sseed = as.integer(rseed[1]) else sseed = rseed[1:m.sim]
	}
	model = spec@model
	ipars = model$pars
	pars = ipars[ipars[,"Fixed"]==1,1]
	parnames = names(pars)
	modelnames = .checkallfixed(spec)
	if(is.na(all(match(modelnames, parnames), 1:length(modelnames)))) {
		cat("\nstarpath-->error: parameters names do not match specification\n")
		cat("Expected Parameters are: ")
		cat(paste(modelnames))
		cat("\n")
		stop("Exiting", call. = FALSE)
	}
	# once more into the spec
	setfixed(spec)<-as.list(pars)
	model = spec@model
	ipars = model$pars
	idx = model$pidx
	modelinc = model$modelinc
	# Enlarge Series:
	n = n.sim + n.start
	N = 0
	xreg = .simregressorspath(model, xregsim = xregsim, vregsim = vregsim, ssim = ssim, n, m.sim)
	mexsim = xreg$mexsimlist
	ssim = xreg$ssimlist
	vexsim = xreg$vexsimlist
	distribution = model$modeldesc$distribution	
	sig = ipars["sigma",1]
	yfun = model$modeldata$fun
	mar = max(c(modelinc[c(2,6,10,14,4,8,12,16,17)], ifelse(modelinc[49]==1, max(model$modeldata$ylags), 0)))
	m = max(c(modelinc[32:33], mar))
	ylags = model$modeldata$ylags
	gfun = switch(modelinc[46],"starxsim1", "starxsim2","starxsim3","starxsim4")
	if(length(sseed) == 1){
		set.seed(sseed)
		z = .custzdist(custom.dist = custom.dist, distribution = model$modeldesc$distribution, n = n, m.sim = m.sim, 
				mu = 0, sigma = 1, lambda = ipars[idx["ghlambda",1], 1], 
				skew = ipars[idx["skew",1], 1], shape = ipars[idx["shape",1],1])
	} else{
		z = matrix(0, ncol = m.sim, nrow = n)
		for(i in 1:m.sim){
			set.seed(sseed[i])
			z[,i] = .custzdist(custom.dist = custom.dist, distribution = model$modeldesc$distribution, n = n, m.sim = 1, 
					mu = 0, sigma = 1, lambda = ipars[idx["ghlambda",1], 1], skew = ipars[idx["skew",1], 1], 
					shape = ipars[idx["shape",1],1])
		}
	}
	z = rbind(matrix(rep(0, m), nrow = m, ncol = m.sim), z)
	
	# create the presample information
	if(!is.na(presigma[1])){
		presigma = as.vector(presigma)
		if(length(presigma)<m) stop(paste("\nstarpath-->error: presigma must be of length ", m, sep=""))
		presigma = tail(presigma,m)
	} else{
		stop("\nstarpath-->error: presigma cannot be NA when using starpath method.")
	}
	if(!is.na(prereturns[1])){
		prereturns = as.vector(prereturns)
		if(length(prereturns)<mar) stop(paste("\nstarpath-->error: prereturns must be of length ", mar, sep=""))
		prereturnsx = tail(prereturns, mar)
		prereturnsy = prereturns
		yfn = length(prereturnsy)
	} else{
		stop("\nstarpath-->error: prereturns cannot be NA when using starpath method.")	
	}
	if(!is.na(preresiduals[1])){
		preresiduals = as.vector(preresiduals)
		if(length(preresiduals)<m) stop(paste("\nstarpath-->error: preresiduals must be of length ", m, sep=""))
		preres = tail(as.vector(preresiduals),m)
		z[1:m,] = preres/presigma
	} else{
		warning("\nstarpath-->warning: preresiduals not provided...setting to zero.")
		preres = rep(0, m)
	}
	
	# input vectors/matrices
	# outpus matrices
	sigmaSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	seriesSim = matrix(0, ncol = m.sim, nrow = n.sim)
	residSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	condmSim = probSim  =  vector(mode="list", length = m.sim)
	z[is.na(z) | is.nan(z) | !is.finite(z)] = 0
	constm = matrix(0, ncol = modelinc[46], nrow = n)
	# 2 cases for the star model:
	# case 1: probability determined by external regressors in which case
	# the simulation is quite fast
	# case 2: probability determined by lag(fun(y)) in which case much slower...
	# ...need to consider re-writing probability dynamics in C and also calling an
	# user defined functions from C.
	arglist = list()
	arglist$ipars = ipars
	arglist$idx = idx
	arglist$model = model
	arglist$probs = matrix(0, ncol=modelinc[38], nrow=n)
	for(i in 1:m.sim){
		h = c(presigma^2, rep(0, n))
		x = c(prereturnsx, rep(0, n))
		
		if(modelinc[50]==2){
			ngrd = which(preres<0)
			tmpr = rep(0, length(preres))
			tmpr[ngrd] = preres[ngrd] * preres[ngrd]
			nres = c(tmpr,   rep(0, n))
		} else{
			nres = 0
		}
		if(modelinc[50]==3){
			kappa = egarchKappa(ipars[idx["ghlambda",1],1], ipars[idx["shape",1],1], ipars[idx["skew",1],1], model$modeldesc$distribution)	
		} else{
			kappa = 1
		}
		res = c(preres, rep(0, n))

		ans1 = try(.C("garchsimC", model = as.integer(modelinc[1:60]), pars = as.double(ipars[,1]), 
						idx = as.integer(idx[,1]-1), h = as.double(h), z = as.double(z[,i]), 
						res = as.double(res), e = as.double(res*res), nres = as.double(nres),
						meanz = as.double(kappa), vexdata = as.double(vexsim[[i]]), 
						T = as.integer(n+m), m = as.integer(m), PACKAGE = "twinkle"), silent = TRUE)
		if(inherits(ans1, "try-error")) stop("\nstarpath-->error: error in calling C function....\n")
		sigmaSim[,i] = ans1$h[(n.start + m + 1):(n+m)]^(1/2)
		residSim[,i] = ans1$res[(n.start + m + 1):(n+m)]
		simres = tail(ans1$res,n+mar)
		# ToDo: change to accomodate modelinc[20]
		for(j in 1:modelinc[46]){
			if(modelinc[c(1,5,9,13)[j]]>0) constm[,j] = as.numeric(ipars[idx[paste("s",j,".phi0",sep=""),1],1])
		}
		if(modelinc[3]>0){
			for(j in 1:modelinc[46]) constm[,j] = constm[,j] + matrix( ipars[idx[paste("s",j,".xi",sep=""),1]:idx[paste("s",j,".xi",sep=""),2], 1], ncol = modelinc[3] ) %*% t( matrix( mexsim[[i]], ncol = modelinc[3] ) )
		}
		constm = rbind(matrix(0, ncol=modelinc[46], nrow=mar), constm)
		if(modelinc[49]==2){
			arglist$XL = ssim[[i]]
			ptmp = switch(modelinc[46], dstar1path(arglist), dstar2path(arglist), dstar3path(arglist), dstar4path(arglist))
			rsim = switch(modelinc[46],
					try(.C("starxsim1", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
									x = as.double(x), s = double(modelinc[46]*(n+mar)), res = as.double(simres),
									prob = as.double(ptmp$probs), constm = as.double(constm), m = as.integer(mar),
									T = as.integer(mar+n), PACKAGE="twinkle"), silent = TRUE),
					try(.C("starxsim2", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
									x = as.double(x), s = double(modelinc[46]*(n+mar)), res = as.double(simres),
									prob = as.double(ptmp$probs), constm = as.double(constm), m = as.integer(mar),
									T = as.integer(mar+n), PACKAGE="twinkle"), silent = TRUE),
					try(.C("starxsim3", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
									x = as.double(x), s = double(modelinc[46]*(n+mar)), res = as.double(simres),
									prob = as.double(ptmp$probs), constm = as.double(constm), m = as.integer(mar),
									T = as.integer(mar+n), PACKAGE="twinkle"), silent = TRUE),
					try(.C("starxsim4", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
									x = as.double(x), s = double(modelinc[46]*(n+mar)), res = as.double(simres),
									prob = as.double(ptmp$probs), constm = as.double(constm), m = as.integer(mar),
									T = as.integer(mar+n), PACKAGE="twinkle"), silent = TRUE)
			)
			seriesSim[,i] = tail(rsim$x, n.sim)
			probSim[[i]]  = tail(ptmp$probs, n.sim)
			condmSim[[i]]  = tail(matrix(rsim$s, ncol = modelinc[46]), n.sim)
			colnames(condmSim[[i]]) = colnames(probSim[[i]]) = paste("state[",1:modelinc[46],"]",sep="")
			rownames(condmSim[[i]]) = rownames(probSim[[i]]) = paste("T+",1:n.sim,sep="")	
		} else{
			# NOTE: would probably benefit from using .Call and performing the whole process in C++ with a call to R for "yfun"...
			x = c(prereturnsx, rep(0, n))
			y = c(prereturnsx, rep(0, n))
			yf = c(prereturnsy, rep(0, n))
			psim = matrix(0, ncol = modelinc[46], nrow = mar+n)
			xcondm = matrix(0, ncol = modelinc[46], nrow = mar+n)
			for(j in 1:n){
				if(modelinc[48]==1){
					ytmp = c(yfun(as.numeric(yf[1:(yfn+j-1)])), NA)
					ytmp = tail(ytmp, mar+j)
				} else{
					ytmp = c(as.numeric(y[1:(mar+j-1)]),NA)	
				}
				XL = matrix(NA, ncol = length(ylags), nrow = mar+j)
				nXL = nrow(XL)
				for(k in 1:ncol(XL)){
					XL[(ylags[k]+1):nXL,k] = ytmp[1:(nXL-ylags[k])]
				}
				XL[is.na(XL)]=0
				arglist$XL = XL
				ptmp = switch(modelinc[46], dstar1path(arglist), dstar2path(arglist), dstar3path(arglist), dstar4path(arglist))
				psim[mar+j,] = tail(ptmp$probs, 1)
				rsim = switch(modelinc[46],
						try(.C("starxsim1", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
										x = as.double(x), s = double(modelinc[46]*(mar+j)), res = as.double(simres[1:(mar+j)]),
										prob = as.double(psim[1:(mar+j),,drop=FALSE]), constm = as.double(constm[1:(mar+j),,drop=FALSE]),
										m = as.integer(mar+j-1), T = as.integer(mar+j), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim2",  model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
										x = as.double(x), s = double(modelinc[46]*(mar+j)), res = as.double(simres[1:(mar+j)]),
										prob = as.double(psim[1:(mar+j),,drop=FALSE]), constm = as.double(constm[1:(mar+j),,drop=FALSE]),
										m = as.integer(mar+j-1), T = as.integer(mar+j), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim3",  model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
										x = as.double(x), s = double(modelinc[46]*(mar+j)), res = as.double(simres[1:(mar+j)]),
										prob = as.double(psim[1:(mar+j),,drop=FALSE]), constm = as.double(constm[1:(mar+j),,drop=FALSE]),
										m = as.integer(mar+j-1), T = as.integer(mar+j), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim4",  model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
										x = as.double(x), s = double(modelinc[46]*(mar+j)), res = as.double(simres[1:(mar+j)]),
										prob = as.double(psim[1:(mar+j),,drop=FALSE]), constm = as.double(constm[1:(mar+j),,drop=FALSE]),
										m = as.integer(mar+j-1), T = as.integer(mar+j), PACKAGE="twinkle"), silent = TRUE),
				)
				x[mar+j] = rsim$x[mar+j]
				y[mar+j] = rsim$x[mar+j]
				yf[yfn+j] = y[mar+j]
				xcondm[mar+j,] = tail(matrix(rsim$s, ncol = modelinc[46]),1)
			}
			seriesSim[,i] = tail(x, n.sim)
			probSim[[i]]  = tail(ptmp$probs, n.sim)
			condmSim[[i]] = tail(xcondm, n.sim)
			colnames(probSim[[i]]) = paste("state[",1:modelinc[46],"]",sep="")
			rownames(probSim[[i]]) = paste("T+",1:n.sim,sep="")
		}
	}
	residSim = tail(residSim, n.sim)
	rownames(residSim) = rownames(seriesSim) = rownames(sigmaSim) = paste("T+",1:n.sim,sep="")
	# check:
	# u=3
	# xtmp = (ipars["s1.phi0",1] + ipars["s1.phi1",1]*x[u-1]+ipars["s1.phi2",1]*x[u-2])*psim[u,1]+ (ipars["s2.phi0",1] + ipars["s2.phi1",1]*x[u-1]+ipars["s2.phi2",1]*x[u-2])*(psim[u,2]) + simres[u]
	sim = list(sigmaSim = sigmaSim, seriesSim = seriesSim, residSim = residSim, probSim = probSim)
	sol = new("STARpath",
			path = sim,
			model = model,
			seed = as.integer(sseed))
	return(sol)
}

.starpath.mixture = function(spec, n.sim = 1000, n.start = 0, m.sim = 1, prereturns = NA, preresiduals = NA, 
		rseed = NA, custom.dist = list(name = NA, distfit = NA), xregsim = NULL, ssim = NULL, 
		probsim = NULL)
{
	# some checks
	if(is.na(rseed[1])){
		sseed = as.integer(runif(1,0,as.integer(Sys.time())))
	} else{
		if(length(rseed) != m.sim) sseed = as.integer(rseed[1]) else sseed = rseed[1:m.sim]
	}
	# Enlarge Series:
	# need to allow for arfima case:
	model = spec@model
	ipars = model$pars
	pars = ipars[ipars[,"Fixed"]==1,1]
	parnames = names(pars)
	modelnames = .checkallfixed(spec)
	if(is.na(all(match(modelnames, parnames), 1:length(modelnames)))) {
		cat("\nstarpath-->error: parameters names do not match specification\n")
		cat("Expected Parameters are: ")
		cat(paste(modelnames))
		cat("\n")
		stop("Exiting", call. = FALSE)
	}
	# once more into the spec
	setfixed(spec)<-as.list(pars)
	model = spec@model
	ipars = model$pars
	idx = model$pidx
	modelinc = model$modelinc
	# Enlarge Series:
	n = n.sim + n.start
	N = 0
	sig1 = ipars["s1.sigma",1]
	sig2 = ipars["s2.sigma",1]
	if(modelinc[46]==3)	sig3 = ipars["s3.sigma",1] else sig3 = NULL
	if(modelinc[46]==4)	sig4 = ipars["s4.sigma",1] else sig4 = NULL
	
	xreg = .simregressorspath(model, xregsim = xregsim, vregsim = NULL, ssim = ssim, n, m.sim)
	mexsim = xreg$mexsimlist
	ssim = xreg$ssimlist
	distribution = model$modeldesc$distribution	
	yfun = model$modeldata$fun
	# max of AR terms and STAR dynamics in the case that y is used
	mar = max(c(modelinc[c(2,6,10,14,4,8,12,16,17)], ifelse(modelinc[49]==1, max(model$modeldata$ylags), 0)))
	ylags = model$modeldata$ylags
	gfun = switch(modelinc[46],"starxsim1", "starxsim2","starxsim3","starxsim4")
	# Random Samples from the Distribution
	# will postmultiply by mixture sigma once we have probs
	if(length(sseed) == 1){
		set.seed(sseed)
		zresidSim = .custzdist(custom.dist = custom.dist, distribution = model$modeldesc$distribution, n = n, m.sim = m.sim, 
				mu = 0, sigma = 1, lambda = ipars[idx["ghlambda",1], 1], 
				skew = ipars[idx["skew",1], 1], shape = ipars[idx["shape",1],1])
	} else{
		residSim = matrix(0, ncol = m.sim, nrow = n)
		for(i in 1:m.sim){
			set.seed(sseed[i])
			zresidSim[,i] = .custzdist(custom.dist = custom.dist, distribution = model$modeldesc$distribution, n = n, m.sim = 1, 
					mu = 0, sigma = 1, lambda = ipars[idx["ghlambda",1], 1], skew = ipars[idx["skew",1], 1], 
					shape = ipars[idx["shape",1],1])
		}
	}
	# create the presample information
	if(!is.na(prereturns[1])){
		prereturns = as.vector(prereturns)
		if(length(prereturns)<mar) stop(paste("\nstarpath-->error: prereturns must be of length ", mar, sep=""))
		prereturnsx = tail(prereturns, mar)
		prereturnsy = prereturns
		yfn = length(prereturnsy)
	} else{
		stop("\nstarpath-->error: prereturns cannot be NA when using starpath method.")
	}
	if(!is.na(preresiduals[1])){
		preresiduals = as.vector(preresiduals)
		if(length(preresiduals)<mar) stop(paste("\nstarpath-->error: preresiduals must be of length ", mar, sep=""))
		preres = tail(preresiduals, mar)
	} else{
		preres = rep(0, mar)
	}
	# input vectors/matrices
	# outpus matrices
	residSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	sigmaSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	seriesSim = matrix(0, ncol = m.sim, nrow = n.sim)
	condmSim = probSim  =  vector(mode="list", length = m.sim)
	constm = matrix(0, ncol = modelinc[46], nrow = n)
	# 2 cases for the star model:
	# case 1: probability determined by external regressors in which case
	# the simulation is quite fast
	# case 2: probability determined by lag(fun(y)) in which case much slower...
	# ...need to consider re-writing probability dynamics in C and also calling an
	# user defined functions from C.
	arglist = list()
	arglist$ipars = ipars
	arglist$idx = idx
	arglist$model = model
	arglist$probs = matrix(0, ncol=modelinc[46], nrow=n)
	for(i in 1:m.sim){
		x = c(prereturnsx, rep(0, n))
		
		for(j in 1:modelinc[46]){
			if(modelinc[c(1,5,9,13)[j]]>0) constm[,j] = as.numeric(ipars[idx[paste("s",j,".phi0",sep=""),1],1])
		}
		if(modelinc[3]>0){
			for(j in 1:modelinc[46]) constm[,j] = constm[,j] + matrix( ipars[idx[paste("s",j,".xi",sep=""),1]:idx[paste("s",j,".xi",sep=""),2], 1], ncol = modelinc[3] ) %*% t( matrix( mexsim[[i]], ncol = modelinc[3] ) )			
		}
		constm = rbind(matrix(0, ncol=modelinc[46], nrow=mar), constm)
		if(modelinc[49]==2){
			arglist$XL = ssim[[i]]
			ptmp = switch(modelinc[46], dstar1path(arglist), dstar2path(arglist), dstar3path(arglist), dstar4path(arglist))
			# This is a 2 state model
			sigsim = switch(modelinc[46],
					tail(ptmp$probs[,1]*sig1 + ptmp$probs[,2]*sig2, n),
					tail(ptmp$probs[,1]*sig1 + ptmp$probs[,2]*sig2, n),
					tail(ptmp$probs[,1]*sig1 + ptmp$probs[,2]*sig2 + ptmp$probs[,3]*sig3, n),
					tail(ptmp$probs[,1]*sig1 + ptmp$probs[,2]*sig2 + ptmp$probs[,3]*sig3 + ptmp$probs[,4]*sig4, n))
				
			simres = c(preres, zresidSim[,i]*sigsim)
			rsim = switch(modelinc[46],
					try(.C("starxsim1", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
									x = as.double(x), s = double(modelinc[46]*(n+mar)), res = as.double(simres),
									prob = as.double(ptmp$probs), constm = as.double(constm), m = as.integer(mar),
									T = as.integer(mar+n), PACKAGE="twinkle"), silent = TRUE),
					try(.C("starxsim2", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
									x = as.double(x), s = double(modelinc[46]*(n+mar)), res = as.double(simres),
									prob = as.double(ptmp$probs), constm = as.double(constm), m = as.integer(mar),
									T = as.integer(mar+n), PACKAGE="twinkle"), silent = TRUE),
					try(.C("starxsim3", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
									x = as.double(x), s = double(modelinc[46]*(n+mar)), res = as.double(simres),
									prob = as.double(ptmp$probs), constm = as.double(constm), m = as.integer(mar),
									T = as.integer(mar+n), PACKAGE="twinkle"), silent = TRUE),
					try(.C("starxsim4",model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
									x = as.double(x), s = double(modelinc[46]*(n+mar)), res = as.double(simres),
									prob = as.double(ptmp$probs), constm = as.double(constm), m = as.integer(mar),
									T = as.integer(mar+n), PACKAGE="twinkle"), silent = TRUE)
			)
			seriesSim[,i] = tail(rsim$x, n.sim)
			residSim[,i] = tail(zresidSim[,i]*sigsim,n.sim)
			sigmaSim[,i] = tail(sigsim,n.sim)
			probSim[[i]]  = tail(ptmp$probs, n.sim)
			condmSim[[i]]  = tail(matrix(rsim$s, ncol = modelinc[46]), n.sim)
			colnames(condmSim[[i]]) = colnames(probSim[[i]]) = paste("state[",1:modelinc[46],"]",sep="")
			rownames(condmSim[[i]]) = rownames(probSim[[i]]) = paste("T+",1:n.sim,sep="")			
		} else{
			# NOTE: would probably benefit from using .Call and performing the whole process in C++ with a call to R for "yfun"...
			# y is padded with pre-returns
			x = c(prereturnsx, rep(0, n))
			y = c(prereturnsx, rep(0, n))
			# NOTE: sometimes the user may need to pass a vector which is longer since yfun may have
			# a lag structure (see base example with vDF dataset).
			yf = c(prereturnsy, rep(0, n))
			psim = matrix(0, ncol = modelinc[46], nrow = mar+n)
			xcondm = matrix(0, ncol = modelinc[46], nrow = mar+n)
			sigsim = simres = c(rep(0, mar+n))
			simres[1:mar] = preres
			for(j in 1:n){
				if(modelinc[48]==1){
					ytmp = c(yfun(as.numeric(yf[1:(yfn+j-1)])), NA)
					ytmp = tail(ytmp, mar+j)
				} else{
					ytmp = c(as.numeric(y[1:(mar+j-1)]),NA)	
				}
				XL = matrix(NA, ncol = length(ylags), nrow = mar+j)
				nXL = nrow(XL)
				for(k in 1:ncol(XL)){
					XL[(ylags[k]+1):nXL,k] = ytmp[1:(nXL-ylags[k])]
				}
				XL[is.na(XL)]=0
				arglist$XL = XL
				ptmp = switch(modelinc[46], dstar1path(arglist), dstar2path(arglist), dstar3path(arglist), dstar4path(arglist))
				psim[mar+j,] = tail(ptmp$probs, 1)
				msimx = switch(modelinc[46],
						psim[mar+j,1]*sig1+psim[mar+j,2]*sig2,
						psim[mar+j,1]*sig1+psim[mar+j,2]*sig2,
						psim[mar+j,1]*sig1+psim[mar+j,2]*sig2+psim[mar+j,3]*sig3,
						psim[mar+j,1]*sig1+psim[mar+j,2]*sig2+psim[mar+j,3]*sig3+psim[mar+j,4]*sig4)
				
				sigsim[mar+j] = msimx
				# zresidSim is size (burnin+n.sim) x m.sim (no mar)
				simres[mar+j] = zresidSim[j,i]*sigsim[mar+j]
				rsim = switch(modelinc[46],
						try(.C("starxsim1", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
										x = as.double(x), s = double(modelinc[46]*(mar+j)),res = as.double(simres[1:(mar+j)]),
										prob = as.double(psim[1:(mar+j),,drop=FALSE]), constm = as.double(constm[1:(mar+j),,drop=FALSE]),
										m = as.integer(mar+j-1), T = as.integer(mar+j), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim2", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
										x = as.double(x), s = double(modelinc[46]*(mar+j)),res = as.double(simres[1:(mar+j)]),
										prob = as.double(psim[1:(mar+j),,drop=FALSE]), constm = as.double(constm[1:(mar+j),,drop=FALSE]),
										m = as.integer(mar+j-1), T = as.integer(mar+j), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim3", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
										x = as.double(x), s = double(modelinc[46]*(mar+j)),res = as.double(simres[1:(mar+j)]),
										prob = as.double(psim[1:(mar+j),,drop=FALSE]), constm = as.double(constm[1:(mar+j),,drop=FALSE]),
										m = as.integer(mar+j-1), T = as.integer(mar+j), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim4", model = as.integer(modelinc), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
										x = as.double(x), s = double(modelinc[46]*(mar+j)),res = as.double(simres[1:(mar+j)]),
										prob = as.double(psim[1:(mar+j),,drop=FALSE]), constm = as.double(constm[1:(mar+j),,drop=FALSE]),
										m = as.integer(mar+j-1), T = as.integer(mar+j), PACKAGE="twinkle"), silent = TRUE)
				)
				x[mar+j] = rsim$x[mar+j]
				y[mar+j] = rsim$x[mar+j]
				yf[yfn+j] = y[mar+j]
				xcondm[mar+j,] = tail(matrix(rsim$s, ncol = modelinc[46]),1)
			}
			seriesSim[,i] = tail(x, n,sim)
			residSim[,i] = tail(simres, n.sim)
			sigmaSim[,i] = tail(sigsim, n.sim)
			probSim[[i]]  = tail(ptmp$probs, n.sim)
			condmSim[[i]] = tail(xcondm, n.sim)
			colnames(probSim[[i]]) = paste("state[",1:modelinc[46],"]",sep="")
			rownames(probSim[[i]]) = paste("T+",1:n.sim,sep="")
		}
	}
	rownames(residSim) = rownames(seriesSim) = paste("T+",1:n.sim,sep="")
	# check:
	# u=3
	# xtmp = (ipars["s1.phi0",1] + ipars["s1.phi1",1]*x[u-1]+ipars["s1.phi2",1]*x[u-2])*psim[u,1]+ (ipars["s2.phi0",1] + ipars["s2.phi1",1]*x[u-1]+ipars["s2.phi2",1]*x[u-2])*(psim[u,2]) + simres[u]
	sim = list(seriesSim = seriesSim, residSim = residSim, probSim = probSim, sigmaSim = sigmaSim)
	sol = new("STARpath",
			path = sim,
			model = model,
			seed = as.integer(sseed))
	return(sol)
}