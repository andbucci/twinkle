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


### ToDo: Need to rank the probabilities every time before inputing them into
# them into the matrix using the intercepts (if they exist as ranks..if they don't
# then assume they are zero).
.rollstar = function(spec, data, n.ahead = 1, forecast.length = 500, 
		n.start = NULL, refit.every = 25, refit.window = c("recursive", "moving"), 
		window.size = NULL, solver = "msoptim", fit.control = list(), solver.control = list(),
		calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05), cluster = NULL,
		keep.coef = TRUE, n=12, ...)
{
	tic = Sys.time()
	if(is.null(solver.control$trace)) trace = 0 else trace = solver.control$trace
	if(is.null(fit.control$stationarity)) fit.control$stationarity = TRUE
	if(is.null(fit.control$fixed.se)) fit.control$fixed.se = FALSE
	if(is.null(fit.control$rec.init)) fit.control$rec.init = 'all'
	mm = match(names(fit.control), c("stationarity", "fixed.se", "rec.init"))
	if(any(is.na(mm))){
		idx = which(is.na(mm))
		enx = NULL
		for(i in 1:length(idx)) enx = c(enx, names(fit.control)[idx[i]])
		warning(paste(c("unidentified option(s) in fit.control:\n", enx), sep="", collapse=" "), call. = FALSE, domain = NULL)
	}
	refit.window = refit.window[1]
	xdata = .extractdata(data)
	data = xdata$data
	dindex = xdata$index
	period = xdata$period
	T = NROW(data)
	modelinc = spec@model$modelinc
	model = spec@model
	if(modelinc[47]==1){
		chk = all.equal(dindex, index(model$fixed.prob))
		if(!is.logical(chk) | chk == FALSE){
			print(paste("\n",chk,sep=""))
			stop("\nrollstar-->error: data and fixed.probs indices do not match\n")
		}
		fprobs = model$fixed.prob
		fex = TRUE
	} else{
		fprobs = NULL
		fex = FALSE
	}
	if(modelinc[3] > 0){
		chk = all.equal(dindex, index(model$modeldata$mexdata))
		if(!is.logical(chk) | chk == FALSE){
			print(paste("\n",chk,sep=""))
			stop("\nrollstar-->error: data and external.regressor (mean) indices do not match\n")
		}
		mexdata = model$modeldata$mexdata
		mex = TRUE
	} else{
		mexdata = NULL
		mex = FALSE
	}
	if(modelinc[49]==2){
		chk = all.equal(dindex, index(model$modeldata$s))
		if(!is.logical(chk) | chk == FALSE){
			print(paste("\n",chk,sep=""))
			stop("\nrollstar-->error: data and 's' probability dynamics regressor indices do not match\n")
		}
		sxdata = model$modeldata$s
		sxex = TRUE
	} else{
		sxdata = NULL
		sxex = FALSE
	}
	if(modelinc[39]>0){
		chk = all.equal(dindex, index(model$modeldata$vexdata))
		if(!is.logical(chk) | chk == FALSE){
			print(paste("\n",chk,sep=""))
			stop("\nrollstar-->error: data and external.regressor (variance) indices do not match\n")
		}
		vexdata = model$modeldata$vexdata
		vex = TRUE
	} else{
		vexdata = NULL
		vex = FALSE
	}
	if(n.ahead>1) stop("\nrollstar:--> n.ahead>1 not supported...try again.")
	if(is.null(n.start)){
		if(is.null(forecast.length)) stop("\nrollstar:--> forecast.length amd n.start are both NULL....try again.")
		n.start = T - forecast.length
	} else{
		forecast.length = T - n.start
	}
	if(T<=n.start) stop("\nrollstar:--> start cannot be greater than length of data")
	# the ending points of the estimation window
	s = seq(n.start+refit.every, T, by = refit.every)
	m = length(s)
	# the rolling forecast length
	out.sample = rep(refit.every, m)
	# adjustment to include all the datapoints from the end
	if(s[m]<T){
		s = c(s,T)
		m = length(s)
		out.sample = c(out.sample, s[m]-s[m-1])
	}
	if(refit.window == "recursive"){
		rollind = lapply(1:m, FUN = function(i) 1:s[i])
	} else{
		if(!is.null(window.size)){
			if(window.size<100) stop("\nrollstar:--> window size must be greater than 100.")
			rollind = lapply(1:m, FUN = function(i) max(1, (s[i]-window.size-out.sample[i])):s[i])
		} else{
			rollind = lapply(1:m, FUN = function(i) (1+(i-1)*refit.every):s[i])
		}
	}
	# distribution
	distribution = spec@model$modeldesc$distribution
	if(any(distribution==c("snorm","sstd","sged","jsu","nig","ghyp","ghst"))) skewed = TRUE else skewed = FALSE
	if(any(distribution==c("std","sstd","ged","sged","jsu","nig","ghyp","ghst"))) shaped = TRUE else shaped = FALSE
	if(any(distribution==c("ghyp"))) ghyp = TRUE else ghyp = FALSE
	
	if( !is.null(cluster) ){
		clusterEvalQ(cl = cluster, library(twinkle))
		clusterExport(cluster, c("data", "dindex", "s","refit.every", 
						"keep.coef", "shaped", "skewed", "ghyp", 
						"rollind", "spec", "out.sample", "mex", "vex", "sxex", "fex",
						"solver", "solver.control", "fit.control", "n"), envir = environment())
		if(mex) clusterExport(cluster, c("mexdata"), envir = environment())
		if(vex) clusterExport(cluster, c("vexdata"), envir = environment())
		if(fex) clusterExport(cluster, c("fprobs"), envir = environment())
		if(sxex) clusterExport(cluster, c("sxdata"), envir = environment())
		
		tmp = parLapply(cl = cluster, 1:m, fun = function(i){
					if(mex) spec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
					if(vex) spec@model$modeldata$vexdata = vexdata[rollind[[i]],,drop=FALSE]
					if(sxex) spec@model$modeldata$s = sxdata[rollind[[i]],,drop=FALSE]
					if(fex) spec@model$fixed.prob = fprobs[rollind[[i]],,drop=FALSE]
					fit = try(starfit(spec, xts::xts(data[rollind[[i]]], dindex[rollind[[i]]]), out.sample = out.sample[i], 
									solver = solver, solver.control = solver.control, 
									fit.control = fit.control, n = n), silent=TRUE)
					# 3 cases: General Error, Failure to Converge, Failure to invert Hessian (bad solution)
					if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
						ans = list(y = NA, cf = NA, converge = FALSE, loglik = NA)
					} else{
						if(mex) fmex = tail(mexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fmex = NULL
						if(vex) fvex = tail(vexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fvex = NULL
						if(sxex) fsxex = tail(sxdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fsxex = NULL
						if(fex) ffex = tail(fprobs[rollind[[i]],,drop=FALSE], out.sample[i]) else ffex = NULL
						f = starforecast(fit, n.ahead = 1, n.roll = out.sample[i]-1, 
								external.forecasts = list(xregfor = fmex, vregfor = fvex, sfor = fsxex,
										probfor = ffex))
						ret = as.numeric( fitted(f) )
						if(spec@model$modelinc[50]>0) sig = as.numeric( sigma(f) ) else sig = rep(coef(fit)["sigma"], length(ret))
						if(shaped) shp = rep(coef(fit)["shape"], out.sample[i]) else shp = rep(0, out.sample[i])
						if(skewed) skw = rep(coef(fit)["skew"], out.sample[i]) else skw = rep(0, out.sample[i])
						if(ghyp) shpgig = rep(coef(fit)["ghlambda"], out.sample[i]) else shpgig = rep(0, out.sample[i])
						if(!fex) fp = coredata(states(f)) else fp = coredata(ffex)
						rlz = tail(data[rollind[[i]]], out.sample[i])
						# use xts for indexing the forecasts
						y = as.data.frame(cbind(ret, sig, skw, shp, shpgig, fp, rlz))
						stx = paste("Prob[State=",1:ncol(fp),"]",sep="")
						rownames(y) = tail(as.character(dindex[rollind[[i]]]), out.sample[i])						
						colnames(y) = c("Mu", "Sigma", "Skew", "Shape", "Shape(GIG)", stx, "Realized")
						if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
						ans = list(y = y, cf = cf, converge = TRUE, loglik = likelihood(fit))
					}
					return(ans)})
	} else{
		tmp = lapply(as.list(1:m), FUN = function(i){
					if(mex) spec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
					if(vex) spec@model$modeldata$vexdata = vexdata[rollind[[i]],,drop=FALSE]
					if(sxex) spec@model$modeldata$s = sxdata[rollind[[i]],,drop=FALSE]
					if(fex) spec@model$fixed.prob = fprobs[rollind[[i]],,drop=FALSE]
					fit = try(starfit(spec, xts(data[rollind[[i]]], dindex[rollind[[i]]]), out.sample = out.sample[i], 
									solver = solver, solver.control = solver.control, 
									fit.control = fit.control, n = n), silent=TRUE)
					if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
						ans = list(y = NA, cf = NA, converge = FALSE, loglik = NA)
					} else{
						if(mex) fmex = tail(mexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fmex = NULL
						if(vex) fvex = tail(vexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fvex = NULL
						if(sxex) fsxex = tail(sxdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fsxex = NULL
						if(fex) ffex = tail(fprobs[rollind[[i]],,drop=FALSE], out.sample[i]) else ffex = NULL
						f = starforecast(fit, n.ahead = 1, n.roll = out.sample[i]-1, 
								external.forecasts = list(xregfor = fmex, vregfor = fvex, sfor = fsxex,
										probfor = ffex))
						ret = as.numeric( fitted(f) )
						if(spec@model$modelinc[50]>0) sig = as.numeric( sigma(f) ) else sig = rep(coef(fit)["sigma"], length(ret))
						if(shaped) shp = rep(coef(fit)["shape"], out.sample[i]) else shp = rep(0, out.sample[i])
						if(skewed) skw = rep(coef(fit)["skew"], out.sample[i]) else skw = rep(0, out.sample[i])
						if(ghyp) shpgig = rep(coef(fit)["ghlambda"], out.sample[i]) else shpgig = rep(0, out.sample[i])
						if(!fex) fp = coredata(states(f)) else fp = coredata(ffex)
						rlz = tail(data[rollind[[i]]], out.sample[i])
						# use xts for indexing the forecasts
						y = as.data.frame(cbind(ret, sig, skw, shp, shpgig, fp, rlz))
						stx = paste("Prob[State=",1:ncol(fp),"]",sep="")
						rownames(y) = tail(as.character(dindex[rollind[[i]]]), out.sample[i])						
						colnames(y) = c("Mu", "Sigma", "Skew", "Shape", "Shape(GIG)", stx, "Realized")
						if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
						ans = list(y = y, cf = cf, converge = TRUE, loglik = likelihood(fit))
					}
					return(ans)})
	}
	conv = sapply(tmp, FUN = function(x) x$converge)
	if(any(!conv)){
		warning("\nnon-converged estimation windows present...resubsmit object with different solver parameters...")
		noncidx = which(!conv)
		model = list()
		model$spec = spec
		model$data = data
		model$index = dindex
		model$period = period
		model$n.ahead = n.ahead
		model$forecast.length = forecast.length 
		model$n.start = n.start
		model$n.refits = m
		model$refit.every = refit.every
		model$refit.window = refit.window
		model$window.size = window.size
		model$calculate.VaR = calculate.VaR
		model$VaR.alpha = VaR.alpha
		model$keep.coef = keep.coef
		model$noncidx = noncidx
		model$rollind = rollind
		model$out.sample = out.sample
		forecast = tmp
		toc = Sys.time()-tic
		model$elapsed = toc
		ans = new("STARroll",
				model = model,
				forecast = forecast)
		return(ans)
	} else{
		noncidx = NULL
		forc = tmp[[1]]$y
		for(i in 2:m){
			forc = rbind(forc, tmp[[i]]$y)
		}
		cf = vector(mode = "list", length = m)
		for(i in 1:m){
			cf[[i]]$index = dindex[tail(rollind[[i]],1) - out.sample[i]]
			cf[[i]]$coef = tmp[[i]]$cf
		}
		if(calculate.VaR){
			if(is.null(VaR.alpha)) VaR.alpha = c(0.01, 0.05)
			VaR.matrix = matrix(NA, ncol = length(VaR.alpha)+1, nrow = NROW(forc))
			for(i in 1:length(VaR.alpha)){
				VaR.matrix[,i] = qdist(VaR.alpha[i], mu = forc[,1], sigma = forc[,2], 
						skew = forc[,3], shape = forc[,4], lambda = forc[,5], 
						distribution = spec@model$modeldesc$distribution)
			}
			VaR.matrix[,length(VaR.alpha)+1] = forc[,ncol(forc)]
			colnames(VaR.matrix) = c(paste("alpha(", round(VaR.alpha,2)*100, "%)",sep=""), "realized")
			VaR.matrix = as.data.frame(VaR.matrix)
			rownames(VaR.matrix) = rownames(forc)
		} else{
			VaR.matrix = NULL
		}
		model = list()
		model$spec = spec
		model$data = data
		model$index = dindex
		model$period = period
		model$n.ahead = n.ahead
		model$forecast.length = forecast.length 
		model$n.start = n.start
		model$refit.every = refit.every
		model$n.refits = m
		model$refit.window = refit.window
		model$window.size = window.size
		model$calculate.VaR = calculate.VaR
		model$VaR.alpha = VaR.alpha
		model$keep.coef = keep.coef
		model$noncidx = noncidx
		model$coef = cf
		model$rollind = rollind
		model$out.sample = out.sample
		model$loglik = sapply(tmp, function(x) x$loglik)
		forecast = list(VaR = VaR.matrix, density = forc)
	}
	toc = Sys.time()-tic
	model$elapsed = toc
	ans = new("STARroll",
			model = model,
			forecast = forecast)
	return( ans )
}


.resumerollstar = function(object, solver = "msoptim", fit.control = list(), solver.control = list(), cluster = NULL, n = 12)
{
	if(!is.null(object@model$noncidx)){
		noncidx = object@model$noncidx
		tic = Sys.time()
		if(is.null(solver.control$trace)) trace = 0 else trace = solver.control$trace
		if(is.null(fit.control$stationarity)) fit.control$stationarity = TRUE
		if(is.null(fit.control$fixed.se)) fit.control$fixed.se = FALSE
		if(is.null(fit.control$rec.init)) fit.control$rec.init = 'all'
		mm = match(names(fit.control), c("stationarity", "fixed.se", "rec.init"))
		if(any(is.na(mm))){
			idx = which(is.na(mm))
			enx = NULL
			for(i in 1:length(idx)) enx = c(enx, names(fit.control)[idx[i]])
			warning(paste(c("unidentified option(s) in fit.control:\n", enx), sep="", collapse=" "), call. = FALSE, domain = NULL)
		}
		model = object@model
		keep.coef = model$keep.coef
		spec = model$spec
		datanames = model$datanames
		data = model$data
		dindex = model$index
		period= model$period
		T = NROW(data)
		modelinc = spec@model$modelinc
		calculate.VaR = model$calculate.VaR
		VaR.alpha = model$VaR.alpha
		if(modelinc[47]==1){
			chk = all.equal(dindex, index(spec@model$fixed.prob))
			if(!is.logical(chk) | chk == FALSE){
				print(paste("\n",chk,sep=""))
				stop("\nrollstar-->error: data and fixed.probs indices do not match\n")
			}
			fprobs = spec@model$fixed.prob
			fex = TRUE
		} else{
			fprobs = NULL
			fex = FALSE
		}
		if(modelinc[3] > 0){
			chk = all.equal(dindex, index(spec@model$modeldata$mexdata))
			if(!is.logical(chk) | chk == FALSE){
				print(paste("\n",chk,sep=""))
				stop("\nrollstar-->error: data and external.regressor (mean) indices do not match\n")
			}
			mexdata = spec@model$modeldata$mexdata
			mex = TRUE
		} else{
			mexdata = NULL
			mex = FALSE
		}
		if(modelinc[49]==2){
			chk = all.equal(dindex, index(spec@model$modeldata$s))
			if(!is.logical(chk) | chk == FALSE){
				print(paste("\n",chk,sep=""))
				stop("\nrollstar-->error: data and 's' probability dynamics regressor indices do not match\n")
			}
			sxdata = spec@model$modeldata$s
			sxex = TRUE
		} else{
			sxdata = NULL
			sxex = FALSE
		}
		if(modelinc[39]>0){
			chk = all.equal(dindex, index(spec@model$modeldata$vexdata))
			if(!is.logical(chk) | chk == FALSE){
				print(paste("\n",chk,sep=""))
				stop("\nrollstar-->error: data and external.regressor (variance) indices do not match\n")
			}
			vexdata = spec@model$modeldata$vexdata
			vex = TRUE
		} else{
			vexdata = NULL
			vex = FALSE
		}
		n.ahead = model$n.ahead
		n.start = model$n.start
		forecast.length = model$forecast.length
		refit.every = model$refit.every
		refit.window = model$refit.window
		window.size = model$window.size
		if(n.ahead>1) stop("\nrollstar:--> n.ahead>1 not supported...try again.")
		if(is.null(n.start)){
			if(is.null(forecast.length)) stop("\nrollstar:--> forecast.length amd n.start are both NULL....try again.")
			n.start = T - forecast.length
		}
		if(T<=n.start) stop("\nrollstar:--> start cannot be greater than length of data")
		# the ending points of the estimation window
		s = seq(n.start+refit.every, T, by = refit.every)
		m = length(s)
		# the rolling forecast length
		out.sample = rep(refit.every, m)
		# adjustment to include all the datapoints from the end
		if(s[m]<T){
			s = c(s,T)
			m = length(s)
			out.sample = c(out.sample, s[m]-s[m-1])
		}
		if(refit.window == "recursive"){
			rollind = lapply(1:m, FUN = function(i) 1:s[i])
		} else{
			if(!is.null(window.size)){
				if(window.size<100) stop("\nrollstar:--> window size must be greater than 100.")
				rollind = lapply(1:m, FUN = function(i) max(1, (s[i]-window.size-out.sample[i])):s[i])
			} else{
				rollind = lapply(1:m, FUN = function(i) (1+(i-1)*refit.every):s[i])
			}
		}
		# distribution
		distribution = spec@model$modeldesc$distribution
		if(any(distribution==c("snorm","sstd","sged","jsu","nig","ghyp","ghst"))) skewed = TRUE else skewed = FALSE
		if(any(distribution==c("std","sstd","ged","sged","jsu","nig","ghyp","ghst"))) shaped = TRUE else shaped = FALSE
		if(any(distribution==c("ghyp"))) ghyp = TRUE else ghyp = FALSE
		if( !is.null(cluster) ){
			clusterEvalQ(cl = cluster, library(twinkle))
			clusterExport(cluster, c("data", "dindex","s","refit.every",
							"keep.coef", "shaped", "skewed", "ghyp", 
							"rollind", "spec", "out.sample", "mex", "vex", "sxex", "fex",
							"solver", "solver.control", "fit.control", "n"), envir = environment())
			if(mex) clusterExport(cluster, c("mexdata"), envir = environment())
			if(vex) clusterExport(cluster, c("vexdata"), envir = environment())
			if(fex) clusterExport(cluster, c("fprobs"), envir = environment())
			if(sxex) clusterExport(cluster, c("sxdata"), envir = environment())
			tmp = parLapply(cl = cluster, as.list(noncidx), fun = function(i){
						if(mex) spec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
						if(vex) spec@model$modeldata$vexdata = vexdata[rollind[[i]],,drop=FALSE]
						if(sxex) spec@model$modeldata$s = sxdata[rollind[[i]],,drop=FALSE]
						if(fex) spec@model$fixed.prob = fprobs[rollind[[i]],,drop=FALSE]
						fit = try(starfit(spec, xts::xts(data[rollind[[i]]], dindex[rollind[[i]]]), out.sample = out.sample[i], 
										solver = solver, solver.control = solver.control, 
										fit.control = fit.control, n = n), silent=TRUE)
						# 3 cases: General Error, Failure to Converge, Failure to invert Hessian (bad solution)
						if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
							ans = list(y = NA, cf = NA, converge = FALSE, loglik = NA)
						} else{
							if(mex) fmex = tail(mexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fmex = NULL
							if(vex) fvex = tail(vexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fvex = NULL
							if(sxex) fsxex = tail(sxdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fsxex = NULL
							if(fex) ffex = tail(fprobs[rollind[[i]],,drop=FALSE], out.sample[i]) else ffex = NULL
							f = starforecast(fit, n.ahead = 1, n.roll = out.sample[i]-1, 
									external.forecasts = list(xregfor = fmex, vregfor = fvex, sfor = fsxex,
											probfor = ffex))
							ret = as.numeric( fitted(f) )
							if(spec@model$modelinc[50]>0) sig = as.numeric( sigma(f) ) else sig = rep(coef(fit)["sigma"], length(ret))
							if(shaped) shp = rep(coef(fit)["shape"], out.sample[i]) else shp = rep(0, out.sample[i])
							if(skewed) skw = rep(coef(fit)["skew"], out.sample[i]) else skw = rep(0, out.sample[i])
							if(ghyp) shpgig = rep(coef(fit)["ghlambda"], out.sample[i]) else shpgig = rep(0, out.sample[i])
							if(!fex) fp = coredata(states(f)) else fp = coredata(ffex)
							rlz = tail(data[rollind[[i]]], out.sample[i])
							# use xts for indexing the forecasts
							y = as.data.frame(cbind(ret, sig, skw, shp, shpgig, fp, rlz))
							stx = paste("Prob[State=",1:ncol(fp),"]",sep="")
							rownames(y) = tail(as.character(dindex[rollind[[i]]]), out.sample[i])						
							colnames(y) = c("Mu", "Sigma", "Skew", "Shape", "Shape(GIG)", stx, "Realized")
							if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
							ans = list(y = y, cf = cf, converge = TRUE, loglik = likelihood(fit))
						}
						return(ans)})
		} else{
			tmp = lapply(as.list(noncidx), FUN = function(i){
						if(mex) spec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
						if(vex) spec@model$modeldata$vexdata = vexdata[rollind[[i]],,drop=FALSE]
						if(sxex) spec@model$modeldata$s = sxdata[rollind[[i]],,drop=FALSE]
						if(fex) spec@model$fixed.prob = fprobs[rollind[[i]],,drop=FALSE]
						fit = try(starfit(spec, xts(data[rollind[[i]]], dindex[rollind[[i]]]), out.sample = out.sample[i], 
										solver = solver, solver.control = solver.control, 
										fit.control = fit.control, n = n), silent=TRUE)
						if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
							ans = list(y = NA, cf = NA, converge = FALSE, loglik = NA)
						} else{
							if(mex) fmex = tail(mexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fmex = NULL
							if(vex) fvex = tail(vexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fvex = NULL
							if(sxex) fsxex = tail(sxdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fsxex = NULL
							if(fex) ffex = tail(fprobs[rollind[[i]],,drop=FALSE], out.sample[i]) else ffex = NULL
							f = starforecast(fit, n.ahead = 1, n.roll = out.sample[i]-1, 
									external.forecasts = list(xregfor = fmex, vregfor = fvex, sfor = fsxex,
											probfor = ffex))
							ret = as.numeric( fitted(f) )
							if(spec@model$modelinc[50]>0) sig = as.numeric( sigma(f) ) else sig = rep(coef(fit)["sigma"], length(ret))
							if(shaped) shp = rep(coef(fit)["shape"], out.sample[i]) else shp = rep(0, out.sample[i])
							if(skewed) skw = rep(coef(fit)["skew"], out.sample[i]) else skw = rep(0, out.sample[i])
							if(ghyp) shpgig = rep(coef(fit)["ghlambda"], out.sample[i]) else shpgig = rep(0, out.sample[i])
							if(!fex) fp = coredata(states(f)) else fp = coredata(ffex)
							rlz = tail(data[rollind[[i]]], out.sample[i])
							# use xts for indexing the forecasts
							y = as.data.frame(cbind(ret, sig, skw, shp, shpgig, fp, rlz))
							stx = paste("Prob[State=",1:ncol(fp),"]",sep="")
							rownames(y) = tail(as.character(dindex[rollind[[i]]]), out.sample[i])						
							colnames(y) = c("Mu", "Sigma", "Skew", "Shape", "Shape(GIG)", stx, "Realized")
							if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
							ans = list(y = y, cf = cf, converge = TRUE, loglik = likelihood(fit))
						}
						return(ans)})
		}		
		forecast = object@forecast
		conv = sapply(tmp, FUN = function(x) x$converge)
		for(i in 1:length(noncidx)){
			if(conv[i]) forecast[[noncidx[i]]] = tmp[[i]]
		}
		if(any(!conv)){
			warning("\nnon-converged estimation windows present...resubsmit object with different solver parameters...")
			noncidx = noncidx[which(!conv)]
			model = list()
			model$spec = spec
			model$data = data
			model$index = dindex
			model$period = period
			model$n.ahead = n.ahead
			model$forecast.length = forecast.length 
			model$n.start = n.start
			model$refit.every = refit.every
			model$n.refits = m
			model$refit.window = refit.window
			model$window.size = window.size
			model$calculate.VaR = calculate.VaR
			model$VaR.alpha = VaR.alpha
			model$keep.coef = keep.coef
			model$noncidx = noncidx
			model$rollind = rollind
			model$out.sample = out.sample
			forecast = forecast
			toc = Sys.time()-tic
			model$elapsed = toc
			ans = new("STARroll",
					model = model,
					forecast = forecast)
			return( ans )			
		} else{
			noncidx = NULL
			forc = forecast[[1]]$y
			for(i in 2:m){
				forc = rbind(forc, forecast[[i]]$y)
			}
			cf = vector(mode = "list", length = m)
			for(i in 1:m){
				cf[[i]]$index = dindex[tail(rollind[[i]],1) - out.sample[i]]
				cf[[i]]$coef = forecast[[i]]$cf
			}
			if(calculate.VaR){
				if(is.null(VaR.alpha)) VaR.alpha = c(0.01, 0.05)
				VaR.matrix = matrix(NA, ncol = length(VaR.alpha)+1, nrow = NROW(forc))
				for(i in 1:length(VaR.alpha)){
					VaR.matrix[,i] = qdist(VaR.alpha[i], mu = forc[,1], sigma = forc[,2], 
							skew = forc[,3], shape = forc[,4], lambda = forc[,5], 
							distribution = spec@model$modeldesc$distribution)
				}
				VaR.matrix[,length(VaR.alpha)+1] = forc[,ncol(forc)]
				colnames(VaR.matrix) = c(paste("alpha(", round(VaR.alpha,2)*100, "%)",sep=""), "realized")
				VaR.matrix = as.data.frame(VaR.matrix)
				rownames(VaR.matrix) = rownames(forc)
			} else{
				VaR.matrix = NULL
			}
			model = list()
			model$spec = spec
			model$data = data
			model$index = dindex
			model$period = period
			model$n.ahead = n.ahead
			model$forecast.length = forecast.length 
			model$n.start = n.start
			model$refit.every = refit.every
			model$n.refits = m
			model$refit.window = refit.window
			model$window.size = window.size
			model$calculate.VaR = calculate.VaR
			model$VaR.alpha = VaR.alpha
			model$keep.coef = keep.coef
			model$noncidx = noncidx
			model$rollind = rollind
			model$out.sample = out.sample
			model$coef = cf
			model$loglik = sapply(tmp, function(x) x$loglik)
			forecast = list(VaR = VaR.matrix, density = forc)
		}
		toc = Sys.time()-tic
		model$elapsed = toc
		ans = new("STARroll",
				model = model,
				forecast = forecast)
	} else{
		# do nothing...all converged
		ans = object
	}
	return( ans )
}
