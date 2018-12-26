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
.stars1sLLH = function(pars, arglist)
{
	modelinc = arglist$model$modelinc
	if(modelinc[47]==1){
		probs = arglist$probs
	} else{
		probs = rep(1, length(arglist$data))
	}
	data = arglist$data
	returnType = arglist$returnType
	starenv = arglist$starenv
	# rejoin fixed and pars
	model = arglist$model
	estidx = arglist$estidx
	idx = model$pidx
	ipars = arglist$ipars	
	ipars[estidx, 1] = pars
	if(arglist$transform==1){
		if(modelinc[1]>0 && ipars["s1.phi0",4]==1) ipars["s1.phi0",1] = logtransform(ipars["s1.phi0",1], ipars["s1.phi0","LB"], ipars["s1.phi0","UB"])
		if(modelinc[2]>0 && ipars["s1.phi1",4]==1) ipars["s1.phi1",1] = logtransform(ipars["s1.phi1",1], ipars["s1.phi1","LB"], ipars["s1.phi1","UB"])
		if(modelinc[4]>0 && ipars["s1.psi1",4]==1) ipars["s1.psi1",1] = logtransform(ipars["s1.psi1",1], ipars["s1.psi1","LB"], ipars["s1.psi1","UB"])		
		if(modelinc[30]>0 && ipars["sigma",4]==1) ipars["sigma",1] = logtransform(ipars["sigma",1], ipars["sigma","LB"], ipars["sigma","UB"])
		if(modelinc[41]>0 && ipars["skew",4]==1) ipars["skew",1] = logtransform(ipars["skew",1], ipars["skew","LB"], ipars["skew","UB"])
		if(modelinc[42]>0 && ipars["shape",4]==1) ipars["shape",1] = logtransform(ipars["shape",1], ipars["shape","LB"], ipars["shape","UB"])
		if(modelinc[43]>0 && ipars["ghlambda",4]==1) ipars["ghlambda",1] = logtransform(ipars["ghlambda",1], ipars["ghlambda","LB"], ipars["ghlambda","UB"])
	}
	trace = arglist$trace
	T = length(data)
	mexdata = arglist$mexdata
	fit.control = arglist$fit.control	
	distribution = model$modeldesc$distribution
	# distribution number
	# 1 = norm, 2=snorm, 3=std, 4=sstd, 5=ged, 6=sged, 7=nig
	dist = model$modeldesc$distno
	if(modelinc[3]>0) mexdata = as.double(as.vector(mexdata)) else mexdata = double(1)
	
	ans = try( .C("starxfilters1s", model = as.integer(modelinc[1:60]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					x = as.double(data), res = double(T), z = double(T), 
					mexdata = mexdata, prob = as.double(probs), 
					constm = double(T), condm = double(T), marx = double(T), T = as.integer(T), 
					LHT = double(T), llh = double(1),  PACKAGE = "twinkle"), silent = TRUE )
	if(inherits(ans, "try-error")){
		if(arglist$pmode!=1){
			assign("twinkle_csol", 1, envir = starenv)
			assign("twinkle_filtermessage", ans, envir = starenv)
			if( trace > 0 ) cat(paste("\nstarfit-->warning: ", get("twinkle_filtermessage", starenv),"\n", sep=""))
			return(llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))))
		} else{
			return(llh = 1e10)
		}
	} else{
		if(arglist$pmode!=1){
			assign("twinkle_csol", 0, envir = starenv)
		}
	}
	z = ans$z
	epsx = ans$res
	llh = ans$llh
	
	if(is.finite(llh) && !is.na(llh) && !is.nan(llh)){
		if(arglist$pmode!=1) assign("star_llh", llh, envir = starenv)
	} else {
		if(arglist$pmode!=1) llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))) else llh = 1e10
	}
	
	# LHT = raw scores
	# LHT = -ans$LHT[(m):T]
	LHT = -ans$LHT
	ans = switch(returnType,
			llh = llh,
			LHT = LHT,
			all = list(llh = llh, res = epsx, z = z, LHT = LHT, probs = probs, condm = matrix(ans$condm, ncol=1),
					constm = matrix(ans$constm, ncol=1), ipars = ipars, pmu = NULL))
	return(ans)
}

.stars2sLLH = function(pars, arglist)
{
	modelinc = arglist$model$modelinc
	data = arglist$data
	returnType = arglist$returnType
	starenv = arglist$starenv
	# rejoin fixed and pars
	model = arglist$model
	estidx = arglist$estidx
	idx = model$pidx
	ipars = arglist$ipars	
	ipars[estidx, 1] = pars
	if(arglist$transform==1){
		iidx = which(ipars[,4]==1 & ipars[,7]==1)
		if(length(iidx)>0){
			ipars[iidx,1] = logtransform(ipars[iidx,1], ipars[iidx,5], ipars[iidx,6])
		}
	}
	# Transform to positive
	ipars["s1.gamma",1] = abs(ipars["s1.gamma",1])
	if(modelinc[47]==1){
		probs = arglist$probs
		pmu = NULL
	} else{
		xprobs = dstar2(ipars, arglist)
		probs = xprobs$probs
		pmu = xprobs$pmu
	}
	trace = arglist$trace
	T = length(data)
	mexdata = arglist$mexdata
	fit.control = arglist$fit.control	
	distribution = model$modeldesc$distribution
	# distribution number
	# 1 = norm, 2=snorm, 3=std, 4=sstd, 5=ged, 6=sged, 7=nig
	dist = model$modeldesc$distno
	if(modelinc[3]>0) mexdata = as.double(as.vector(mexdata)) else mexdata = double(1)
	
	ans = try( .C("starxfilters2s", model = as.integer(modelinc[1:60]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					x = as.double(data), res = double(T), z = double(T), 
					mexdata = mexdata, prob = as.double(probs), 
					constm = double(2*T), condm = double(2*T), marx = double(T),
					T = as.integer(T), LHT = double(T), llh = double(1),  PACKAGE = "twinkle"), silent = TRUE )
	if(inherits(ans, "try-error")){
		if(arglist$pmode!=1){
			assign("twinkle_csol", 1, envir = starenv)
			assign("twinkle_filtermessage", ans, envir = starenv)
			if( trace > 0 ) cat(paste("\nstarfit-->warning: ", get("twinkle_filtermessage", starenv),"\n", sep=""))
			return(llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))))
		} else{
			return(llh = 1e10)
		}
	} else{
		if(arglist$pmode!=1){
			assign("twinkle_csol", 0, envir = starenv)
		}
	}
	z = ans$z
	epsx = ans$res
	llh = ans$llh
	
	if(is.finite(llh) && !is.na(llh) && !is.nan(llh)){
		if(arglist$pmode!=1) assign("star_llh", llh, envir = starenv)
	} else {
		if(arglist$pmode!=1) llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))) else llh = 1e10
	}
	
	# LHT = raw scores
	# LHT = -ans$LHT[(m):T]
	LHT = -ans$LHT
	ans = switch(returnType,
			llh = llh,
			LHT = LHT,
			all = list(llh = llh, res = epsx, z = z, LHT = LHT, probs = probs, 
					condm = matrix(ans$condm, ncol=2), constm = matrix(ans$constm, ncol=2), 
					ipars = ipars, pmu = pmu))
	return(ans)
}

.stars3sLLH = function(pars, arglist)
{
	modelinc = arglist$model$modelinc
	data = arglist$data
	returnType = arglist$returnType
	starenv = arglist$starenv
	# rejoin fixed and pars
	model = arglist$model
	estidx = arglist$estidx
	idx = model$pidx
	ipars = arglist$ipars	
	ipars[estidx, 1] = pars
	if(arglist$transform==1){
		iidx = which(ipars[,4]==1 & ipars[,7]==1)
		if(length(iidx)>0){
			ipars[iidx,1] = logtransform(ipars[iidx,1], ipars[iidx,5], ipars[iidx,6])
		}
	}
	# Transform to positive
	ipars["s1.gamma",1] = abs(ipars["s1.gamma",1])
	ipars["s2.gamma",1] = abs(ipars["s2.gamma",1])
	if(modelinc[47]==1){
		probs = arglist$probs
		pmu = NULL
	} else{
		xprobs = dstar3(ipars, arglist)
		probs = xprobs$probs
		pmu = xprobs$pmu
	}
	trace = arglist$trace
	T = length(data)
	mexdata = arglist$mexdata
	fit.control = arglist$fit.control
	distribution = model$modeldesc$distribution
	# distribution number
	# 1 = norm, 2=snorm, 3=std, 4=sstd, 5=ged, 6=sged, 7=nig
	dist = model$modeldesc$distno
	if(modelinc[3]>0) mexdata = as.double(as.vector(mexdata)) else mexdata = double(1)
	
	ans = try( .C("starxfilters3s", model = as.integer(modelinc[1:60]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					x = as.double(data), res = double(T), z = double(T), 
					mexdata = mexdata, prob = as.double(probs), 
					constm = double(3*T), condm = double(3*T),  marx = double(T),
					T = as.integer(T), LHT = double(T), llh = double(1),  PACKAGE = "twinkle"), silent = TRUE )
	if(inherits(ans, "try-error")){
		if(arglist$pmode!=1){
			assign("twinkle_csol", 1, envir = starenv)
			assign("twinkle_filtermessage", ans, envir = starenv)
			if( trace > 0 ) cat(paste("\nstarfit-->warning: ", get("twinkle_filtermessage", starenv),"\n", sep=""))
			return(llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))))
		} else{
			return(llh = 1e10)
		}
	} else{
		if(arglist$pmode!=1){
			assign("twinkle_csol", 0, envir = starenv)
		}
	}
	z = ans$z
	epsx = ans$res
	llh = ans$llh
	
	if(is.finite(llh) && !is.na(llh) && !is.nan(llh)){
		if(arglist$pmode!=1) assign("star_llh", llh, envir = starenv)
	} else {
		if(arglist$pmode!=1) llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))) else llh = 1e10
	}
	
	# LHT = raw scores
	# LHT = -ans$LHT[(m):T]
	LHT = -ans$LHT
	ans = switch(returnType,
			llh = llh,
			LHT = LHT,
			all = list(llh = llh, res = epsx, z = z, LHT = LHT, probs = probs, condm = matrix(ans$condm, ncol = 3), 
					constm = matrix(ans$constm, ncol=3), ipars = ipars, pmu = pmu))
	return(ans)
}

.stars4sLLH = function(pars, arglist)
{
	modelinc = arglist$model$modelinc
	data = arglist$data
	returnType = arglist$returnType
	starenv = arglist$starenv
	# rejoin fixed and pars
	model = arglist$model
	estidx = arglist$estidx
	idx = model$pidx
	ipars = arglist$ipars	
	ipars[estidx, 1] = pars
	if(arglist$transform==1){
		iidx = which(ipars[,4]==1 & ipars[,7]==1)
		if(length(iidx)>0){
			ipars[iidx,1] = logtransform(ipars[iidx,1], ipars[iidx,5], ipars[iidx,6])
		}
	}
	# Transform to positive
	ipars["s1.gamma",1] = abs(ipars["s1.gamma",1])
	ipars["s2.gamma",1] = abs(ipars["s2.gamma",1])
	ipars["s3.gamma",1] = abs(ipars["s3.gamma",1])
	
	if(modelinc[47]==1){
		probs = arglist$probs
		pmu = NULL
	} else{
		xprobs = dstar4(ipars, arglist)
		probs = xprobs$probs
		pmu = xprobs$pmu
	}
	trace = arglist$trace
	T = length(data)
	mexdata = arglist$mexdata
	fit.control = arglist$fit.control
	distribution = model$modeldesc$distribution
	# distribution number
	# 1 = norm, 2=snorm, 3=std, 4=sstd, 5=ged, 6=sged, 7=nig
	dist = model$modeldesc$distno
	if(modelinc[3]>0) mexdata = as.double(as.vector(mexdata)) else mexdata = double(1)
	
	ans = try( .C("starxfilters4s", model = as.integer(modelinc[1:60]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					x = as.double(data), res = double(T), z = double(T), 
					mexdata = mexdata, prob = as.double(probs), 
					constm = double(4*T), condm = double(4*T),  marx = double(T),
					T = as.integer(T), LHT = double(T), llh = double(1),  PACKAGE = "twinkle"), silent = TRUE )
	if(inherits(ans, "try-error")){
		if(arglist$pmode!=1){
			assign("twinkle_csol", 1, envir = starenv)
			assign("twinkle_filtermessage", ans, envir = starenv)
			if( trace > 0 ) cat(paste("\nstarfit-->warning: ", get("twinkle_filtermessage", starenv),"\n", sep=""))
			return(llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))))
		} else{
			return(llh = 1e10)
		}
	} else{
		if(arglist$pmode!=1){
			assign("twinkle_csol", 0, envir = starenv)
		}
	}
	z = ans$z
	epsx = ans$res
	llh = ans$llh
	
	if(is.finite(llh) && !is.na(llh) && !is.nan(llh)){
		if(arglist$pmode!=1) assign("star_llh", llh, envir = starenv)
	} else {
		if(arglist$pmode!=1) llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))) else llh = 1e10
	}
	
	# LHT = raw scores
	# LHT = -ans$LHT[(m):T]
	LHT = -ans$LHT
	ans = switch(returnType,
			llh = llh,
			LHT = LHT,
			all = list(llh = llh, res = epsx, z = z, LHT = LHT, probs = probs, condm = matrix(ans$condm, ncol = 4),
					constm = matrix(ans$constm, ncol=4), ipars = ipars, pmu = pmu))
	return(ans)
}
#--------------------------------------------------------------------------
# Dynamics (STARX-GARCH)
.stars1dLLH = function(pars, arglist)
{
	modelinc = arglist$model$modelinc
	if(modelinc[47]==1){
		probs = arglist$probs
		pmu = NULL
	} else{
		probs = rep(1, length(arglist$data))
	}
	data = arglist$data
	returnType = arglist$returnType
	starenv = arglist$starenv
	# rejoin fixed and pars
	model = arglist$model
	estidx = arglist$estidx
	idx = model$pidx
	ipars = arglist$ipars	
	ipars[estidx, 1] = pars
	if(arglist$transform==1){
		iidx = which(ipars[,4]==1 & ipars[,7]==1)
		if(length(iidx)>0){
			ipars[iidx,1] = logtransform(ipars[iidx,1], ipars[iidx,5], ipars[iidx,6])
		}
	}
	trace = arglist$trace
	T = length(data)
	mexdata = arglist$mexdata
	fit.control = arglist$fit.control	
	m = max(modelinc[32:33])
	distribution = model$modeldesc$distribution
	dist = model$modeldesc$distno
	
	sres = try( .C("starxfilters1s", model = as.integer(modelinc[1:60]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					x = as.double(data), res = double(T), z = double(T), 
					mexdata = mexdata, prob = as.double(probs), 
					constm = double(T), condm = double(T), marx = double(T),
					T = as.integer(T), LHT = double(T), llh = double(1),  PACKAGE = "twinkle"), silent = TRUE )
	if(inherits(sres, "try-error")){
		if(arglist$pmode!=1){
			assign("twinkle_csol", 1, envir = starenv)
			assign("twinkle_filtermessage", sres, envir = starenv)
			if( trace > 0 ) cat(paste("\nstarfit-->warning: ", get("twinkle_filtermessage", starenv),"\n", sep=""))
			return(llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))))
		} else{
			return(llh = 1e10)
		}
	} else{
		if(arglist$pmode!=1){
			assign("twinkle_csol", 0, envir = starenv)
		}
	}
	res = sres$res
	
	sumalpha = sum(ipars[idx["alpha",1]:idx["alpha",2],1])
	sumbeta  = sum(ipars[idx["beta",1]:idx["beta",2],1])
	lst = switch(modelinc[50],
			{
				list(persist = sumalpha+sumbeta, nres = double(1), meanz = double(1))
			},
			{
				persist = sumalpha + sumbeta + sum(apply(as.data.frame(ipars[idx["gamma",1]:idx["gamma",2],1]), 1 , FUN=function(x) 
									x * pneg(ipars[idx["ghlambda",1],1], ipars[idx["shape",1],1], ipars[idx["skew",1],1], distribution)))
				list(persist = persist, nres = double(T), meanz = double(1))
			},
			{
				meanz = as.double(egarchKappa(ipars[idx["ghlambda",1],1], ipars[idx["shape",1],1], ipars[idx["skew",1],1], distribution))
				list(persist = sumbeta, nres = double(1), meanz = meanz)
			})
	persist = lst$persist
	res[is.na(res) | !is.finite(res) | is.nan(res)] = 0

	mvar = ifelse(arglist$recinit$type==1, mean(res[1:arglist$recinit$n]*res[1:arglist$recinit$n]), backcastv(res, T, arglist$recinit$n))
	if(modelinc[39]>0) {
		mv = sum(apply(matrix(vexdata, ncol = modelinc[39]), 2, "mean")*ipars[idx["vxreg",1]:idx["vxreg",2],1])
	} else{
		mv = 0
	}
	hEst = mvar
	# variance targeting
	if(modelinc[50]==3){
		if(modelinc[31]==0){
			mvar2 = ifelse(modelinc[45]>=0, modelinc[45], mvar)
			ipars[idx["omega",1],1] = log(mvar2) * max(1 - persist, 0.001) - mv
			assign("omega", ipars[idx["omega",1],1], starenv)
		}
		if(is.na(hEst) | !is.finite(hEst) | is.nan(hEst)) hEst = var(data)
	} else {
		if(modelinc[31]>0){
			ipars[idx["omega",1],1] = max(eps, ipars[idx["omega",1],1]) 
		} else{
			mvar2 = ifelse(modelinc[45]>=0, modelinc[45], mvar)
			ipars[idx["omega",1],1] = mvar2 * (1 - persist) - mv
			assign("omega", ipars[idx["omega",1],1], starenv)
		}
		if(is.na(hEst) | !is.finite(hEst) | is.nan(hEst)) hEst = var(data) * (1 - persist)
	}
	# if we have external regressors in variance equation we cannot have
	# stationarity checks in likelihood. Stationarity conditions not valid for
	# solnp solver which implements them internally as constraints.
	if(fit.control$stationarity == 1 && modelinc[39] == 0){
		if(!is.na(persist) && persist >= 1){
			if(arglist$pmode!=1){
				return(llh = get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv))))
			} else{
				return(llh = 1e10)
			}
		}
	}
	if(modelinc[39]>0) vexdata = as.double(as.vector(vexdata)) else vexdata = double(1)
	
	ans = try( .C("garchfilterC", model = as.integer(modelinc[1:60]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					prob = double(1), res = as.double(res), e = double(T), z = double(T), vexdata = vexdata, 
					hEst = as.double(hEst), h = double(T), nres = lst$nres, 
					meanz = lst$meanz, m = as.integer(m), T = as.integer(T), 
					LHT = double(T), llh = double(1),  PACKAGE = "twinkle"), silent = TRUE )
	if(inherits(ans, "try-error")){
		if(arglist$pmode!=1){
			assign("twinkle_csol", 1, envir = starenv)
			assign("twinkle_filtermessage", ans, envir = starenv)
			if( trace > 0 ) cat(paste("\nstarfit-->warning: ", get("twinkle_filtermessage", starenv),"\n", sep=""))
			return(llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))))
		} else{
			return(llh = 1e10)
		}
	} else{
		if(arglist$pmode!=1){
			assign("twinkle_csol", 0, envir = starenv)
		}
	}
	h = ans$h
	z = ans$z
	epsx = ans$res
	llh = ans$llh
	
	if(is.finite(llh) && !is.na(llh) && !is.nan(llh)){
		if(arglist$pmode!=1) assign("star_llh", llh, envir = starenv)
	} else {
		if(arglist$pmode!=1) llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))) else llh = 1e10
	}
	LHT = -ans$LHT
	ans = switch(returnType,
			llh = llh,
			LHT = LHT,
			all = list(llh = llh, h = h, res = epsx, z = z, LHT = LHT, 
					probs = probs, condm = matrix(sres$condm, ncol=1), constm = matrix(sres$constm, ncol=1), ipars = ipars, pmu = pmu))
	return(ans)
}

.stars2dLLH = function(pars, arglist)
{
	modelinc = arglist$model$modelinc
	data = arglist$data
	returnType = arglist$returnType
	starenv = arglist$starenv
	# rejoin fixed and pars
	model = arglist$model
	estidx = arglist$estidx
	idx = model$pidx
	ipars = arglist$ipars	
	ipars[estidx, 1] = pars
	if(arglist$transform==1){
		iidx = which(ipars[,4]==1 & ipars[,7]==1)
		if(length(iidx)>0){
			ipars[iidx,1] = logtransform(ipars[iidx,1], ipars[iidx,5], ipars[iidx,6])
		}
	}
	# Transform to positive
	ipars["s1.gamma",1] = abs(ipars["s1.gamma",1])
	
	if(modelinc[47]==1){
		probs = arglist$probs
		pmu = NULL
	} else{
		xprobs = dstar2(ipars, arglist)
		probs = xprobs$probs
		pmu = xprobs$pmu
	}
	trace = arglist$trace
	T = length(data)
	mexdata = arglist$mexdata
	fit.control = arglist$fit.control	
	m = max(modelinc[32:33])
	distribution = model$modeldesc$distribution
	dist = model$modeldesc$distno
	if(modelinc[3]>0) mexdata = as.double(as.vector(mexdata)) else mexdata = double(1)
	
	sres = try( .C("starxfilters2sx", model = as.integer(modelinc[1:60]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					x = as.double(data), res = double(T), z = double(T), 
					mexdata = mexdata, prob = as.double(probs), 
					constm = double(2*T), condm = double(2*T), marx = double(T),
					T = as.integer(T), PACKAGE = "twinkle"), silent = TRUE )
	if(inherits(sres, "try-error")){
		if(arglist$pmode!=1){
			assign("twinkle_csol", 1, envir = starenv)
			assign("twinkle_filtermessage", sres, envir = starenv)
			if( trace > 0 ) cat(paste("\nstarfit-->warning: ", get("twinkle_filtermessage", starenv),"\n", sep=""))
			return(llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))))
		} else{
			return(llh = 1e10)
		}
	} else{
		if(arglist$pmode!=1){
			assign("twinkle_csol", 0, envir = starenv)
		}
	}
	res = sres$res
	
	sumalpha = sum(ipars[idx["alpha",1]:idx["alpha",2],1])
	sumbeta  = sum(ipars[idx["beta",1]:idx["beta",2],1])
	lst = switch(modelinc[50],
			{
				list(persist = sumalpha+sumbeta, nres = double(1), meanz = double(1))
			},
			{
				persist = sumalpha + sumbeta + sum(apply(as.data.frame(ipars[idx["gamma",1]:idx["gamma",2],1]), 1 , FUN=function(x) 
									x * pneg(ipars[idx["ghlambda",1],1], ipars[idx["shape",1],1], ipars[idx["skew",1],1], distribution)))
				list(persist = persist, nres = double(T), meanz = double(1))
			},
			{
				meanz = as.double(egarchKappa(ipars[idx["ghlambda",1],1], ipars[idx["shape",1],1], ipars[idx["skew",1],1], distribution))
				list(persist = sumbeta, nres = double(1), meanz = meanz)
			})
	persist = lst$persist
	res[is.na(res) | !is.finite(res) | is.nan(res)] = 0
	
	mvar = ifelse(arglist$recinit$type==1, mean(res[1:arglist$recinit$n]*res[1:arglist$recinit$n]), backcastv(res, T, arglist$recinit$n))
	if(modelinc[39]>0) {
		mv = sum(apply(matrix(vexdata, ncol = modelinc[39]), 2, "mean")*ipars[idx["vxreg",1]:idx["vxreg",2],1])
	} else{
		mv = 0
	}
	hEst = mvar
	# variance targeting
	if(modelinc[50]==3){
		if(modelinc[31]==0){
			mvar2 = ifelse(modelinc[45]>=0, modelinc[45], mvar)
			ipars[idx["omega",1],1] = log(mvar2) * max(1 - persist, 0.001) - mv
			assign("omega", ipars[idx["omega",1],1], starenv)
		}
		if(is.na(hEst) | !is.finite(hEst) | is.nan(hEst)) hEst = var(data)
	} else {
		if(modelinc[31]>0){
			ipars[idx["omega",1],1] = max(eps, ipars[idx["omega",1],1]) 
		} else{
			mvar2 = ifelse(modelinc[45]>=0, modelinc[45], mvar)
			ipars[idx["omega",1],1] = mvar2 * (1 - persist) - mv
			assign("omega", ipars[idx["omega",1],1], starenv)
		}
		if(is.na(hEst) | !is.finite(hEst) | is.nan(hEst)) hEst = var(data) * (1 - persist)
	}
	# if we have external regressors in variance equation we cannot have
	# stationarity checks in likelihood. Stationarity conditions not valid for
	# solnp solver which implements them internally as constraints.
	if(fit.control$stationarity == 1 && modelinc[39] == 0){
		if(!is.na(persist) && persist >= 1){
			if(arglist$pmode!=1){
				return(llh = get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv))))
			} else{
				return(llh = 1e10)
			}
		}
	}
	if(modelinc[39]>0) vexdata = as.double(as.vector(vexdata)) else vexdata = double(1)
	
	ans = try( .C("garchfilterC", model = as.integer(modelinc[1:60]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					prob = double(1), res = as.double(res), e = double(T), 
					z = double(T), vexdata = vexdata, 
					hEst = as.double(hEst), h = double(T), nres = lst$nres, 
					meanz = lst$meanz, m = as.integer(m), T = as.integer(T), 
					LHT = double(T), llh = double(1),  PACKAGE = "twinkle"), silent = TRUE )
	if(inherits(ans, "try-error")){
		if(arglist$pmode!=1){
			assign("twinkle_csol", 1, envir = starenv)
			assign("twinkle_filtermessage", ans, envir = starenv)
			if( trace > 0 ) cat(paste("\nstarfit-->warning: ", get("twinkle_filtermessage", starenv),"\n", sep=""))
			return(llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))))
		} else{
			return(llh = 1e10)
		}
	} else{
		if(arglist$pmode!=1){
			assign("twinkle_csol", 0, envir = starenv)
		}
	}
	h = ans$h
	z = ans$z
	epsx = ans$res
	llh = ans$llh
	
	if(is.finite(llh) && !is.na(llh) && !is.nan(llh)){
		if(arglist$pmode!=1) assign("star_llh", llh, envir = starenv)
	} else {
		if(arglist$pmode!=1) llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))) else llh = 1e10
	}
	LHT = -ans$LHT
	ans = switch(returnType,
			llh = llh,
			LHT = LHT,
			all = list(llh = llh, h = h, res = epsx, z = z, LHT = LHT, 
					probs = probs, condm = matrix(sres$condm, ncol=2), constm = matrix(sres$constm, ncol=2), ipars = ipars, pmu = pmu))
	return(ans)
}


.stars2dLLH2mix = function(pars, arglist)
{	
	modelinc = arglist$model$modelinc
	data = arglist$data
	returnType = arglist$returnType
	starenv = arglist$starenv
	# rejoin fixed and pars
	model = arglist$model
	estidx = arglist$estidx
	idx = model$pidx
	ipars = arglist$ipars
	ipars[estidx, 1] = pars
	if(arglist$transform==1){
		iidx = which(ipars[,4]==1 & ipars[,7]==1)
		if(length(iidx)>0){
			ipars[iidx,1] = logtransform(ipars[iidx,1], ipars[iidx,5], ipars[iidx,6])
		}
	}
	# Transform to positive
	ipars["s1.gamma",1] = abs(ipars["s1.gamma",1])
	
	if(modelinc[47]==1){
		probs = arglist$probs
		pmu = NULL
	} else{
		xprobs = dstar2(ipars, arglist)
		probs = xprobs$probs
		pmu = xprobs$pmu
	}
	trace = arglist$trace
	T = length(data)
	mexdata = arglist$mexdata
	fit.control = arglist$fit.control	
	m = 0
	distribution = model$modeldesc$distribution
	dist = model$modeldesc$distno
	if(modelinc[3]>0) mexdata = as.double(as.vector(mexdata)) else mexdata = double(1)
	
	sres = try( .C("starxfilters2sx", model = as.integer(modelinc[1:60]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					x = as.double(data), res = double(T), z = double(T), 
					mexdata = mexdata, prob = as.double(probs), 
					constm = double(2*T), condm = double(2*T), marx = double(T),
					T = as.integer(T), PACKAGE = "twinkle"), silent = TRUE )
	if(inherits(sres, "try-error")){
		if(arglist$pmode!=1){
			assign("twinkle_csol", 1, envir = starenv)
			assign("twinkle_filtermessage", sres, envir = starenv)
			if( trace > 0 ) cat(paste("\nstarfit-->warning: ", get("twinkle_filtermessage", starenv),"\n", sep=""))
			return(llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))))
		} else{
			return(llh = 1e10)
		}
	} else{
		if(arglist$pmode!=1){
			assign("twinkle_csol", 0, envir = starenv)
		}
	}
	res = sres$res	
	res[is.na(res) | !is.finite(res) | is.nan(res)] = 0		
	ans = try( .C("mixturefilterC", model = as.integer(modelinc[1:60]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					prob = as.double(probs), res = as.double(res), z = double(T), 
					h = double(T), T = as.integer(T), 
					LHT = double(T), llh = double(1),  PACKAGE = "twinkle"), silent = TRUE )
	if(inherits(ans, "try-error")){
		if(arglist$pmode!=1){
			assign("twinkle_csol", 1, envir = starenv)
			assign("twinkle_filtermessage", ans, envir = starenv)
			if( trace > 0 ) cat(paste("\nstarfit-->warning: ", get("twinkle_filtermessage", starenv),"\n", sep=""))
			return(llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))))
		} else{
			return(llh = 1e10)
		}
	} else{
		if(arglist$pmode!=1){
			assign("twinkle_csol", 0, envir = starenv)
		}
	}
	h = ans$h
	z = ans$z
	epsx = ans$res
	llh = ans$llh
	if(is.finite(llh) && !is.na(llh) && !is.nan(llh)){
		if(arglist$pmode!=1) assign("star_llh", llh, envir = starenv)
	} else {
		if(arglist$pmode!=1) llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))) else llh = 1e10
	}
	LHT = -ans$LHT
	ans = switch(returnType,
			llh = llh,
			LHT = LHT,
			all = list(llh = llh, h = h, res = epsx, z = z, LHT = LHT, 
					probs = probs, condm = matrix(sres$condm, ncol=2), 
					constm = matrix(sres$constm, ncol=2), ipars = ipars, pmu = pmu))
	return(ans)
}

.stars3dLLH3mix = function(pars, arglist)
{	
	modelinc = arglist$model$modelinc
	data = arglist$data
	returnType = arglist$returnType
	starenv = arglist$starenv
	# rejoin fixed and pars
	model = arglist$model
	estidx = arglist$estidx
	idx = model$pidx
	ipars = arglist$ipars
	ipars[estidx, 1] = pars
	if(arglist$transform==1){
		iidx = which(ipars[,4]==1 & ipars[,7]==1)
		if(length(iidx)>0){
			ipars[iidx,1] = logtransform(ipars[iidx,1], ipars[iidx,5], ipars[iidx,6])
		}
	}
	# Transform to positive
	ipars["s1.gamma",1] = abs(ipars["s1.gamma",1])
	ipars["s2.gamma",1] = abs(ipars["s2.gamma",1])
	
	if(modelinc[47]==1){
		probs = arglist$probs
		pmu = NULL
	} else{
		xprobs = dstar3(ipars, arglist)
		probs = xprobs$probs
		pmu = xprobs$pmu
	}
	trace = arglist$trace
	T = length(data)
	mexdata = arglist$mexdata
	fit.control = arglist$fit.control	
	m = 0
	distribution = model$modeldesc$distribution
	dist = model$modeldesc$distno
	if(modelinc[3]>0) mexdata = as.double(as.vector(mexdata)) else mexdata = double(1)
	
	sres = try( .C("starxfilters3sx", model = as.integer(modelinc[1:60]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					x = as.double(data), res = double(T), z = double(T), 
					mexdata = mexdata, prob = as.double(probs), 
					constm = double(3*T), condm = double(3*T), marx = double(T),
					T = as.integer(T), PACKAGE = "twinkle"), silent = TRUE )
	if(inherits(sres, "try-error")){
		if(arglist$pmode!=1){
			assign("twinkle_csol", 1, envir = starenv)
			assign("twinkle_filtermessage", sres, envir = starenv)
			if( trace > 0 ) cat(paste("\nstarfit-->warning: ", get("twinkle_filtermessage", starenv),"\n", sep=""))
			return(llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))))
		} else{
			return(llh = 1e10)
		}
	} else{
		if(arglist$pmode!=1){
			assign("twinkle_csol", 0, envir = starenv)
		}
	}
	res = sres$res	
	res[is.na(res) | !is.finite(res) | is.nan(res)] = 0		
	ans = try( .C("mixturefilterC", model = as.integer(modelinc[1:60]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					prob = as.double(probs), res = as.double(res), z = double(T), 
					h = double(T), T = as.integer(T), 
					LHT = double(T), llh = double(1),  PACKAGE = "twinkle"), silent = TRUE )
	if(inherits(ans, "try-error")){
		if(arglist$pmode!=1){
			assign("twinkle_csol", 1, envir = starenv)
			assign("twinkle_filtermessage", ans, envir = starenv)
			if( trace > 0 ) cat(paste("\nstarfit-->warning: ", get("twinkle_filtermessage", starenv),"\n", sep=""))
			return(llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))))
		} else{
			return(llh = 1e10)
		}
	} else{
		if(arglist$pmode!=1){
			assign("twinkle_csol", 0, envir = starenv)
		}
	}
	h = ans$h
	z = ans$z
	epsx = ans$res
	llh = ans$llh
	if(is.finite(llh) && !is.na(llh) && !is.nan(llh)){
		if(arglist$pmode!=1) assign("star_llh", llh, envir = starenv)
	} else {
		if(arglist$pmode!=1) llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))) else llh = 1e10
	}
	LHT = -ans$LHT
	ans = switch(returnType,
			llh = llh,
			LHT = LHT,
			all = list(llh = llh, h = h, res = epsx, z = z, LHT = LHT, 
					probs = probs, condm = matrix(sres$condm, ncol=3), 
					constm = matrix(sres$constm, ncol=3), ipars = ipars, pmu = pmu))
	return(ans)
}

.stars4dLLH4mix = function(pars, arglist)
{	
	modelinc = arglist$model$modelinc
	data = arglist$data
	returnType = arglist$returnType
	starenv = arglist$starenv
	# rejoin fixed and pars
	model = arglist$model
	estidx = arglist$estidx
	idx = model$pidx
	ipars = arglist$ipars
	ipars[estidx, 1] = pars
	if(arglist$transform==1){
		iidx = which(ipars[,4]==1 & ipars[,7]==1)
		if(length(iidx)>0){
			ipars[iidx,1] = logtransform(ipars[iidx,1], ipars[iidx,5], ipars[iidx,6])
		}
	}
	# Transform to positive
	ipars["s1.gamma",1] = abs(ipars["s1.gamma",1])
	ipars["s2.gamma",1] = abs(ipars["s2.gamma",1])
	ipars["s3.gamma",1] = abs(ipars["s3.gamma",1])
	
	
	if(modelinc[47]==1){
		probs = arglist$probs
		pmu = NULL
	} else{
		xprobs = dstar4(ipars, arglist)
		probs = xprobs$probs
		pmu = xprobs$pmu
	}
	trace = arglist$trace
	T = length(data)
	mexdata = arglist$mexdata
	fit.control = arglist$fit.control	
	m = 0
	distribution = model$modeldesc$distribution
	dist = model$modeldesc$distno
	if(modelinc[3]>0) mexdata = as.double(as.vector(mexdata)) else mexdata = double(1)
	
	sres = try( .C("starxfilters4sx", model = as.integer(modelinc[1:60]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					x = as.double(data), res = double(T), z = double(T), 
					mexdata = mexdata, prob = as.double(probs), 
					constm = double(4*T), condm = double(4*T), marx = double(T),
					T = as.integer(T), PACKAGE = "twinkle"), silent = TRUE )
	if(inherits(sres, "try-error")){
		if(arglist$pmode!=1){
			assign("twinkle_csol", 1, envir = starenv)
			assign("twinkle_filtermessage", sres, envir = starenv)
			if( trace > 0 ) cat(paste("\nstarfit-->warning: ", get("twinkle_filtermessage", starenv),"\n", sep=""))
			return(llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))))
		} else{
			return(llh = 1e10)
		}
	} else{
		if(arglist$pmode!=1){
			assign("twinkle_csol", 0, envir = starenv)
		}
	}
	res = sres$res	
	res[is.na(res) | !is.finite(res) | is.nan(res)] = 0		
	ans = try( .C("mixturefilterC", model = as.integer(modelinc[1:60]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					prob = as.double(probs), res = as.double(res), z = double(T), 
					h = double(T), T = as.integer(T), 
					LHT = double(T), llh = double(1),  PACKAGE = "twinkle"), silent = TRUE )
	if(inherits(ans, "try-error")){
		if(arglist$pmode!=1){
			assign("twinkle_csol", 1, envir = starenv)
			assign("twinkle_filtermessage", ans, envir = starenv)
			if( trace > 0 ) cat(paste("\nstarfit-->warning: ", get("twinkle_filtermessage", starenv),"\n", sep=""))
			return(llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))))
		} else{
			return(llh = 1e10)
		}
	} else{
		if(arglist$pmode!=1){
			assign("twinkle_csol", 0, envir = starenv)
		}
	}
	h = ans$h
	z = ans$z
	epsx = ans$res
	llh = ans$llh
	if(is.finite(llh) && !is.na(llh) && !is.nan(llh)){
		if(arglist$pmode!=1) assign("star_llh", llh, envir = starenv)
	} else {
		if(arglist$pmode!=1) llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))) else llh = 1e10
	}
	LHT = -ans$LHT
	ans = switch(returnType,
			llh = llh,
			LHT = LHT,
			all = list(llh = llh, h = h, res = epsx, z = z, LHT = LHT, 
					probs = probs, condm = matrix(sres$condm, ncol=4), 
					constm = matrix(sres$constm, ncol=4), ipars = ipars, pmu = pmu))
	return(ans)
}

.stars3dLLH = function(pars, arglist)
{
	modelinc = arglist$model$modelinc
	data = arglist$data
	returnType = arglist$returnType
	starenv = arglist$starenv
	# rejoin fixed and pars
	model = arglist$model
	estidx = arglist$estidx
	idx = model$pidx
	ipars = arglist$ipars	
	ipars[estidx, 1] = pars
	if(arglist$transform==1){
		iidx = which(ipars[,4]==1 & ipars[,7]==1)
		if(length(iidx)>0){
			ipars[iidx,1] = logtransform(ipars[iidx,1], ipars[iidx,5], ipars[iidx,6])
		}
	}
	# Transform to positive
	ipars["s1.gamma",1] = abs(ipars["s1.gamma",1])
	ipars["s2.gamma",1] = abs(ipars["s2.gamma",1])
	
	if(modelinc[47]==1){
		probs = arglist$probs
		pmu = NULL
	} else{
		xprobs = dstar3(ipars, arglist)
		probs = xprobs$probs
		pmu = xprobs$pmu
	}
	trace = arglist$trace
	T = length(data)
	mexdata = arglist$mexdata
	fit.control = arglist$fit.control	
	m = max(modelinc[32:33])
	distribution = model$modeldesc$distribution
	dist = model$modeldesc$distno
	if(modelinc[3]>0) mexdata = as.double(as.vector(mexdata)) else mexdata = double(1)
	
	sres = try( .C("starxfilters3s", model = as.integer(modelinc[1:60]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					x = as.double(data), res = double(T), z = double(T), 
					mexdata = mexdata, prob = as.double(probs), 
					constm = double(3*T), condm = double(3*T), marx = double(T),
					T = as.integer(T), LHT = double(T), llh = double(1),  PACKAGE = "twinkle"), silent = TRUE )
	
	if(inherits(sres, "try-error")){
		if(arglist$pmode!=1){
			assign("twinkle_csol", 1, envir = starenv)
			assign("twinkle_filtermessage", sres, envir = starenv)
			if( trace > 0 ) cat(paste("\nstarfit-->warning: ", get("twinkle_filtermessage", starenv),"\n", sep=""))
			return(llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))))
		} else{
			return(llh = 1e10)
		}
	} else{
		if(arglist$pmode!=1){
			assign("twinkle_csol", 0, envir = starenv)
		}
	}
	res = sres$res
	
	sumalpha = sum(ipars[idx["alpha",1]:idx["alpha",2],1])
	sumbeta  = sum(ipars[idx["beta",1]:idx["beta",2],1])
	lst = switch(modelinc[50],
			{
				list(persist = sumalpha+sumbeta, nres = double(1), meanz = double(1))
			},
			{
				persist = sumalpha + sumbeta + sum(apply(as.data.frame(ipars[idx["gamma",1]:idx["gamma",2],1]), 1 , FUN=function(x) 
									x * pneg(ipars[idx["ghlambda",1],1], ipars[idx["shape",1],1], ipars[idx["skew",1],1], distribution)))
				list(persist = persist, nres = double(T), meanz = double(1))
			},
			{
				meanz = as.double(egarchKappa(ipars[idx["ghlambda",1],1], ipars[idx["shape",1],1], ipars[idx["skew",1],1], distribution))
				list(persist = sumbeta, nres = double(1), meanz = meanz)
			})
	persist = lst$persist
	res[is.na(res) | !is.finite(res) | is.nan(res)] = 0
	
	mvar = ifelse(arglist$recinit$type==1, mean(res[1:arglist$recinit$n]*res[1:arglist$recinit$n]), backcastv(res, T, arglist$recinit$n))
	if(modelinc[39]>0) {
		mv = sum(apply(matrix(vexdata, ncol = modelinc[39]), 2, "mean")*ipars[idx["vxreg",1]:idx["vxreg",2],1])
	} else{
		mv = 0
	}
	hEst = mvar
	# variance targeting
	if(modelinc[50]==3){
		if(modelinc[31]==0){
			mvar2 = ifelse(modelinc[45]>=0, modelinc[45], mvar)
			ipars[idx["omega",1],1] = log(mvar2) * max(1 - persist, 0.001) - mv
			assign("omega", ipars[idx["omega",1],1], starenv)
		}
		if(is.na(hEst) | !is.finite(hEst) | is.nan(hEst)) hEst = var(data)
	} else {
		if(modelinc[31]>0){
			ipars[idx["omega",1],1] = max(eps, ipars[idx["omega",1],1]) 
		} else{
			mvar2 = ifelse(modelinc[45]>=0, modelinc[45], mvar)
			ipars[idx["omega",1],1] = mvar2 * (1 - persist) - mv
			assign("omega", ipars[idx["omega",1],1], starenv)
		}
		if(is.na(hEst) | !is.finite(hEst) | is.nan(hEst)) hEst = var(data) * (1 - persist)
	}
	# if we have external regressors in variance equation we cannot have
	# stationarity checks in likelihood. Stationarity conditions not valid for
	# solnp solver which implements them internally as constraints.
	if(fit.control$stationarity == 1 && modelinc[39] == 0){
		if(!is.na(persist) && persist >= 1){
			if(arglist$pmode!=1){
				return(llh = get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv))))
			} else{
				return(llh = 1e10)
			}
		}
	}
	if(modelinc[39]>0) vexdata = as.double(as.vector(vexdata)) else vexdata = double(1)
	
	ans = try( .C("garchfilterC", model = as.integer(modelinc[1:60]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					prob = double(1), res = as.double(res), e = double(T), 
					z = double(T), vexdata = vexdata, 
					hEst = as.double(hEst), h = double(T), nres = lst$nres, 
					meanz = lst$meanz, m = as.integer(m), T = as.integer(T), 
					LHT = double(T), llh = double(1),  PACKAGE = "twinkle"), silent = TRUE )
	if(inherits(ans, "try-error")){
		if(arglist$pmode!=1){
			assign("twinkle_csol", 1, envir = starenv)
			assign("twinkle_filtermessage", ans, envir = starenv)
			if( trace > 0 ) cat(paste("\nstarfit-->warning: ", get("twinkle_filtermessage", starenv),"\n", sep=""))
			return(llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))))
		} else{
			return(llh = 1e10)
		}
	} else{
		if(arglist$pmode!=1){
			assign("twinkle_csol", 0, envir = starenv)
		}
	}
	h = ans$h
	z = ans$z
	epsx = ans$res
	llh = ans$llh
	
	if(is.finite(llh) && !is.na(llh) && !is.nan(llh)){
		if(arglist$pmode!=1) assign("star_llh", llh, envir = starenv)
	} else {
		if(arglist$pmode!=1) llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))) else llh = 1e10
	}
	LHT = -ans$LHT
	ans = switch(returnType,
			llh = llh,
			LHT = LHT,
			all = list(llh = llh, h = h, res = epsx, z = z, LHT = LHT, 
					probs = probs, condm = matrix(sres$condm, ncol=3), constm = matrix(sres$constm, ncol=3), ipars = ipars, pmu = pmu))
	return(ans)
}

.stars4dLLH = function(pars, arglist)
{
	modelinc = arglist$model$modelinc
	data = arglist$data
	returnType = arglist$returnType
	starenv = arglist$starenv
	# rejoin fixed and pars
	model = arglist$model
	estidx = arglist$estidx
	idx = model$pidx
	ipars = arglist$ipars	
	ipars[estidx, 1] = pars
	if(arglist$transform==1){
		iidx = which(ipars[,4]==1 & ipars[,7]==1)
		if(length(iidx)>0){
			ipars[iidx,1] = logtransform(ipars[iidx,1], ipars[iidx,5], ipars[iidx,6])
		}
	}
	# Transform to positive
	ipars["s1.gamma",1] = abs(ipars["s1.gamma",1])
	ipars["s2.gamma",1] = abs(ipars["s2.gamma",1])
	ipars["s3.gamma",1] = abs(ipars["s3.gamma",1])
	
	if(modelinc[47]==1){
		probs = arglist$probs
		pmu = NULL
	} else{
		xprobs = dstar4(ipars, arglist)
		probs = xprobs$probs
		pmu = xprobs$pmu
	}
	trace = arglist$trace
	T = length(data)
	mexdata = arglist$mexdata
	fit.control = arglist$fit.control	
	m = max(modelinc[32:33])
	distribution = model$modeldesc$distribution
	dist = model$modeldesc$distno
	if(modelinc[3]>0) mexdata = as.double(as.vector(mexdata)) else mexdata = double(1)
	
	sres = try( .C("starxfilters4s", model = as.integer(modelinc[1:60]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					x = as.double(data), res = double(T), z = double(T), 
					mexdata = mexdata, prob = as.double(probs), 
					constm = double(4*T), condm = double(4*T), marx = double(T),
					T = as.integer(T), LHT = double(T), llh = double(1),  PACKAGE = "twinkle"), silent = TRUE )
	
	if(inherits(sres, "try-error")){
		if(arglist$pmode!=1){
			assign("twinkle_csol", 1, envir = starenv)
			assign("twinkle_filtermessage", sres, envir = starenv)
			if( trace > 0 ) cat(paste("\nstarfit-->warning: ", get("twinkle_filtermessage", starenv),"\n", sep=""))
			return(llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))))
		} else{
			return(llh = 1e10)
		}
	} else{
		if(arglist$pmode!=1){
			assign("twinkle_csol", 0, envir = starenv)
		}
	}
	res = sres$res
	
	sumalpha = sum(ipars[idx["alpha",1]:idx["alpha",2],1])
	sumbeta  = sum(ipars[idx["beta",1]:idx["beta",2],1])
	lst = switch(modelinc[50],
			{
				list(persist = sumalpha+sumbeta, nres = double(1), meanz = double(1))
			},
			{
				persist = sumalpha + sumbeta + sum(apply(as.data.frame(ipars[idx["gamma",1]:idx["gamma",2],1]), 1 , FUN=function(x) 
									x * pneg(ipars[idx["ghlambda",1],1], ipars[idx["shape",1],1], ipars[idx["skew",1],1], distribution)))
				list(persist = persist, nres = double(T), meanz = double(1))
			},
			{
				meanz = as.double(egarchKappa(ipars[idx["ghlambda",1],1], ipars[idx["shape",1],1], ipars[idx["skew",1],1], distribution))
				list(persist = sumbeta, nres = double(1), meanz = meanz)
			})
	persist = lst$persist
	res[is.na(res) | !is.finite(res) | is.nan(res)] = 0
	
	mvar = ifelse(arglist$recinit$type==1, mean(res[1:arglist$recinit$n]*res[1:arglist$recinit$n]), backcastv(res, T, arglist$recinit$n))
	if(modelinc[39]>0) {
		mv = sum(apply(matrix(vexdata, ncol = modelinc[39]), 2, "mean")*ipars[idx["vxreg",1]:idx["vxreg",2],1])
	} else{
		mv = 0
	}
	hEst = mvar
	# variance targeting
	if(modelinc[50]==3){
		if(modelinc[31]==0){
			mvar2 = ifelse(modelinc[45]>=0, modelinc[45], mvar)
			ipars[idx["omega",1],1] = log(mvar2) * max(1 - persist, 0.001) - mv
			assign("omega", ipars[idx["omega",1],1], starenv)
		}
		if(is.na(hEst) | !is.finite(hEst) | is.nan(hEst)) hEst = var(data)
	} else {
		if(modelinc[31]>0){
			ipars[idx["omega",1],1] = max(eps, ipars[idx["omega",1],1]) 
		} else{
			mvar2 = ifelse(modelinc[45]>=0, modelinc[45], mvar)
			ipars[idx["omega",1],1] = mvar2 * (1 - persist) - mv
			assign("omega", ipars[idx["omega",1],1], starenv)
		}
		if(is.na(hEst) | !is.finite(hEst) | is.nan(hEst)) hEst = var(data) * (1 - persist)
	}
	# if we have external regressors in variance equation we cannot have
	# stationarity checks in likelihood. Stationarity conditions not valid for
	# solnp solver which implements them internally as constraints.
	if(fit.control$stationarity == 1 && modelinc[39] == 0){
		if(!is.na(persist) && persist >= 1){
			if(arglist$pmode!=1){
				return(llh = get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv))))
			} else{
				return(llh = 1e10)
			}
		}
	}
	if(modelinc[39]>0) vexdata = as.double(as.vector(vexdata)) else vexdata = double(1)
	
	ans = try( .C("garchfilterC", model = as.integer(modelinc[1:60]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					prob = double(1), res = as.double(res), e = double(T), 
					z = double(T), vexdata = vexdata, 
					hEst = as.double(hEst), h = double(T), nres = lst$nres, 
					meanz = lst$meanz, m = as.integer(m), T = as.integer(T), 
					LHT = double(T), llh = double(1),  PACKAGE = "twinkle"), silent = TRUE )
	if(inherits(ans, "try-error")){
		if(arglist$pmode!=1){
			assign("twinkle_csol", 1, envir = starenv)
			assign("twinkle_filtermessage", ans, envir = starenv)
			if( trace > 0 ) cat(paste("\nstarfit-->warning: ", get("twinkle_filtermessage", starenv),"\n", sep=""))
			return(llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))))
		} else{
			return(llh = 1e10)
		}
	} else{
		if(arglist$pmode!=1){
			assign("twinkle_csol", 0, envir = starenv)
		}
	}
	h = ans$h
	z = ans$z
	epsx = ans$res
	llh = ans$llh
	
	if(is.finite(llh) && !is.na(llh) && !is.nan(llh)){
		if(arglist$pmode!=1) assign("star_llh", llh, envir = starenv)
	} else {
		if(arglist$pmode!=1) llh = (get("star_llh", starenv) + 0.1*(abs(get("star_llh", starenv)))) else llh = 1e10
	}
	LHT = -ans$LHT
	ans = switch(returnType,
			llh = llh,
			LHT = LHT,
			all = list(llh = llh, h = h, res = epsx, z = z, LHT = LHT, 
					probs = probs, condm = matrix(sres$condm, ncol=4), constm = matrix(sres$constm, ncol=4), ipars = ipars, pmu = pmu))
	return(ans)
}

dstar2 = function(ipars, arglist)
{
	probs = arglist$probs
	model = arglist$model
	estidx = arglist$estidx
	idx = model$pidx
	XL = arglist$XL
	N = arglist$N
	modelinc = model$modelinc
	if(modelinc[21]>0){
		beta = ipars[idx["s1.beta",1],1]
	} else{
		beta = 0
	}
	gamma = ipars[idx["s1.gamma",1],1]
	alpha = ipars[idx["s1.alpha",1]:idx["s1.alpha",2],1]
	cnst = ipars[idx["s1.c",1],1]
	initp = gamma*(sum(alpha * colMeans(XL[1:N,,drop=FALSE]))-cnst) /(1-beta)
	pmu = gamma*(as.numeric(XL%*%alpha)-cnst)
	if(modelinc[21]>0) pmu = .recfilter(as.double(pmu), as.double(beta), init = as.double(initp))
	probs[,1] = 1/(1+exp(-pmu))
	probs[,2] = 1 - probs[,1]
	return(list(probs = probs, pmu = matrix(pmu, ncol=1), initp = initp))
}

dstar3 = function(ipars, arglist)
{
	probs = arglist$probs
	model = arglist$model
	estidx = arglist$estidx
	idx = model$pidx
	XL = arglist$XL
	N = arglist$N
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
	
	initp1 = gamma1*(sum(alpha1 * colMeans(XL[1:N,,drop=FALSE]))-cnst1)/(1-beta1)
	initp2 = gamma2*(sum(alpha2 * colMeans(XL[1:N,,drop=FALSE]))-cnst2)/(1-beta2)
	pmu1 = gamma1*(as.numeric(XL%*%alpha1)-cnst1)
	pmu2 = gamma2*(as.numeric(XL%*%alpha2)-cnst2)
	if(modelinc[21]>0) pmu1 = .recfilter(as.double(pmu1), as.double(beta1), init = as.double(initp1))
	if(modelinc[25]>0) pmu2 = .recfilter(as.double(pmu2), as.double(beta2), init = as.double(initp2))
	expmu1 = pmin(exp(100), exp(pmu1))
	expmu2 = pmin(exp(100), exp(pmu2))
	probs[,1] = expmu1/(1+expmu1+expmu2)
	probs[,2] = expmu2/(1+expmu1+expmu2)
	probs[,3] = 1-probs[,1]-probs[,2]
	return(list(probs = probs, pmu = cbind(pmu1, pmu2), initp = c(initp1, initp2)))
}

dstar4 = function(ipars, arglist)
{
	probs = arglist$probs
	model = arglist$model
	estidx = arglist$estidx
	idx = model$pidx
	XL = arglist$XL
	N = arglist$N
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
	
	initp1 = gamma1*(sum(alpha1 * colMeans(XL[1:N,,drop=FALSE]))-cnst1) /(1-beta1)
	initp2 = gamma2*(sum(alpha2 * colMeans(XL[1:N,,drop=FALSE]))-cnst2) /(1-beta2)
	initp3 = gamma3*(sum(alpha3 * colMeans(XL[1:N,,drop=FALSE]))-cnst3) /(1-beta3)
	
	pmu1 = gamma*(as.numeric(XL%*%alpha1)-cnst1)
	pmu2 = gamma*(as.numeric(XL%*%alpha2)-cnst2)
	pmu3 = gamma*(as.numeric(XL%*%alpha3)-cnst3)
	
	if(modelinc[21]>0) pmu1 = .recfilter(as.double(pmu1), as.double(beta1), init = as.double(initp1))
	if(modelinc[25]>0) pmu2 = .recfilter(as.double(pmu2), as.double(beta2), init = as.double(initp2))
	if(modelinc[29]>0) pmu3 = .recfilter(as.double(pmu3), as.double(beta3), init = as.double(initp3))
	
	expmu1 = pmin(exp(60), exp(pmu1))
	expmu2 = pmin(exp(60), exp(pmu2))
	expmu3 = pmin(exp(60), exp(pmu3))
	
	probs[,1] = expmu1/(1+expmu1+expmu2+expmu3)
	probs[,2] = expmu2/(1+expmu1+expmu2+expmu3)
	probs[,3] = expmu3/(1+expmu1+expmu2+expmu3)
	probs[,4] = 1-probs[,1]-probs[,2]-probs[,3]
	return(list(probs = probs, pmu = cbind(pmu1, pmu2, pmu3), initp = c(initp1, initp2, initp3)))
}


safefun = function(fun, p)
{
	ans = try(do.call(fun, list(p)), silent=TRUE)
	if(inherits(ans, 'try-error')){
		warning("\nstarfit-->warning: custom fun returned error...not applying fun on unconstrained probability dynamics..")
		ans = p
	} else{
		if(length(ans!=length(p))) 	warning("\nstarfit-->warning: custom fun returned wrong length...not applying fun on unconstrained probability dynamics..")
		ans = p
	}
	if(any(is.na(ans))){
		warning("\nstarfit-->warning: custom fun returned NAs...using na.approx to fill.")
		ans = na.approx(ans)
	}
	if(any(is.na(ans))){
		warning("\nstarfit-->warning: custom fun returned NAs...na.approx could not to fill....not applying fun to unconstrained probability dynamics.")
		ans = p
	}
	if(any(!is.finite(ans))){
		warning("\nstarfit-->warning: custom fun returned non-finite values....not applying fun to unconstrained probability dynamics.")
		ans = p
	}
	return(ans)
}


.recfilter = function(x, filter, init = NULL) 
{
	sides = 2L
	storage.mode(x) <- "double"
	n <- as.integer(length(x))
	if (is.na(n)) stop("invalid value of nrow(x)", domain = NA)
	nser <- 1
	filter <- as.double(filter)
	nfilt <- as.integer(length(filter))
	if(is.na(n)) stop("invalid value of length(filter)", domain = NA)
	if(any(is.na(filter))) stop("missing values in 'filter'")
	ni <- length(init)
	dim(init) <- c(nfilt, nser)
	ind <- seq_len(nfilt)
	y <- try(.C("rfilter", x = x, filter = filter, out = c(rev(init[, 1L]), double(n)), n = as.integer(c(length(x), length(filter))), PACKAGE="twinkle"), silent=TRUE)
	return(y$out[-ind])
}