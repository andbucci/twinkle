#################################################################################
##
##   R package twinkle by Alexios Ghalanos Copyright (C) 2013.
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

#------------------------------------------------------------------------------------
# Model Specification
#------------------------------------------------------------------------------------
starspec = function(
		mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(1,1), 
				maOrder = c(0,0), matype = "linear", statevar = c("y","s"), s = NULL, 
				ylags = 1, statear = FALSE, yfun = NULL, xreg = NULL,
				transform = "logis"), 
		variance.model = list(dynamic = FALSE, model = "sGARCH", 
				garchOrder = c(1, 1), submodel = NULL, vreg = NULL,
				variance.targeting = FALSE), 
		distribution.model = "norm", 
		start.pars = list(), fixed.pars = list(), fixed.prob = NULL, ...)
{
	UseMethod("starspec")
}


# modelinc[1:4] state 1 (mu, phi, xi, ma)
# modelinc[5:8] state 2 (mu, phi, xi, ma)
# modelinc[9:12] state 3 (mu, phi, xi, ma)
# modelinc[13:16] state 3 (mu, phi, xi, ma)
# modelinc[17] = linear MA (>0 = linear, 0 = state)
# modelinc[18:21] star (gamma, c, alpha, beta) state 1:2
# modelinc[22:25] star (gamma, c, alpha, beta) state 3
# modelinc[26:29] star (gamma, c, alpha, beta) state 4
# modelinc[30] sigma (if 1 then no GARCH else GARCH)
# modelinc[30] = sigma1 if mixture
# modelinc[31] = sigma2 if mixture
# modelinc[32] = sigma3 if mixture
# modelinc[33] = sigma4 if mixture
# modelinc[31:40] GARCH (omega, alpha, beta, gamma, eta1, eta2, delta, lambda, vxreg, xi)
# modelinc[41:43] Distribution (skew, shape, ghlambda)
# modelinc[44] distribution c("norm", "snorm", "std", "sstd","ged", "sged", "nig", "ghyp", "jsu", "ghst")
# modelinc[45] variance targeting (FALSE/TRUE)
# modelinc[46] states (max 4)
# modelinc[47] fixed.prob (FALSE/TRUE)
# modelinc[48] yfun (TRUE/FALSE)
# modelinc[49] type (1=y, 2=s)
# modelinc[50] garchmodel (1=sGARCH, 2=gjrGARCH, 3=eGARCH, 4=mixture)
# if(modelinc[50]==0) then NOT dynamic
# modelinc[51] AR type (1: \phi*(y_t-\mu-t), 0: \phi*(y_t) )
# modelinc[52] = lag type (1 = logis, 2 = exp)
.xstarspec = function( 
		mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(1,1), 
				maOrder = c(0,0), matype = "linear", statevar = c("y","s"), s = NULL, 
				ylags = 1, statear = FALSE, yfun = NULL, xreg = NULL, 
				transform = "logis"), 
		variance.model = list(dynamic = FALSE, model = "sGARCH", 
				garchOrder = c(1, 1), submodel = NULL, vreg = NULL,
				variance.targeting = FALSE), 
		distribution.model = "norm", 
		start.pars = list(), fixed.pars = list(), fixed.prob = NULL, ...)
{
	# external.regressors is an xts matrix aligned to the dates of the
	# data matrix (which is passed in the starfit routine), pre-lagged.
	# "lags" is a vector:
	# for the y and v variables this represents the lags to use
	# yfun is a user supplied function which, given y[1:(t-1)] (or x[1:(t-1),]), 
	# returns a transformed value(s) for use at time t.
	intmodel = 0
	if(is.null(variance.model$dynamic)) variance.model$dynamic = FALSE
	if(is.null(variance.model$model)) variance.model$model = "sGARCH"
	modelinc = rep(0, 60)
	# set the custom variance target to 0
	modelinc[45] = 0
	# we allow a maximum of 4 states (in future expansion...now only 2)
	names(modelinc) = c(
			"s1.phi0", "s1.phi", "s1.xi", "s1.psi",
			"s2.phi0", "s2.phi", "s2.xi", "s2.psi", 
			"s3.phi0", "s3.phi", "s3.xi", "s3.psi", 
			"s4.phi0", "s4.phi", "s4.xi", "s4.psi",
			"psi",
			"s1.gamma","s1.c", "s1.alpha", "s1.beta",
			"s2.gamma","s2.c", "s2.alpha", "s2.beta",
			"s3.gamma","s3.c", "s3.alpha", "s3.beta",
			"sigma", "omega", "alpha", "beta", "gamma", "eta1", "eta2", 
			"delta", "lambda", "vxreg", "xi", 
			"skew", "shape", "ghlambda", rep("aux", 17))
	if(variance.model$dynamic && variance.model$model=="mixture"){
		#if(mean.model$states!=2) stop("\nstarspec-->error: mixture model only valid for 2 state case.\n")
		names(modelinc) = c(
				"s1.phi0", "s1.phi", "s1.xi", "s1.psi",
				"s2.phi0", "s2.phi", "s2.xi", "s2.psi", 
				"s3.phi0", "s3.phi", "s3.xi", "s3.psi", 
				"s4.phi0", "s4.phi", "s4.xi", "s4.psi",
				"psi",
				"s1.gamma","s1.c", "s1.alpha", "s1.beta",
				"s2.gamma","s2.c", "s2.alpha", "s2.beta",
				"s3.gamma","s3.c", "s3.alpha", "s3.beta",
				"s1.sigma", "s2.sigma", "s3.sigma", "s4.sigma", "gamma", "eta1", "eta2", 
				"delta", "lambda", "vxreg", "xi", 
				"skew", "shape", "ghlambda", rep("aux", 17))
	}
	if(intmodel==0) modelinc[51] = 0 else modelinc[51] = 1
	modeldesc = list()
	modeldata = list()
	modelinc[46] = min(4, max(1, as.integer(mean.model$states[1])))
	if(as.integer(mean.model$states[1])>4)  warning("\nstarspec-->warning: states must be an integer between 1 and 4 (setting to 4).\n")
	if(as.integer(mean.model$states[1])==0) warning("\nstarspec-->warning: states must be an integer between 1 and 4 (setting to 1).\n")
	# check the option parameters specified and warn on error
	#--------------------------------------------------------
	mm = match(names(mean.model), c("states", "include.intercept", "arOrder", "maOrder", "matype", "statevar", "s", "ylags", 
					"xreg", "statear", "yfun", "transform"))
	if(any(is.na(mm))){
		idx = which(is.na(mm))
		enx = NULL
		for(i in 1:length(idx)) enx = c(enx, names(mean.model)[idx[i]])
		warning(paste(c("\nstarspec-->warning: unidentified option(s) in mean.model:\n", enx), sep="", collapse=" "), call. = FALSE, domain = NULL)
	}
	vm = match(names(variance.model), c("dynamic", "model", "garchOrder", "submodel", "vreg", "variance.targeting"))
	if(any(is.na(vm))){
		idx = which(is.na(vm))
		enx = NULL
		for(i in 1:length(idx)) enx = c(enx, names(variance.model)[idx[i]])
		warning(paste(c("\nstarspec-->warning: unidentified option(s) in variance.model:\n", enx), sep="", collapse=" "), call. = FALSE, domain = NULL)
	}
	# distribution model
	#---------------------------------
	if(is.null(distribution.model)) modeldesc$distribution = "norm"
	valid.distribution = c("norm", "snorm", "std", "sstd","ged", "sged", "nig", "ghyp", "jsu", "ghst")
	dmodel = list()
	distribution = distribution.model[1]
	if(!is.character(distribution)) stop("\nstarspec-->error: cond.distribution argument must be a character")
	if(!any(distribution == valid.distribution)) stop("\nstarspec-->error: the cond.distribution does not appear to be a valid choice.")
	modeldesc$distribution = distribution
	modeldesc$distno = which(distribution == valid.distribution)
	di = .DistributionBounds(distribution)
	modelinc[41] = di$include.skew
	modelinc[42] = di$include.shape
	modelinc[43] = di$include.ghlambda
	modelinc[44] = modeldesc$distno
	# variance model:
	#---------------------------------
	
	if(is.null(variance.model$garchOrder)) variance.model$garchOrder = c(1,1)
	if(is.null(variance.model$variance.targeting)) variance.model$variance.targeting = FALSE
	exn=NULL
	
	if(variance.model$dynamic){
		valid.model = c("sGARCH", "eGARCH", "gjrGARCH","mixture")
		if(is.null(variance.model$model)){
			modeldesc$vmodel = "sGARCH"
		} else{
			modeldesc$vmodel = variance.model$model[1]
			if(!is.character(modeldesc$vmodel))
				stop("\nstarspec-->error: garch model argument must be a character.\n", call. = FALSE)
			if(!any(modeldesc$vmodel == valid.model)) 
				stop("\nstarspec-->error: the garch model does not appear to be a valid choice.\n", call. = FALSE)
		}
		modelinc[50] = which(modeldesc$vmodel==valid.model)
		if(modeldesc$vmodel!="mixture"){
			modeldesc$vsubmodel = NULL
			# depending on model include the additional parameters
			if(is.null(variance.model$garchOrder)){
				modelinc[32] = 1
				modelinc[33] = 1
			} else{
				modelinc[32] = variance.model$garchOrder[1]
				modelinc[33] = variance.model$garchOrder[2]
			}
			
			if( modeldesc$vmodel == "gjrGARCH" ) modelinc[34] = modelinc[32]
			if( modeldesc$vmodel == "eGARCH" ) modelinc[34] = modelinc[32]
			
			if(!is.null(variance.model$vreg)){
				if(!is.xts(variance.model$vreg)) stop("\nstarspec-->error: vreg must be an xts object.\n", call. = FALSE)
				modeldata$vexdata = variance.model$vreg
				modelinc[39] = dim( variance.model$vreg )[2]
				exn = nrow(modeldata$vexdata)
			} else{
				exn = NULL
			}
			if(is.null(variance.model$variance.targeting)){
				modelinc[31] = 1
				modeldesc$variance.targeting = FALSE
			} else{
				if(is.logical(variance.model$variance.targeting)){
					modelinc[31] = as.integer( 1 - variance.model$variance.targeting )
					modeldesc$variance.targeting = variance.model$variance.targeting
				} else{
					modelinc[31] = 0
					modelinc[45] = as.numeric(variance.model$variance.targeting)
					modeldesc$variance.targeting = TRUE
				}
			}
		} else{
			# default 2-state mixture
			if(modelinc[46]<2) stop("\nstarspec-->error: mixture model requires at least 2 states!")
			modelinc[30:31] = 1
			if(modelinc[46]==3) modelinc[32]=1
			if(modelinc[46]==4) modelinc[33]=1
		}
	} else{
		modelinc[30]=1
	}

	# mean model:
	if(!is.null(mean.model$xreg))
	{
		if(!is.xts(mean.model$xreg)) stop("\nstarspec-->error: xreg must be an xts object.\n", call. = FALSE)
		xreg = mean.model$xreg
		# check for matching indices
		if(modelinc[39]>0){
			if(!all.equal(index(xreg), index(modeldata$vexdata))){
				stop("\nstarspec-->error: vreg and xreg indices do not match!")
			}
		}
		
		modeldata$mexdata = xreg
		# if not vexdata then exn is NULL and overwritten here
		# if vexdata, then exn is the same since indices match
		exn = nrow(modeldata$mexdata)
		ncx = ncol(xreg)
		modelinc[3] = ncx
		if(modelinc[46]>1) modelinc[7] = ncx
		if(modelinc[46]>2) modelinc[11] = ncx
		if(modelinc[46]>3) modelinc[15] = ncx
	} else{
		modeldata$mexdata = NULL
	}
	
	if(is.null(mean.model$arOrder)){
		modeldata$arOrder = rep(1, modelinc[46]) 
	} else{
		if(length(mean.model$arOrder)==1) mean.model$arOrder = rep(mean.model$arOrder, modelinc[46])
		modeldata$arOrder = mean.model$arOrder[1:modelinc[46]]
		if(any(is.na(modeldata$arOrder))) stop("\nstarspec-->error: arOrder length must be either NULL, 1 or equal to states.\n")
	}
	modelinc[2] = modeldata$arOrder[1]
	if(modelinc[46]>1) modelinc[6] =  modeldata$arOrder[2]
	if(modelinc[46]>2) modelinc[10] = modeldata$arOrder[3]
	if(modelinc[46]>3) modelinc[14] = modeldata$arOrder[4]
	
	
	if(is.null(mean.model$maOrder)){
		modeldata$maOrder = rep(0, modelinc[46])
	} else{
		if(is.null(mean.model$matype)){
			modelinc[17] = mean.model$maOrder[1]
			modeldata$maOrder = mean.model$maOrder[1:modelinc[46]]
		} else{
			matype = match(tolower(mean.model$matype), c("linear","state"))
			if(is.na(matype)) stop("\nstarspec-->error: matype must be either 'linear' or 'state'.\n")
			if(matype==1){
				modelinc[17] = mean.model$maOrder[1]
			} else{
				if(length(mean.model$maOrder)==1) mean.model$maOrder = rep(mean.model$maOrder, modelinc[43])
				if(length(mean.model$maOrder)!=modelinc[46]) stop("\nstarspec-->error: maOrder length must be either NULL, 1 or equal to states.\n")
				modeldata$maOrder = mean.model$maOrder[1:modelinc[46]]
				modelinc[4] = mean.model$maOrder[1]
				if(modelinc[46]>1) modelinc[8] = mean.model$maOrder[2]
				if(modelinc[46]>2) modelinc[12] = mean.model$maOrder[3]
				if(modelinc[46]>3) modelinc[16] = mean.model$maOrder[4]
			}
		}
	}
	if(is.null(mean.model$include.intercept)) mean.model$include.intercept = rep(1,modelinc[46]) else mean.model$include.intercept = as.integer(mean.model$include.intercept)
	modelinc[1] = as.integer(mean.model$include.intercept[1])
	if(modelinc[46]>1) modelinc[5] = as.integer(mean.model$include.intercept[2])
	if(modelinc[46]>2) modelinc[9] = as.integer(mean.model$include.intercept[3])
	if(modelinc[46]>3) modelinc[13] = as.integer(mean.model$include.intercept[4])
	
	if(modelinc[46]>1) if(is.na(match(mean.model$statevar[1], c("y","s")))) stop("\nstarspec-->error: unrecognized choice for statevar\n")
	if(is.null(mean.model$statevar)){
		modelinc[49]=1
	} else{
		modelinc[49]=which(mean.model$statevar[1]==c("y","s"))
	}
	if(modelinc[49]==2){
		if(!is.xts(mean.model$s)) stop("\nstarspec-->error: s must be an xts object (aligned to index of data)\n")
		modeldata$s = mean.model$s
		if(modelinc[39]>0){
			if(!all.equal(index(modeldata$s), index(modeldata$vexdata))){
				stop("\nstarspec-->error: vreg and s indices do not match!")
			}
		}
		if(modelinc[3]>0){
			if(!all.equal(index(modeldata$s), index(modeldata$mexdata))){
				stop("\nstarspec-->error: xreg and s indices do not match!")
			}
		}
		xm = ncol(modeldata$s)
		exn = nrow(modeldata$s)
		modelinc[20] = xm
	} else{
		xm = 1
		modeldata$ylags = as.integer(mean.model$ylags)
		modelinc[20] = length(modeldata$ylags)
	}
	# check fun
	if(!is.null(mean.model$yfun)){
		modelinc[48] = 1
		if(!is.function(mean.model$yfun)) stop("\nstarspec-->error: yfun must be a function\n")
		chk = check_fun(mean.model$yfun, n=xm)
		modeldata$fun = mean.model$yfun
	} else{
		modeldata$fun = NULL
	}
	if(is.null(mean.model$statear)){
		modeldesc$statear = FALSE 
	} else{
		modeldesc$statear = as.logical(mean.model$statear)
	}
	modelinc[18:19] = 1
	if(modeldesc$statear){
		modelinc[21] = 1
	} else{
		modelinc[21] = 0
	}
	if(modelinc[46]>2){
		modelinc[22:23] = 1
		modelinc[24] = modelinc[20]
		if(modeldesc$statear) modelinc[25] = 1
	}
	if(modelinc[46]>3){
		modelinc[26:27] = 1
		modelinc[28] = modelinc[20]
		if(modeldesc$statear) modelinc[29] = 1
	}
	
	if(!is.null(fixed.prob)){
		modelinc[18:29] = 0
		modeldesc$statear=FALSE
	}
	# now check whether the supplied external data has enough rows to satisfy
	# max(AR, MA, ylags)...this is important for methods which dispatch on the
	# the STARspec class (filter and path simulation).
	mar = max(c(modelinc[c(2,6,10,14,4,8,12,16,17)], ifelse(modelinc[49]==1, max(modeldata$ylags), 0)))
	if(!is.null(exn)){
		if(exn<mar) stop(paste("\nstarspec-->error: model has ", mar, " lags but externally provided dataset has only ", exn, " rows", sep=""))
	}
	pos = 1
	pos.matrix = matrix(0, ncol = 3, nrow = 60)
	colnames(pos.matrix) = c("start", "stop", "include")
	rownames(pos.matrix) = c(
			"s1.phi0", "s1.phi", "s1.xi", "s1.psi",
			"s2.phi0", "s2.phi", "s2.xi", "s2.psi", 
			"s3.phi0", "s3.phi", "s3.xi", "s3.psi", 
			"s4.phi0", "s4.phi", "s4.xi", "s4.psi",
			"psi",
			"s1.gamma","s1.c", "s1.alpha", "s1.beta",
			"s2.gamma","s2.c", "s2.alpha", "s2.beta",
			"s3.gamma","s3.c", "s3.alpha", "s3.beta",
			"sigma", "omega", "alpha", "beta", "gamma", "eta1", "eta2", 
			"delta", "lambda", "vxreg", "xi", 
			"skew", "shape", "ghlambda", rep("aux", 17))
	if(modelinc[50]==4){
		rownames(pos.matrix) = c(
				"s1.phi0", "s1.phi", "s1.xi", "s1.psi",
				"s2.phi0", "s2.phi", "s2.xi", "s2.psi", 
				"s3.phi0", "s3.phi", "s3.xi", "s3.psi", 
				"s4.phi0", "s4.phi", "s4.xi", "s4.psi",
				"psi",
				"s1.gamma","s1.c", "s1.alpha", "s1.beta",
				"s2.gamma","s2.c", "s2.alpha", "s2.beta",
				"s3.gamma","s3.c", "s3.alpha", "s3.beta",
				"s1.sigma", "s2.sigma", "s3.sigma", "s4.sigma", "gamma", "eta1", "eta2", 
				"delta", "lambda", "vxreg", "xi", 
				"skew", "shape", "ghlambda", rep("aux", 17))
	}
	# check if there are starting or fixed
	# check that it is included in the optimization
	# check that it is needed in the model
	for(i in 1:43){
		if( modelinc[i] > 0 ){
			pos.matrix[i,1:3] = c(pos, pos+modelinc[i]-1, 1)
			pos = max(pos.matrix[1:i,2]+1)
		}
	}
	modelnames = .makestarnames(modelinc)
	
	nn = length(modelnames)
	modelmatrix = matrix(0, ncol = 3, nrow = nn)
	rownames(modelmatrix) = modelnames
	colnames(modelmatrix) = c("opt", "fixed", "start")
	fixed.names = names(fixed.pars)
	fp = charmatch(fixed.names, modelnames)
	
	if(!is.null(fixed.names) && any(!is.na(fp))){
		fixed = fp[!is.na(fp)]
		modelmatrix[fixed,2] = 1
		fz = charmatch(modelnames, fixed.names)
		fz = fz[!is.na(fz)]
		fixed.pars = fixed.pars[fz]
		names(fixed.pars) = fixed.names[fz]
	} else{
		fixed.pars = NULL
	}
	modelmatrix[,1] = 1 - modelmatrix[,2]
	start.names = names(start.pars)
	sp = charmatch(start.names, modelnames)
	if(!is.null(start.names) && any(!is.na(sp))){
		start = sp[!is.na(sp)]
		modelmatrix[start,3] = 1
		sz = charmatch(modelnames, start.names)
		sz = sz[!is.na(sz)]
		start.pars = start.pars[sz]
	} else{
		start.pars = NULL
	}
	##################################################################
	# Parameter Matrix
	if(modelinc[50]==4){
		mm = sum(modelinc[c(2:4,6:8,10:12,14:16,17,20,24,28,34:36,39)])
		mm = mm - length( which(modelinc[c(2:4,6:8,10:12,14:16,17,20,24,28,34:36,39)]>0) )
		pars = matrix(0, ncol = 8, nrow = 43 + mm)
	} else{
		mm = sum(modelinc[c(2:4,6:8,10:12,14:16,17,20,24,28,32:36,39)])
		mm = mm - length( which(modelinc[c(2:4,6:8,10:12,14:16,17,20,24,28,32:36,39)]>0) )
		pars = matrix(0, ncol = 8, nrow = 43 + mm)
	}
	colnames(pars) = c("Level", "Fixed", "Include", "Estimate", "LB", "UB", "Transformed", "Type")
	# Type: [1=Linear, 2=State, 3 = Other]
	pidx = matrix(NA, nrow = 43, ncol = 2)
	colnames(pidx) = c("begin", "end")
	rownames(pidx) =   c(
			"s1.phi0", "s1.phi", "s1.xi", "s1.psi",
			"s2.phi0", "s2.phi", "s2.xi", "s2.psi", 
			"s3.phi0", "s3.phi", "s3.xi", "s3.psi", 
			"s4.phi0", "s4.phi", "s4.xi", "s4.psi",
			"psi",
			"s1.gamma","s1.c", "s1.alpha", "s1.beta",
			"s2.gamma","s2.c", "s2.alpha", "s2.beta",
			"s3.gamma","s3.c", "s3.alpha", "s3.beta",
			"sigma", "omega", "alpha", "beta", "gamma", "eta1", "eta2", 
			"delta", "lambda", "vxreg", "xi", 
			"skew", "shape", "ghlambda")
	if(modelinc[50]==4){
		rownames(pidx) = c(
				"s1.phi0", "s1.phi", "s1.xi", "s1.psi",
				"s2.phi0", "s2.phi", "s2.xi", "s2.psi", 
				"s3.phi0", "s3.phi", "s3.xi", "s3.psi", 
				"s4.phi0", "s4.phi", "s4.xi", "s4.psi",
				"psi",
				"s1.gamma","s1.c", "s1.alpha", "s1.beta",
				"s2.gamma","s2.c", "s2.alpha", "s2.beta",
				"s3.gamma","s3.c", "s3.alpha", "s3.beta",
				"s1.sigma", "s2.sigma", "s3.sigma", "s4.sigma", "gamma", "eta1", "eta2", 
				"delta", "lambda", "vxreg", "xi", 
				"skew", "shape", "ghlambda")
	}
	fixed.names = names(fixed.pars)
	pnames = NULL
	nx = 0

	#----------------------------------------
	if(pos.matrix[1,3]==1){
		pars[1, 3] = 1
		pars[1, 1] = 0
		if(any(substr(fixed.names, 1, 7)=="s1.phi0")) pars[1,2] = 1 else pars[1,4] = 1
		pars[1,7] = 1
		pars[1,8] = 1
	}
	pidx[1,1] = 1
	pidx[1,2] = 1
	pnames = c(pnames, "s1.phi0")
	nx = 1
	pn = 1
	#----------------------------------------
	pidx[2,1] = 2
	if(pos.matrix[2,3] == 1){
		pn = length( seq(pos.matrix[2,1], pos.matrix[2,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			pars[(nx+i), 8] = 1
			nnx = paste("s1.phi", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
		# only transform the fist ar parameter to be bounded
		pars[(nx+1),7] = 1
	} else{
		pnames = c(pnames, "s1.phi")
	}
	pidx[2,2] = 1+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[3,1] = nx+1
	if(pos.matrix[3,3] == 1){
		pn = length( seq(pos.matrix[3,1], pos.matrix[3,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			pars[(nx+i), 8] = 1
			nnx = paste("s1.xi", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "s1.xi")
	}
	pidx[3,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[4,1] = nx+1
	if(pos.matrix[4,3] == 1){
		pn = length( seq(pos.matrix[4,1], pos.matrix[4,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			pars[(nx+i), 8] = 1
			nnx = paste("s1.psi", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
		pars[(nx+1),7] = 1
	} else{
		pnames = c(pnames, "s1.psi")
	}
	pidx[4,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[5,1] = nx+1
	if(pos.matrix[5,3] == 1){
		pars[nx+pn, 1] = 0
		pars[nx+pn, 3] = 1
		pars[nx+pn, 8] = 1
		if(any(substr(fixed.names, 1, 7)=="s2.phi0")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
		pars[nx+pn,7] = 1
	}
	pidx[5,2] = nx+pn
	pnames = c(pnames, "s2.phi0")
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[6,1] = nx+1
	if(pos.matrix[6,3] == 1){
		pn = length( seq(pos.matrix[6,1], pos.matrix[6,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			pars[(nx+i), 8] = 1
			nnx = paste("s2.phi", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
		pars[nx+1,7] = 1
	} else{
		pnames = c(pnames, "s2.phi")
	}
	pidx[6,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[7,1] = nx+1
	if(pos.matrix[7,3] == 1){
		pn = length( seq(pos.matrix[7,1], pos.matrix[7,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			pars[(nx+i), 8] = 1
			nnx = paste("s2.xi", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "s2.xi")
	}
	pidx[7,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[8,1] = nx+1
	if(pos.matrix[8,3] == 1){
		pn = length( seq(pos.matrix[8,1], pos.matrix[8,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			pars[(nx+i), 8] = 1
			nnx = paste("s2.psi", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
		pars[(nx+1),7] = 1
	} else{
		pnames = c(pnames, "s2.psi")
	}
	pidx[8,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[9,1] = nx+1
	if(pos.matrix[9,3] == 1){
		pars[nx+pn, 1] = 0
		pars[nx+pn, 3] = 1
		pars[nx+pn, 8] = 1
		if(any(substr(fixed.names, 1, 7)=="s3.phi0")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
		pars[nx+pn,7] = 1
	}
	pidx[9,2] = nx+pn
	pnames = c(pnames, "s3.phi0")
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[10,1] = nx+1
	if(pos.matrix[10,3] == 1){
		pn = length( seq(pos.matrix[10,1], pos.matrix[10,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			pars[(nx+i), 8] = 1
			nnx = paste("s3.phi", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
		pars[nx+1,7] = 1
	} else{
		pnames = c(pnames, "s3.phi")
	}
	pidx[10,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[11,1] = nx+1
	if(pos.matrix[11,3] == 1){
		pn = length( seq(pos.matrix[11,1], pos.matrix[11,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			pars[(nx+i), 8] = 1
			nnx = paste("s3.xi", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "s3.xi")
	}
	pidx[11,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[12,1] = nx+1
	if(pos.matrix[12,3] == 1){
		pn = length( seq(pos.matrix[12,1], pos.matrix[12,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			pars[(nx+i), 8] = 1
			nnx = paste("s3.psi", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
		pars[(nx+1),7] = 1
	} else{
		pnames = c(pnames, "s3.psi")
	}
	pidx[12,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[13,1] = nx+1
	if(pos.matrix[13,3] == 1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		pars[nx+pn, 8] = 1
		if(any(substr(fixed.names, 1, 7)=="s4.phi0")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
		pars[nx+pn,7] = 1
	}
	pidx[13,2] = nx+pn
	pnames = c(pnames, "s4.phi0")
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[14,1] = nx+1
	if(pos.matrix[14,3] == 1){
		pn = length( seq(pos.matrix[14,1], pos.matrix[14,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			pars[(nx+i), 8] = 1
			nnx = paste("s4.phi", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
		pars[nx+1,7] = 1
	} else{
		pnames = c(pnames, "s4.phi")
	}
	pidx[14,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[15,1] = nx+1
	if(pos.matrix[15,3] == 1){
		pn = length( seq(pos.matrix[15,1], pos.matrix[15,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			pars[(nx+i), 8] = 1
			nnx = paste("s4.xi", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "s4.xi")
	}
	pidx[15,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[16,1] = nx+1
	if(pos.matrix[16,3] == 1){
		pn = length( seq(pos.matrix[16,1], pos.matrix[16,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			pars[(nx+i), 8] = 1
			nnx = paste("s4.psi", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
		pars[(nx+1),7] = 1
	} else{
		pnames = c(pnames, "s4.psi")
	}
	pidx[16,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[17,1] = nx+1
	if(pos.matrix[17,3] == 1){
		pn = length( seq(pos.matrix[17,1], pos.matrix[17,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			pars[(nx+i), 8] = 1
			nnx = paste("psi", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "psi")
	}
	pidx[17,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[18,1] = nx+1
	if(pos.matrix[18,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		pars[nx+pn, 8] = 2
		if(any(substr(fixed.names, 1, 8)=="s1.gamma")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pnames = c(pnames, "s1.gamma")
	pidx[18,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[19,1] = nx+1
	if(pos.matrix[19,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		pars[nx+pn, 8] = 2
		if(any(substr(fixed.names, 1, 4)=="s1.c")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pnames = c(pnames, "s1.c")
	pidx[19,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[20,1] = nx+1
	if(pos.matrix[20,3] == 1){
		pn = length( seq(pos.matrix[20,1], pos.matrix[20,2], by = 1) )
		for(i in 1:pn){
			if(i==1){
				pars[(nx+i), 1] = 1
				pars[(nx+i), 3] = 1
				pars[(nx+i), 2] = 1
				pars[(nx+i), 4] = 0
				pars[(nx+i), 8] = 2
				nnx = paste("s1.alpha", i, sep="")
				pnames = c(pnames, nnx)
			} else{
				pars[(nx+i), 1] = 0
				pars[(nx+i), 3] = 1
				pars[(nx+i), 8] = 2
				nnx = paste("s1.alpha", i, sep="")
				sp = na.omit(match(fixed.names, nnx))
				if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
				pnames = c(pnames, nnx)
			}
		}
	} else{
		pnames = c(pnames, "s1.alpha")
	}
	pidx[20,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[21,1] = nx+1
	if(pos.matrix[21,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		pars[nx+pn, 8] = 2
		if(any(substr(fixed.names, 1, 7)=="s1.beta")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
		pars[nx+pn,7] = 1
	}
	pnames = c(pnames, "s1.beta")
	pidx[21,2] = nx+pn
	nx = nx + pn
	pn = 1
	
	#----------------------------------------
	pidx[22,1] = nx+1
	if(pos.matrix[22,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		pars[nx+pn, 8] = 2
		if(any(substr(fixed.names, 1, 8)=="s2.gamma")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pnames = c(pnames, "s2.gamma")
	pidx[22,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[23,1] = nx+1
	if(pos.matrix[23,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		pars[nx+pn, 8] = 2
		if(any(substr(fixed.names, 1, 4)=="s2.c")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pnames = c(pnames, "s2.c")
	pidx[23,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[24,1] = nx+1
	if(pos.matrix[24,3] == 1){
		pn = length( seq(pos.matrix[24,1], pos.matrix[24,2], by = 1) )
		for(i in 1:pn){
			if(i==1){
				pars[(nx+i), 1] = 1
				pars[(nx+i), 3] = 1
				pars[(nx+i), 2] = 1
				pars[(nx+i), 4] = 0
				pars[(nx+i), 8] = 2
				nnx = paste("s2.alpha", i, sep="")
				pnames = c(pnames, nnx)
			} else{
				pars[(nx+i), 1] = 0
				pars[(nx+i), 3] = 1
				pars[(nx+i), 8] = 2
				nnx = paste("s2.alpha", i, sep="")
				sp = na.omit(match(fixed.names, nnx))
				if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
				pnames = c(pnames, nnx)
			}
		}
	} else{
		pnames = c(pnames, "s2.alpha")
	}
	pidx[24,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[25,1] = nx+1
	if(pos.matrix[25,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		pars[nx+pn, 8] = 2
		if(any(substr(fixed.names, 1, 7)=="s2.beta")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
		pars[nx+pn,7] = 1
	}
	pnames = c(pnames, "s2.beta")
	pidx[25,2] = nx+pn
	nx = nx + pn
	pn = 1
	
	
	#----------------------------------------
	pidx[26,1] = nx+1
	if(pos.matrix[26,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		pars[nx+pn, 8] = 2
		if(any(substr(fixed.names, 1, 8)=="s3.gamma")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pnames = c(pnames, "s3.gamma")
	pidx[26,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[27,1] = nx+1
	if(pos.matrix[27,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		pars[nx+pn, 8] = 2
		if(any(substr(fixed.names, 1, 4)=="s3.c")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pnames = c(pnames, "s3.c")
	pidx[27,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[28,1] = nx+1
	if(pos.matrix[28,3] == 1){
		pn = length( seq(pos.matrix[28,1], pos.matrix[28,2], by = 1) )
		for(i in 1:pn){
			if(i==1){
				pars[(nx+i), 1] = 1
				pars[(nx+i), 3] = 1
				pars[(nx+i), 2] = 1
				pars[(nx+i), 4] = 0
				pars[(nx+i), 8] = 2
				nnx = paste("s3.alpha", i, sep="")
				pnames = c(pnames, nnx)
			} else{
				pars[(nx+i), 1] = 0
				pars[(nx+i), 3] = 1
				pars[(nx+i), 8] = 2
				nnx = paste("s3.alpha", i, sep="")
				sp = na.omit(match(fixed.names, nnx))
				if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
				pnames = c(pnames, nnx)
			}
		}
	} else{
		pnames = c(pnames, "s3.alpha")
	}
	pidx[28,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[29,1] = nx+1
	if(pos.matrix[29,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		pars[nx+pn, 8] = 2
		if(any(substr(fixed.names, 1, 7)=="s3.beta")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
		pars[nx+pn,7] = 1
	}
	pnames = c(pnames, "s3.beta")
	pidx[29,2] = nx+pn
	nx = nx + pn
	pn = 1
	
	#----------------------------------------
	pidx[30,1] = nx+1
	if(pos.matrix[30,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		pars[nx+pn, 8] = 3
		if(modelinc[50]==4){
			if(any(substr(fixed.names, 1, 8)=="s1.sigma")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
			pars[nx+pn,7] = 1
		} else{
			if(any(substr(fixed.names, 1, 5)=="sigma")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
			pars[nx+pn,7] = 1
		}
	}
	if(modelinc[50]==4) pnames = c(pnames, "s1.sigma") else pnames = c(pnames, "sigma")
	pidx[30,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	
	pidx[31,1] = nx+1
	if(pos.matrix[31,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		pars[nx+pn, 8] = 3
		if(modelinc[50]==4){
			if(any(substr(fixed.names, 1, 8)=="s2.sigma")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
			pars[nx+pn,7] = 1
		} else{
			if(any(substr(fixed.names, 1, 5)=="omega")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
			pars[nx+pn,7] = 1
		}		
	}
	if(modelinc[50]==4) pnames = c(pnames, "s2.sigma") else pnames = c(pnames, "omega")
	pidx[31,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------	
	pidx[32,1] = nx+1
	if(pos.matrix[32,3]==1){
		pn = length( seq(pos.matrix[32,1], pos.matrix[32,2], by = 1) )
		if(modelinc[50]==4){
			pars[nx+pn, 3] = 1
			pars[nx+pn, 1] = 0
			pars[nx+pn, 8] = 3
			if(any(substr(fixed.names, 1, 8)=="s3.sigma")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
			pars[nx+pn,7] = 1
			pnames = c(pnames, "s3.sigma")			
		} else{
			for(i in 1:pn){
				pars[(nx+i), 1] = 0
				pars[(nx+i), 3] = 1
				pars[(nx+i), 8] = 3
				nnx = paste("alpha", i, sep="")
				sp = na.omit(match(fixed.names, nnx))
				if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
				pnames = c(pnames, nnx)
			}
			pars[nx+1,7] = 1
		}
	} else{
		if(modelinc[50]==4) pnames = c(pnames, "s3.sigma") else pnames = c(pnames, "alpha")
	}
	pidx[32,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[33,1] = nx+1
	if(pos.matrix[33,3]==1){
		pn = length( seq(pos.matrix[33,1], pos.matrix[33,2], by = 1) )
		if(modelinc[50]==4){
			pars[nx+pn, 3] = 1
			pars[nx+pn, 1] = 0
			pars[nx+pn, 8] = 3
			if(any(substr(fixed.names, 1, 8)=="s4.sigma")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
			pars[nx+pn,7] = 1
			pnames = c(pnames, "s4.sigma")
		} else{
			for(i in 1:pn){
				pars[(nx+i), 1] = 0
				pars[(nx+i), 3] = 1
				pars[(nx+i), 8] = 3
				nnx = paste("beta", i, sep="")
				sp = na.omit(match(fixed.names, nnx))
				if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
				pnames = c(pnames, nnx)
			}
			pars[nx+1,7] = 1
			#-------------------------------------------
			# special consideration for the iGARCH model
			#-------------------------------------------
			if(modeldesc$vmodel == "iGARCH"){
				# last beta not estimated
				pars[nx+pn, 4] = 0
				nnx = paste("beta", pn, sep="")
				# do not allow the last beta to be fixed
				if(any(substr(fixed.names, 1, nchar(nnx))==nnx)) pars[(nx+pn), 2] = 0
			}
		}
	} else{
		if(modelinc[50]==4) pnames = c(pnames, "s4.sigma") else pnames = c(pnames, "beta")
	}
	pidx[33,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[34,1] = nx+1
	if(pos.matrix[34,3]==1){
		pn = length( seq(pos.matrix[34,1], pos.matrix[34,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			pars[(nx+i), 8] = 3
			nnx = paste("gamma", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
		pars[nx+1,7] = 1
	} else{
		pnames = c(pnames, "gamma")
	}
	pidx[34,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[35,1] = nx+1
	if(pos.matrix[35,3]==1){
		pn = length( seq(pos.matrix[35,1], pos.matrix[35,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			pars[(nx+i), 8] = 3
			nnx = paste("eta1", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
		pars[nx+1,7] = 1
	} else{
		pnames = c(pnames, "eta1")
	}
	pidx[35,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[36,1] = nx+1
	if(pos.matrix[36,3]==1){
		pn = length( seq(pos.matrix[36,1], pos.matrix[36,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			pars[(nx+i), 8] = 3
			nnx = paste("eta2", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
		pars[nx+1,7] = 1
	} else{
		pnames = c(pnames, "eta2")
	}
	pidx[36,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[37,1] = nx+1
	
	if(pos.matrix[37,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		pars[nx+pn, 8] = 3
		if(any(substr(fixed.names, 1, 5)=="delta")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
		pars[nx+pn,7] = 1
	} else{
		
	}
	pidx[37,2] = nx+pn
	pnames = c(pnames, "delta")
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[38,1] = nx+1
	
	if(pos.matrix[38,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		pars[nx+pn, 8] = 3		
		if(any(substr(fixed.names, 1, 6)=="lambda")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
		pars[nx+pn,7] = 1
	} else{
	}
	pidx[38,2] = nx+pn
	pnames = c(pnames, "lambda")
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[39,1] = nx+1
	if(pos.matrix[39,3]==1){
		pn = length( seq(pos.matrix[39,1], pos.matrix[39,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			pars[(nx+i), 8] = 3
			nnx = paste("vxreg", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "vxreg")
		
	}
	pidx[39,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[40,1] = nx+1
	
	if(pos.matrix[40,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		pars[nx+pn, 8] = 3
		if(any(substr(fixed.names, 1, 2)=="xi")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
		pars[nx+pn,7] = 1
	}
	pidx[40,2] = nx+pn
	pnames = c(pnames, "xi")
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[41,1] = nx+1
	
	if(pos.matrix[41,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		pars[nx+pn, 8] = 3
		if(any(substr(fixed.names, 1, 4)=="skew")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
		pars[nx+pn,7] = 1
	}
	pidx[41,2] = nx+pn
	pnames = c(pnames, "skew")
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[42,1] = nx+1
	if(pos.matrix[42,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		pars[nx+pn, 8] = 3
		if(any(substr(fixed.names, 1, 5)=="shape")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
		pars[nx+pn,7] = 1
	}
	pnames = c(pnames, "shape")
	pidx[42,2] = nx+pn
	nx = nx + pn
	pn = 1
	#----------------------------------------
	pidx[43,1] = nx+1
	
	if(pos.matrix[43,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		pars[nx+pn, 8] = 3
		if(any(substr(fixed.names, 1, 8)=="ghlambda")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
		pars[nx+pn,7] = 1
	}
	pidx[43,2] = nx+pn
	#----------------------------------------
	
	# Once more for fgarch model pars
	pnames = c(pnames, "ghlambda")
	rownames(pars) = pnames
	
	zf = match(fixed.names, rownames(pars))
	if( length(zf)>0 ) pars[zf, 1] = unlist(fixed.pars)
	pars[,"LB"] = NA
	pars[,"UB"] = NA
	
	if(!is.null(fixed.prob))
	{
		modelinc[47] = 1
		if(!is.xts(fixed.prob)) stop("\nstarspec-->error: fixed.prob must be an xts object")
		if(modelinc[46]==1){
			# In the case that states=1 and the use passes a set of probabilities, these
			# are in fact weights (weighted ARX model), else default is to set them to 1.
			fixed.prob = fixed.prob[,1]
			colnames(fixed.prob) = "Weights"
		}
		if(modelinc[46]==2){
			fixed.prob = cbind(fixed.prob[,1], 1-fixed.prob[,1])
			colnames(fixed.prob) = c("Prob[S=1]","Prob[S=2]")
		}
		if(modelinc[46]==3){
			if(ncol(fixed.prob)==2){
				fixed.prob = cbind(fixed.prob[,1], fixed.prob[,2], 1-fixed.prob[,1]-fixed.prob[,2])
				colnames(fixed.prob) = c("Prob[S=1]","Prob[S=2]", "Prob[S=3]")
				test = rowSums(fixed.prob)
				if(any(test!=1.0)) warning("\nstarspec-->warning: probabilities do not all add upto 1 (suggest to check and resubmit).\n")			
			} else if(ncol(fixed.prob)==3){
				colnames(fixed.prob) = c("Prob[S=1]","Prob[S=2]", "Prob[S=3]")
				test = rowSums(fixed.prob)
				if(any(test!=1.0)) warning("\nstarspec-->warning: probabilities do not all add upto 1 (suggest to check and resubmit).\n")			
			} else{
				stop("\nstarspec-->error: fixed.prob must have either 2 or 3 columns for the 3 state case.\n")
			}
		}
		if(modelinc[46]==4){
			if(ncol(fixed.prob)==3){
				fixed.prob = cbind(fixed.prob[,1], fixed.prob[,2], fixed.prob[,3], 1-fixed.prob[,1]-fixed.prob[,2]-fixed.prob[,3])
				colnames(fixed.prob) = c("Prob[S=1]","Prob[S=2]", "Prob[S=3]", "Prob[S=4]")
				test = rowSums(fixed.prob)
				if(any(test!=1.0)) warning("\nstarspec-->warning: probabilities do not all add upto 1 (suggest to check and resubmit).\n")			
			} else if(ncol(fixed.prob)==4){
				colnames(fixed.prob) = c("Prob[S=1]","Prob[S=2]","Prob[S=3]","Prob[S=4]")
				test = rowSums(fixed.prob)
				if(any(test!=1.0)) warning("\nstarspec-->warning: probabilities do not all add upto 1 (suggest to check and resubmit).\n")			
			} else{
				stop("\nstarspec-->error: fixed.prob must have either 3 or 4 columns for the 4 state case.\n")
			}
		}
	} else{
		if(modelinc[46]==1){
			stop("\nstarspec-->error: fixed.prob must be provided for the 1-state case.\n")
		}
	}
	
	model = list(modelinc = modelinc, modeldesc = modeldesc, modeldata = modeldata, pars = pars, 
			start.pars = start.pars, fixed.pars = fixed.pars,  pos.matrix = pos.matrix, 
			pidx = pidx, fixed.prob = fixed.prob)
	ans = new("STARspec", model = model)
	return(ans)

}

setMethod(f = "starspec", definition = .xstarspec)


# CONTINUE HERE
#------------------------------------------------------------------------------------
# get spec from fitted object
#------------------------------------------------------------------------------------
.getstarspec = function(object)
{
	model = object@model
	modelinc = model$modelinc
	if(modelinc[17]>0) maOrder=modelinc[17] else maOrder=c(modelinc[4],modelinc[8],modelinc[12], modelinc[16])[1:modelinc[46]]
	vt = model$modeldesc$variance.targeting
	
	spec = starspec(mean.model = list(states = modelinc[46],
					include.intercept = c(modelinc[1],modelinc[5], modelinc[9], modelinc[13]), 
					arOrder = c(modelinc[2],modelinc[6],modelinc[10], modelinc[14])[1:modelinc[46]],
					maOrder = maOrder,
					matype = ifelse(modelinc[17]>0, "linear","state"),
					statevar = c("y","s")[modelinc[49]], 
					statear = as.logical(model$modeldesc$statear), 
    				ylags = model$modeldata$ylags, yfun = model$modeldata$fun,
    				s = model$modeldata$s, xreg = model$modeldata$mexdata), 
			variance.model = list(dynamic = as.logical(modelinc[50]>0), 
    				model = model$modeldesc$vmodel, garchOrder = modelinc[32:33], 
					submodel = model$modeldesc$vsubmodel, 
					vreg = model$modeldata$vexdata,
					variance.targeting = vt), 
    		distribution.model = model$modeldesc$distribution, 
			start.pars = model$start.pars, fixed.pars = model$fixed.pars,
			fixed.prob = model$fixed.prob)
	idx = which(is.na(spec@model$pars[,"LB"]))
	spec@model$pars[idx,"LB"] = object@model$pars[idx,"LB"]
	idx = which(is.na(spec@model$pars[,"UB"]))
	spec@model$pars[idx,"UB"] = object@model$pars[idx,"UB"]
	return(spec)
}

setMethod(f = "getspec", signature(object = "STARfit"), definition = .getstarspec)

#------------------------------------------------------------------------------------
# fixed parameters
#------------------------------------------------------------------------------------
.setstarfixed = function(object, value){
	# get parameter values
	model = object@model
	modelinc = model$modelinc
	ipars = model$pars
	pars = unlist(value)
	names(pars) = parnames = tolower(names(pars))
	# included parameters in model
	modelnames = rownames(ipars[which(ipars[,4]==1 | ipars[,2]==1), ,drop=FALSE])
	inc = NULL
	for(i in seq_along(parnames)){
		if(is.na(match(parnames[i], modelnames))){
			warning( (paste("\nUnrecognized Parameter in Fixed Values: ", parnames[i], "...Ignored", sep = "")))
		} else{
			inc = c(inc, i)
		}
	}
	fixed.pars = pars[inc]
	names(fixed.pars) = tolower(names(pars[inc]))
	# set parameter values
	if(modelinc[17]>0) maOrder=modelinc[17] else maOrder=c(modelinc[4],modelinc[8],modelinc[12], modelinc[16])[1:modelinc[46]]
	vt = model$modeldesc$variance.targeting
	
	spec = starspec(mean.model = list(states = modelinc[46],
					include.intercept = c(modelinc[1],modelinc[5], modelinc[9], modelinc[13]), 
					arOrder = c(modelinc[2],modelinc[6],modelinc[10], modelinc[14])[1:modelinc[46]],
					maOrder = maOrder,
					matype = ifelse(modelinc[17]>0, "linear","state"),
					statevar = c("y","s")[modelinc[49]], 
					statear = as.logical(model$modeldesc$statear), 
					ylags = model$modeldata$ylags, yfun = model$modeldata$fun,
					s = model$modeldata$s, xreg = model$modeldata$mexdata), 
			variance.model = list(dynamic = as.logical(modelinc[50]>0), 
					model = model$modeldesc$vmodel, garchOrder = modelinc[32:33], 
					submodel = model$modeldesc$vsubmodel, 
					vreg = model$modeldata$vexdata, variance.targeting =vt), 
			distribution.model = model$modeldesc$distribution, 
			start.pars = model$start.pars, fixed.pars = as.list(fixed.pars),
			fixed.prob = model$fixed.prob)
	idx = which(is.na(spec@model$pars[,"LB"]))
	spec@model$pars[idx,"LB"] = object@model$pars[idx,"LB"]
	spec@model$pars[idx,"Transformed"] = object@model$pars[idx,"Transformed"]
	idx = which(is.na(spec@model$pars[,"UB"]))
	spec@model$pars[idx,"UB"] = object@model$pars[idx,"UB"]
	spec@model$pars[idx,"Transformed"] = object@model$pars[idx,"Transformed"]
	return(spec)
}

setReplaceMethod(f="setfixed", signature= c(object = "STARspec", value = "vector"), definition = .setstarfixed)

#------------------------------------------------------------------------------------
# starting parameters
#------------------------------------------------------------------------------------
.setstarstart = function(object, value){
	# get parameter values
	model = object@model
	modelinc = model$modelinc
	ipars = model$pars
	pars = unlist(value)
	names(pars) = parnames = tolower(names(pars))
	# included parameters in model
	modelnames = rownames(ipars[which(ipars[,4]==1 | ipars[,2]==1), ,drop=FALSE])
	inc = NULL
	for(i in seq_along(parnames)){
		if(is.na(match(parnames[i], modelnames))){
			warning( (paste("\nUnrecognized Parameter in Start Values: ", parnames[i], "...Ignored", sep = "")))
		} else{
			inc = c(inc, i)
		}
	}
	start.pars = pars[inc]
	names(start.pars) = tolower(names(pars[inc]))
	# set parameter values
	if(modelinc[17]>0) maOrder=modelinc[17] else maOrder=c(modelinc[4],modelinc[8],modelinc[12], modelinc[16])[1:modelinc[46]]
	vt = model$modeldesc$variance.targeting
	
	spec = starspec(mean.model = list(states = modelinc[46],
					include.intercept = c(modelinc[1],modelinc[5], modelinc[9], modelinc[13]), 
					arOrder = c(modelinc[2],modelinc[6],modelinc[10], modelinc[14])[1:modelinc[46]],
					maOrder = maOrder,
					matype = ifelse(modelinc[17]>0, "linear","state"),
					statevar = c("y","s")[modelinc[49]], 
					statear = as.logical(model$modeldesc$statear), 
					ylags = model$modeldata$ylags, yfun = model$modeldata$fun,
					s = model$modeldata$s, xreg = model$modeldata$mexdata), 
			variance.model = list(dynamic = as.logical(modelinc[50]>0), 
					model = model$modeldesc$vmodel, garchOrder = modelinc[32:33], 
					submodel = model$modeldesc$vsubmodel, 
					vreg = model$modeldata$vexdata, variance.targeting =  vt), 
			distribution.model = model$modeldesc$distribution, 
			fixed.pars = model$fixed.pars, start.pars = as.list(start.pars),
			fixed.prob = model$fixed.prob)
	# ToDo : Check that the starting pars are not outside the upper and lower bounds
	idx = which(is.na(spec@model$pars[,"LB"]))
	spec@model$pars[idx,"LB"] = object@model$pars[idx,"LB"]
	spec@model$pars[idx,"Transformed"] = object@model$pars[idx,"Transformed"]
	idx = which(is.na(spec@model$pars[,"UB"]))
	spec@model$pars[idx,"UB"] = object@model$pars[idx,"UB"]
	spec@model$pars[idx,"Transformed"] = object@model$pars[idx,"Transformed"]
	return(spec)
}

setReplaceMethod(f="setstart", signature= c(object = "STARspec", value = "vector"), definition = .setstarstart)

.checkallfixed = function( spec ){
	# check that a given spec with fixed parameters
	model = spec@model
	pars = model$pars
	pnames = rownames(pars)
	estpars = pnames[as.logical(pars[,2] * pars[,3] + pars[,3] * pars[,4])]
	return( estpars )
}

#------------------------------------------------------------------------------------
# Set parameter bounds
#------------------------------------------------------------------------------------
.setstarbounds = function(object, value){
	model = object@model
	ipars = model$pars
	parnames = tolower(names(value))
	# included parameters in model
	modelnames = rownames(ipars[which(ipars[,4] == 1), ])
	sp = na.omit(match(parnames, modelnames))
	if(length(sp)>0){
		for(i in 1:length(sp)){
			#if(length(value[[modelnames[sp[i]]]])!=2)
			ipars[modelnames[sp[i]], 5] = as.numeric(value[[modelnames[sp[i]]]][1])
			ipars[modelnames[sp[i]], 6] = as.numeric(value[[modelnames[sp[i]]]][2])
			ipars[modelnames[sp[i]], 7] = 1
		}
	}
	object@model$pars = ipars
	return(object)
}
setReplaceMethod(f="setbounds", signature= c(object = "STARspec", value = "vector"), definition = .setstarbounds)

#------------------------------------------------------------------------------------
# Fitting (ML estimation)
#------------------------------------------------------------------------------------
starfit = function(spec, data, out.sample = 0, solver = "optim", solver.control = list(), 
		fit.control = list(stationarity = 0, fixed.se = 0, rec.init = 'all'), cluster = NULL, 
		n = 25, ...)
{
	UseMethod("starfit")
}

setMethod("starfit", signature(spec = "STARspec"), .starfit)

#------------------------------------------------------------------------------------
# Filter
#------------------------------------------------------------------------------------
starfilter = function(spec, data, out.sample = 0, n.old = NULL, rec.init = 'all', ...)
{
	UseMethod("starfilter")
}

setMethod("starfilter", signature(spec = "STARspec"), .starfilter)

#------------------------------------------------------------------------------------
# Forecast (1 and n-ahead)
#------------------------------------------------------------------------------------
starforecast = function(fitORspec, data = NULL, n.ahead = 1, n.roll = 0, out.sample = 0, 
		external.forecasts = list(xregfor = NULL, vregfor = NULL, sfor = NULL, probfor = NULL), 
		method = c("an.parametric", "an.kernel", "mc.empirical", "mc.parametric", "mc.kernel"), 
		mc.sims = NULL, ...)
{
	UseMethod("starforecast")
}

setMethod("starforecast", signature(fitORspec = "STARfit"), .starforecast.fit)
setMethod("starforecast", signature(fitORspec = "STARspec"), .starforecast.spec)

#------------------------------------------------------------------------------------
# Simulation
#------------------------------------------------------------------------------------
starsim = function(fit, n.sim = 1000, n.start = 0, m.sim = 1, presigma = NA, prereturns = NA, preresiduals = NA, rseed = NA, 
		custom.dist = list(name = NA, distfit = NA), xregsim = NULL, vregsim = NULL, ssim = NULL, probsim = NULL, ...)
{
	UseMethod("starsim")
}

.starsim = function(fit, n.sim = 1000, n.start = 0, m.sim = 1, presigma = NA, 
		prereturns = NA, preresiduals = NA, rseed = NA, 
		custom.dist = list(name = NA, distfit = NA), xregsim = NULL, 
		vregsim = NULL, ssim = NULL, probsim = NULL)
{
	if(fit@model$modelinc[50]>0){
		if(fit@model$modelinc[50]==4){
			ans = .starsim.mixture(fit = fit, n.sim = n.sim, n.start = n.start, m.sim = m.sim,
					prereturns = prereturns, preresiduals = preresiduals, rseed = rseed, 
					custom.dist = custom.dist, xregsim = xregsim, ssim = ssim, 
					probsim = probsim)
		} else{
			ans = .starsim.dynamic(fit = fit, n.sim = n.sim, n.start = n.start, m.sim = m.sim,
					presigma = presigma, prereturns = prereturns, preresiduals = preresiduals, 
					rseed = rseed, custom.dist = custom.dist, xregsim = xregsim, 
					vregsim = vregsim, ssim = ssim, probsim = probsim)
		}
	} else{
		ans = .starsim.static(fit = fit, n.sim = n.sim, n.start = n.start, m.sim = m.sim, 
				prereturns = prereturns, preresiduals = preresiduals, rseed = rseed, 
				custom.dist = custom.dist, xregsim = xregsim, ssim = ssim, 
				probsim = probsim)
	}
	return(ans)
}

setMethod("starsim", signature(fit = "STARfit"), .starsim)


#------------------------------------------------------------------------------------
# path simulation method
#------------------------------------------------------------------------------------
starpath = function(spec, n.sim = 1000, n.start = 0, m.sim = 1, presigma = NA, 
		prereturns = NA, preresiduals = NA, rseed = NA, 
		custom.dist = list(name = NA, distfit = NA), xregsim = NULL, vregsim = NULL, 
		ssim = NULL, probsim = NULL, ...)

{
	UseMethod("starpath")
}

.starpath = function(spec, n.sim = 1000, n.start = 0, m.sim = 1, presigma = NA, 
		prereturns = NA, preresiduals = NA, rseed = NA, 
		custom.dist = list(name = NA, distfit = NA), xregsim = NULL, 
		vregsim= NULL, ssim = NULL, probsim = NULL)
{
	if(spec@model$modelinc[50]>0){
		if(spec@model$modelinc[50]==4){
			ans = .starpath.mixture(spec = spec, n.sim = n.sim, n.start = n.start, m.sim = m.sim,
					prereturns = prereturns, preresiduals = preresiduals, rseed = rseed,
					custom.dist = custom.dist, xregsim = xregsim, ssim = ssim, 
					probsim = probsim)
		} else{
			ans = .starpath.dynamic(spec = spec, n.sim = n.sim, n.start = n.start, m.sim = m.sim,
					presigma = presigma, prereturns = prereturns, preresiduals = preresiduals, 
					rseed = rseed, custom.dist = custom.dist, xregsim = xregsim, 
					vregsim = vregsim, ssim = ssim, probsim = probsim)
		}
	} else{
		ans = .starpath.static(spec = spec, n.sim = n.sim, n.start = n.start, m.sim = m.sim, 
				prereturns = prereturns, preresiduals = preresiduals, rseed = rseed, 
				custom.dist = custom.dist, xregsim = xregsim, ssim = ssim, 
				probsim = probsim)
	}
	return(ans)
}

setMethod("starpath", signature(spec = "STARspec"), .starpath)

#------------------------------------------------------------------------------------
# rolling backtest method
#------------------------------------------------------------------------------------
rollstar = function(spec, data, n.ahead = 1, forecast.length = 500, 
		n.start = NULL, refit.every = 25, refit.window = c("recursive", "moving"), 
		window.size = NULL, solver = "msoptim", fit.control = list(), solver.control = list(),
		calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05), cluster = NULL,
		keep.coef = TRUE, ...)
{
	UseMethod("rollstar")
}

setMethod("rollstar", signature(spec = "STARspec"),  definition = .rollstar)

#------------------------------------------------------------------------------------
# residuals method
#------------------------------------------------------------------------------------
.starresids = function(object, standardize = FALSE)
{
	if(class(object)[1] == "STARfit" | class(object)[1] == "STARfilter"){
		D = object@model$modeldata$index[1:object@model$modeldata$T]
	}
	ans = switch(class(object)[1],
			STARfit = xts(object@fit$residuals, D),
			STARfilter = xts(object@filter$residuals, D))
	if(standardize) ans = ans/sigma(object)
	return(ans)
}

setMethod("residuals", signature(object = "STARfit"), .starresids)
setMethod("residuals", signature(object = "STARfilter"), .starresids)

#------------------------------------------------------------------------------------
# fitted method
#------------------------------------------------------------------------------------
.starfitted = function(object)
{
	if(class(object)[1] == "STARfit" | class(object)[1] == "STARfilter"){
		D = object@model$modeldata$index[1:object@model$modeldata$T]
	}
	switch(class(object)[1],
			STARfit = xts(object@fit$fitted.values, D), 
			STARfilter = xts(object@model$modeldata$data[1:object@model$modeldata$T] - object@filter$residuals, D),
			STARsim = {
				ans = object@simulation$seriesSim
				rownames(ans) = paste("T+",1:NROW(object@simulation$seriesSim), sep="")
				return(ans)
			},
			STARpath ={
				ans = object@path$seriesSim
				rownames(ans) = paste("T+",1:NROW(object@path$seriesSim), sep="")
				return(ans)
			},
			STARforecast = object@forecast$seriesFor
	)
}
setMethod("fitted", signature(object = "STARfit"), .starfitted)
setMethod("fitted", signature(object = "STARfilter"), .starfitted)
setMethod("fitted", signature(object = "STARsim"), .starfitted)
setMethod("fitted", signature(object = "STARpath"), .starfitted)
setMethod("fitted", signature(object = "STARforecast"), .starfitted)

#------------------------------------------------------------------------------------
# states (probabilities and dynamics)
#------------------------------------------------------------------------------------
states = function(object, ...)
{
	UseMethod("states")
}

.starstate = function(object, type="prob")
{
	if(class(object)[1] == "STARfit" | class(object)[1] == "STARfilter"){
		D = object@model$modeldata$index[1:object@model$modeldata$T]
	}
	if(type=="prob"){
		ans = switch(class(object)[1],
				STARfit = xts(object@fit$probability[1:object@model$modeldata$T,], D),
				STARfilter = xts(object@filter$probability[1:object@model$modeldata$T,], D))
		colnames(ans)<-paste("state[",1:ncol(ans),"]",sep="")
	} else if(type=="condm"){
		ans = switch(class(object)[1],
				STARfit = xts(object@fit$condm[1:object@model$modeldata$T,], D),
				STARfilter = xts(object@filter$condm[1:object@model$modeldata$T,], D))
		colnames(ans)<-paste("state[",1:ncol(ans),"]",sep="")
	} else{
		if(object@model$modelinc[46]>1 & object@model$modelinc[47]==0){
			ans = switch(class(object)[1],
					STARfit = xts(object@fit$pmu[1:object@model$modeldata$T,], D),
					STARfilter = xts(object@filter$pmu[1:object@model$modeldata$T,], D))			
		} else{
			stop("\nNo dynamics for fixed.probs or 1-state model.")
		}
	}
	return(ans)
}

setMethod("states", signature(object = "STARfit"), .starstate)
setMethod("states", signature(object = "STARfilter"), .starstate)

.starstatef = function(object)
{
	P = object@forecast$probFor
	if(object@forecast$method=="n.ahead-1"){
		f = matrix(P, ncol = dim(P)[3])
		rownames(f)<-colnames(P[,,1,drop=FALSE])
		colnames(f)<-paste("State_",1:dim(P)[3],"[T+1]",sep="")
	} else{
		f = P
	}
	return(f)
}
setMethod("states", signature(object = "STARforecast"), .starstatef)


.starstatesim = function(object, type="prob", sim=1)
{
	sim = as.integer(max(sim,1))
	m.sim = ncol(object@simulation$seriesSim)
	if(sim>m.sim) stop("\ntwinkle-->error: sim>m.sim!")
	if(type=="prob"){
		ans = object@simulation$probSim[[sim]]
		colnames(ans)<-paste("state[",1:ncol(ans),"]",sep="")
		rownames(ans) = paste("T+",1:NROW(ans), sep="")
	} else if(type=="condm"){
		ans = object@simulation$condmSim[[sim]]
		colnames(ans)<-paste("state[",1:ncol(ans),"]",sep="")
		rownames(ans) = paste("T+",1:NROW(ans), sep="")
	} else{
			stop("\nOnly prob and condm type allowed for simulation objects.")
	}
	return(ans)
}

.starstatepath = function(object, type="prob", sim=1)
{
	sim = as.integer(max(sim,1))
	m.sim = ncol(object@path$seriesSim)
	if(sim>m.sim) stop("\ntwinkle-->error: sim>m.sim!")
	if(type=="prob"){
		ans = object@path$probSim[[sim]]
		colnames(ans)<-paste("state[",1:ncol(ans),"]",sep="")
		rownames(ans) = paste("T+",1:NROW(ans), sep="")
	} else if(type=="condm"){
		ans = object@path$condmSim[[sim]]
		colnames(ans)<-paste("state[",1:ncol(ans),"]",sep="")
		rownames(ans) = paste("T+",1:NROW(ans), sep="")
	} else{
		stop("\nOnly prob and condm type allowed for simulation objects.")
	}
	return(ans)
}
setMethod("states", signature(object = "STARsim"), .starstatesim)
setMethod("states", signature(object = "STARpath"), .starstatepath)



#------------------------------------------------------------------------------------
# coef
#------------------------------------------------------------------------------------
.starfitcoef = function(object)
{
	ans = switch(class(object)[1],
			STARfit = object@model$pars[object@model$pars[,4]==1, 1],
			STARfilter = object@model$pars[object@model$pars[,2]==1, 1])
	return(ans)
}

setMethod("coef", signature(object = "STARfit"), .starfitcoef)
setMethod("coef", signature(object = "STARfilter"), .starfitcoef)

.starrollcoef = function(object)
{
	if(!is.null(object@model$noncidx)) stop("\nObject contains non-converged estimation windows.")
	return(object@model$coef) 
}

setMethod("coef", signature(object = "STARroll"), .starrollcoef)

#------------------------------------------------------------------------------------
# Log-Likelihood
.starLikelihood = function(object)
{
	switch(class(object)[1],
			STARfit = object@fit$LLH,
			STARfilter = object@filter$LLH)
}

setMethod("likelihood", signature(object = "STARfit"), .starLikelihood)
setMethod("likelihood", signature(object = "STARfilter"), .starLikelihood)

#------------------------------------------------------------------------------------
# Information Criteria (standardized for data size)
#------------------------------------------------------------------------------------
.starinfocriteria = function(object)
{
	# indicator object@fit$ipars[,4] denotes the estimated parameters
	if(is(object, "STARfilter")){
		np = 0
	} else{
		np = sum(object@fit$ipars[,4])
	}
	itest = .information.test(likelihood(object), nObs = length(fitted(object)), 
			nPars = np)
	itestm = matrix(0, ncol = 1, nrow = 4)
	itestm[1,1] = itest$AIC
	itestm[2,1] = itest$BIC
	itestm[3,1] = itest$SIC
	itestm[4,1] = itest$HQIC
	colnames(itestm) = ""
	rownames(itestm) = c("Akaike", "Bayes", "Shibata", "Hannan-Quinn")
	return(itestm)
}

setMethod("infocriteria", signature(object = "STARfit"), .starinfocriteria)
setMethod("infocriteria", signature(object = "STARfilter"), .starinfocriteria)


.rollstarinfocriteria = function(object)
{
	# indicator object@fit$ipars[,4] denotes the estimated parameters
	np = as.numeric(sum(object@model$spec@model$pars[,4]))
	lik = as.numeric(object@model$loglik)
	n = length(lik)
	N = as.numeric(sapply(object@model$rollind, "length") - object@model$out.sample)
	itestm = sapply(1:n, function(i) unlist(.information.test(lik[i], 
						nObs = N[i], nPars = np)))
	colnames(itestm) = paste("window_",1:n, sep="")
	rownames(itestm) = c("Akaike", "Bayes", "Shibata", "Hannan-Quinn")
	return(itestm)
}

setMethod("infocriteria", signature(object = "STARroll"), .rollstarinfocriteria)

#------------------------------------------------------------------------------------
# vcov
.vcovstar = function(object, robust = FALSE)
{
	if(robust){
		vc = object@fit$robust.cvar
	} else{
		vc = object@fit$cvar
	}
	colnames(vc) = rownames(vc) = names(object@model$pars[object@model$pars[,4]==1,1])
	return(vc)
}

setMethod("vcov", signature(object = "STARfit"), .vcovstar)


#------------------------------------------------------------------------------------
# solver convergence
#------------------------------------------------------------------------------------
.starconvergence = function(object){
	return( object@fit$convergence )
}
setMethod("convergence", signature(object = "STARfit"),  definition = .starconvergence)

.convergenceroll = function(object){
	nonc = object@model$noncidx
	if(is.null(nonc)){
		ans = 0 
	} else{
		ans = 1
		attr(ans, 'nonconverged')<-nonc
	}
	return(ans)
}
setMethod("convergence", signature(object = "STARroll"),  definition = .convergenceroll)
#------------------------------------------------------------------------------------
# score (jacobian)
#------------------------------------------------------------------------------------
score = function(object, ... ) { UseMethod("score") }

.starscore = function(object)
{
	return(object@fit$scores)
}
setMethod("score", signature(object = "STARfit"),  definition = .starscore)

#------------------------------------------------------------------------------------
# S4 model.matrix
#------------------------------------------------------------------------------------
# ToDo: provide option for returning the model.matrix of the state equation
# together with a function to evaluate it.
setOldClass("xts")

modelmatrix = function(object, data, linear = TRUE, ... ) { UseMethod("modelmatrix") }

.modelmatrix = function(object, data, linear = TRUE, ...)
{
	model = object@model
	modelinc = model$modelinc
	inc = 0
	if(is(object, "STARspec")){
		if(is.null(data)) stop("\nmodelmatrix error:--> you need to supply an xts data series when using modelmatrix with a STARspec object")
		if(!is.xts(data)) stop("\nmodelmatrix error:--> you need to supply a xts data series when using modelmatrix with a STARspec object")
		y = coredata(data)[,1]
		T = length(y)
		yindex = index(data)
		res = rep(0, T)
	} else{
		T = model$modeldata$T
		y = model$modeldata$data[1:T]
		yindex = model$modeldata$index[1:T]
		res = as.numeric(residuals(object))
	}
	if(linear){
		# Linear part
		r = modelinc[46]
		iconst = modelinc[c(1,5,9,13)]
		ipars= model$pars
		idx  = model$pidx
		iar  = modelinc[c(2,6,10,14)]
		ireg = modelinc[c(3,7,11,15)]
		ima  = modelinc[c(4,8,12,16)]
		lma = modelinc[17]
		xreg = coredata(model$modeldata$mexdata[1:T,,drop=FALSE])
		X = NULL
		cnames = NULL
		for(i in 1:r){
			if(iconst[i]>0){
				X = cbind(X, matrix(1, nrow=T, ncol=1))
				cp = ipars[idx["s1.phi0",1],1]
				cnames  = c(cnames, paste("s",i,".phi0",sep=""))						
			} else{
				cp = 0
			}
			if(iar[i]>0){
				for(j in 1:iar[i]){
					X = cbind(X, lagf.numeric(y, j)-inc*modelinc[46]*cp)
					cnames  = c(cnames, paste("s",i,".phi",j,sep=""))
				}
				X[is.na(X)] = 0
			}
			if(ireg[i]>0){
				X = cbind(X, xreg)
				cnames  = c(cnames, paste("s",i,".xi",1:ncol(xreg),sep=""))
			}
			if(ima[i]>0){
				for(j in 1:ima[i]){
					X = cbind(X, lagf.numeric(res, j))
					cnames  = c(cnames, paste("s",i,".psi",j,sep=""))
				}
				X[is.na(X)] = 0
			}
		}
		if(lma>0){
			for(j in 1:lma){
				X = cbind(X, lagf.numeric(res, j))
				cnames  = c(cnames, paste(".psi",j,sep=""))
			}
		}
		colnames(X) = cnames
		rownames(X) = as.character(yindex)
	} else{
		# State part
		if(modelinc[47]==0){
			if(modelinc[49]==2){
				chk = all(yindex==index(model$modeldata$s))
				if(!is.logical(chk) | chk == FALSE){
					print(paste("\n",chk,sep=""))
					stop("\ntwinkle-->error: data and s indices do not match\n")
				}
				X = model$modeldata$s
			} else{
				# we apply yfun here if it is needed for efficiency:
				if(modelinc[48]==1){
					ytmp = xts(model$modeldata$fun(y), yindex)
				} else{
					ytmp = xts(y, yindex)
				}
				X = NULL
				for(i in 1:length(model$modeldata$ylags)){
					if(i==1) X = lag(ytmp, model$modeldata$ylags[i]) else X = cbind(X, lag(ytmp, model$modeldata$ylags[i]))
				}
				X[is.na(X)]=0
				colnames(X) = paste("s[t-",model$modeldata$ylags,"]",sep="")
			}
			# convert to matrix
		} else{
			warning("\nmodel based on fixed probabilities has not state modelmatrix!")
			X = NULL
		}
	}
	return(X)
}
setMethod("modelmatrix", signature(object = "STARfit", data="missing"),  definition = .modelmatrix)
setMethod("modelmatrix", signature(object = "STARspec", data="xts"),  definition = .modelmatrix)

#------------------------------------------------------------------------------------
# sigma (GARCH)
#------------------------------------------------------------------------------------
# for 2 state return the blended sigma
.starsigma = function(object)
{
	if(object@model$modelinc[50]>0){
		if(class(object)[1] == "STARfit" | class(object)[1] == "STARfilter"){
			D = object@model$modeldata$index[1:object@model$modeldata$T]
		}
		switch(class(object)[1],
				STARfit = xts(object@fit$sigma, D),
				STARfilter = xts(object@filter$sigma, D),
				STARsim = {
					ans = object@simulation$sigmaSim
					rownames(ans) = paste("T+",1:NROW(object@simulation$sigmaSim), sep="")
					return(ans)
				},
				STARpath = {
					ans = object@path$sigmaSim
					rownames(ans) = paste("T+",1:NROW(object@path$sigmaSim), sep="")
					return(ans)
				},
				STARforecast = object@forecast$sigmaFor
		)
	} else{
		object@model$pars["sigma",1]
	}
}

setMethod("sigma", signature(object = "STARfit"), .starsigma)
setMethod("sigma", signature(object = "STARfilter"), .starsigma)
setMethod("sigma", signature(object = "STARsim"), .starsigma)
setMethod("sigma", signature(object = "STARpath"), .starsigma)
setMethod("sigma", signature(object = "STARforecast"), .starsigma)

#------------------------------------------------------------------------------------
# as.data.frame (for roll object)
#------------------------------------------------------------------------------------
.starrolldf = function(x, row.names = NULL, optional = FALSE, which = "density")
{
	if(!is.null(x@model$noncidx)) stop("\nObject contains non-converged estimation windows.")
	if(which == "density") ans =  x@forecast$density else ans = x@forecast$VaR
	return(ans)
}
setMethod("as.data.frame", signature(x = "STARroll"), .starrolldf)

#------------------------------------------------------------------------------------
# resume rollstar
#------------------------------------------------------------------------------------

setMethod("resume", signature(object = "STARroll"),  definition = .resumerollstar)


#------------------------------------------------------------------------------------
# quantile
#------------------------------------------------------------------------------------
.starquantile = function(x, probs=c(0.01, 0.05))
{	
	if(class(x)=="STARroll"){
		d = x@model$spec@model$modeldesc$distribution
	} else{
		d = x@model$modeldesc$distribution
		s = sigma(x)
		m = fitted(x)
	}
	di = .DistributionBounds(d)
	if(di$include.skew)  skew  = x@model$pars["skew",1] else skew = 0
	if(di$include.shape) shape = x@model$pars["shape",1] else shape = 0
	if(di$include.ghlambda) ghlambda = x@model$pars["ghlambda",1] else ghlambda = 0
	
	if(class(x)=="STARforecast"){
		ans = matrix(NA, dim(m)[1], dim(m)[2])
		if(length(probs)>1) stop("\nprobs must be a scalar for a STARforecast object")
		for(i in 1:NCOL(ans)) ans[,i] = qdist(d, probs[1], mu = m[,i], if(x@model$modelinc[30]==0) sigma = s[,i] else sigma = s, 
					skew = skew, shape = shape, lambda = ghlambda)
		colnames(ans) = colnames(m)
		rownames(ans) = rownames(m)
	} else if(class(x)=="STARsim"){
		if(length(probs)>1) stop("\nprobs must be a scalar for a STARsim object\n")
		ans = matrix(NA, dim(m)[1], dim(m)[2])
		for(i in 1:NCOL(ans)) ans[,i] = qdist(d, probs[1], mu = m[,i], if(x@model$modelinc[30]==0) sigma = s[,i] else sigma = s, 
					skew = skew, shape = shape, lambda = ghlambda)
		colnames(ans) = colnames(m)
		rownames(ans) = rownames(m)
	} else if(class(x)=="STARpath"){
		if(length(probs)>1) stop("\nprobs must be a scalar for a STARpath object\n")
		ans = matrix(NA, dim(m)[1], dim(m)[2])
		for(i in 1:NCOL(ans)) ans[,i] = qdist(d, probs[1], mu = m[,i], if(x@model$modelinc[30]==0) sigma = s[,i] else sigma = s, 
					skew = skew, shape = shape, lambda = ghlambda)
		colnames(ans) = colnames(m)
		rownames(ans) = rownames(m)
	} else if(class(x)=="STARroll"){
		skew = x@forecast$density[,"Skew"]
		shape = x@forecast$density[,"Shape"]
		lambda = x@forecast$density[,"Shape(GIG)"]
		s = x@forecast$density[,"Sigma"]
		m = x@forecast$density[,"Mu"]
		ans = matrix(NA, ncol = length(probs), nrow = length(m))
		for(i in 1:length(probs)) ans[,i] = as.numeric(m) + qdist(d, probs[i], skew = skew, shape = shape, 
					lambda = lambda) * as.numeric(s)
		colnames(ans) = paste("q[", probs,"]", sep="")
		ans = xts(ans, as.POSIXct(rownames(x@forecast$density)))
	} else{
		ans = matrix(NA, ncol = length(probs), nrow = NROW(m))
		for(i in 1:length(probs)) ans[,i] = as.numeric(m) + qdist(d, probs[i], skew = skew, shape = shape, 
					lambda = ghlambda) * as.numeric(s)
		colnames(ans) = paste("q[", probs,"]", sep="")
		ans = xts(ans, index(m))
	}
	return(ans)
}
setMethod("quantile", signature(x = "STARfit"), .starquantile)
setMethod("quantile", signature(x = "STARfilter"), .starquantile)
setMethod("quantile", signature(x = "STARforecast"), .starquantile)
setMethod("quantile", signature(x = "STARsim"), .starquantile)
setMethod("quantile", signature(x = "STARpath"), .starquantile)
setMethod("quantile", signature(x = "STARroll"), .starquantile)
#------------------------------------------------------------------------------------
# probability integral transformation (pit) method
#------------------------------------------------------------------------------------
.starpit = function(object)
{	
	if(class(object)=="STARroll"){
		d = object@model$spec@model$modeldesc$distribution
		skew = object@forecast$density[,"Skew"]
		shape = object@forecast$density[,"Shape"]
		lambda = object@forecast$density[,"Shape(GIG)"]
		s = object@forecast$density[,"Sigma"]
		m = object@forecast$density[,"Mu"]
		r = object@forecast$density[,"Realized"]
		ans =  pdist(d, q = r, mu = as.numeric(m), sigma = as.numeric(s), 
				skew = skew, shape = shape, lambda = ghlambda)
		ans = xts(ans, as.POSIXct(rownames(object@forecast$density)))
		colnames(ans) = "pit"
	} else{
		d = object@model$modeldesc$distribution
		di = .DistributionBounds(d)
		if(di$include.skew)  skew  = object@model$pars["skew",1] else skew = 0
		if(di$include.shape) shape = object@model$pars["shape",1] else shape = 0
		if(di$include.ghlambda) ghlambda = object@model$pars["ghlambda",1] else ghlambda = 0
		s = sigma(object)
		m = fitted(object)
		r = object@model$modeldata$data[1:object@model$modeldata$T]
		ans =  pdist(d, q = r, mu = as.numeric(m), sigma = as.numeric(s), 
				skew = skew, shape = shape, lambda = ghlambda)
		ans = xts(ans, index(s))
		colnames(ans) = "pit"
	}
	return(ans)
}
setMethod("pit", signature(object = "STARfit"), .starpit)
setMethod("pit", signature(object = "STARfilter"), .starpit)
setMethod("pit", signature(object = "STARroll"), .starpit)
#------------------------------------------------------------------------------------
# Show method
#------------------------------------------------------------------------------------
setMethod("show",
		signature(object = "STARspec"),
		function(object){
			vmodel = object@model$modeldesc$vmodel
			model = object@model
			modelinc = object@model$modelinc
			dist = c("norm", "snorm", "std", "sstd","ged", "sged", "nig", "ghyp", "jsu", "ghst")[modelinc[44]]
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\n*          STAR Model Spec        *", sep = ""))
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\nstates       : ", modelinc[46], sep = ""))
			if(modelinc[46]>1) cat(paste("\nstatevar     : ", c("y","s")[modelinc[49]], sep = ""))
			if(modelinc[46]>1) cat(paste("\nstatear      : ", as.logical(modelinc[21]), sep = ""))
			cnt = modelinc[c(1,5,9,13)]
			cat(paste("\nconstant     :"), "[",as.logical(cnt[1:modelinc[46]]),"]")
			ar = modelinc[c(1,5,9,13)+1]
			cat(paste("\narOrder      :"), "[",as.integer(ar[1:modelinc[46]]),"]")
			ma = modelinc[c(1,5,9,13)+3]
			if(modelinc[17]>0){
			cat(paste("\nmaOrder (Linear) : ", modelinc[17], sep = ""))
			} else{
				if(sum(ma)>0) cat(paste("\nmaOrder (state) :"), "[",as.integer(ma[1:modelinc[46]]),"]")
			}
			cat(paste("\nregressors   : ", modelinc[3], sep=""))
			vty = c("dynamic","mixture","static")[ifelse(modelinc[50]>0, ifelse(modelinc[50]==4,2,1),3)]
			cat(paste("\nvariance     : ", vty, sep=""))
			cat(paste("\ndistribution : ", dist, "\n",sep = ""))
			invisible(object)
		})
setMethod("show",
		signature(object = "STARfit"),
		function(object){
			vmodel = object@model$modeldesc$vmodel
			model = object@model
			modelinc = object@model$modelinc
			dist = c("norm", "snorm", "std", "sstd","ged", "sged", "nig", "ghyp", "jsu", "ghst")[modelinc[44]]
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\n*          STAR Model Fit         *", sep = ""))
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\nstates       : ", modelinc[46], sep = ""))
			cat(paste("\nstatevar     : ", c("y","s")[modelinc[49]], sep = ""))
			cat(paste("\nstatear      : ", as.logical(modelinc[21]), sep = ""))
			vty = c("dynamic","mixture","static")[ifelse(modelinc[50]>0, ifelse(modelinc[50]==4,2,1),3)]
			cat(paste("\nvariance     : ", vty, sep=""))
			cat(paste("\ndistribution : ", dist, sep = ""))
			if(convergence(object)==0){
				cat("\n\n")
				cat("\nOptimal Parameters (Robust Standard Errors)")
				cat(paste("\n------------------------------------\n",sep=""))
				print(round(object@fit$robust.matcoef,6), digits = 5)
				cat("\nLogLikelihood :", likelihood(object), "\n")
				itestm = infocriteria(object)
				print(itestm,digits=5)
				res = as.numeric(residuals(object))
				rss = sum(res^2)
				N = object@model$modeldata$T
				tot = sum(scale(object@model$modeldata$data[1:N], scale=FALSE)^2)
				Rsq = 1 - rss/tot
				npars = sum(object@fit$ipars[,4])
				Rsqadj = Rsq - (1-Rsq)*(npars/(N-npars-1))
				cat("\nr.squared         : ", round(Rsq,4))
				cat("\nr.squared (adj)   : ", round(Rsqadj,4))
				cat("\nRSS               : ", round(rss, 5))
				cat("\nskewness (res)    : ", round(.skewness(res), 5))
				cat("\nex.kurtosis (res) : ", round(.kurtosis(res), 5))
				armat = ar_root(object)
				if(armat$use){
					cat("\n")
					cat("\nAR roots\n")
					print(armat$rmat)
				}
				mamat = ma_root(object)
				if(mamat$use){
					cat("\n")
					cat("\nMA roots\n")
					print(mamat$rmat)
				}
				cat("\n")
				cat("\n")
			}
			invisible(object)
		})

setMethod("show",
		signature(object = "STARfilter"),
		function(object){
			vmodel = object@model$modeldesc$vmodel
			model = object@model
			modelinc = object@model$modelinc
			dist = c("norm", "snorm", "std", "sstd","ged", "sged", "nig", "ghyp", "jsu", "ghst")[modelinc[44]]
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\n*          STAR Model Fit         *", sep = ""))
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\nstates       : ", modelinc[46], sep = ""))
			cat(paste("\nstatevar     : ", c("y","s")[modelinc[49]], sep = ""))
			cat(paste("\nstatear      : ", as.logical(modelinc[21]), sep = ""))
			vty = c("dynamic","mixture","static")[ifelse(modelinc[50]>0, ifelse(modelinc[50]==4,2,1),3)]
			cat(paste("\nvariance     : ", vty, sep=""))
			cat(paste("\ndistribution : ", dist, sep = ""))
			cat("\n\n")
			cat("\nFilter Parameters")
			cat(paste("\n---------------------------------------\n",sep=""))
			print(matrix(coef(object), ncol=1, dimnames = list(names(coef(object)), "")), digits = 5)
			cat("\nLogLikelihood :", object@filter$LLH, "\n")
			itestm = infocriteria(object)
			cat("\nInformation Criteria")
			cat(paste("\n---------------------------------------\n",sep=""))
			print(itestm,digits=5)
			invisible(object)
		})


setMethod("show",
		signature(object = "STARforecast"),
		function(object){
			model = object@model
			cat(paste("\n*------------------------------------*", sep = ""))
			cat(paste("\n*        STAR Model Forecast         *", sep = ""))
			cat(paste("\n*------------------------------------*", sep = ""))
			n.ahead = object@forecast$n.ahead
			cat(paste("\nHorizon        : ", n.ahead, sep = ""))
			cat(paste("\nRoll Steps     : ",object@forecast$n.roll, sep = ""))
			cat(paste("\nSTAR forecast  : ",object@forecast$method, sep = ""))
			n.start = object@forecast$n.start
			if(n.start>0) infor = ifelse(n.ahead>n.start, n.start, n.ahead) else infor = 0
			cat(paste("\nOut of Sample  : ", infor, "\n", sep = ""))
			cat(paste("\n0-roll forecast [T0=", as.character(object@model$modeldata$index[object@model$modeldata$T]), "]:\n", sep=""))
			if(model$modelinc[47]>0){
				zz = cbind(object@forecast$seriesFor[,1], object@forecast$sigmaFor[,1])
				colnames(zz) = c("Series", "Sigma")
				rownames(zz) = paste("T+",1:NROW(object@forecast$seriesFor), sep="")
			} else{
				zz = object@forecast$seriesFor[,1,drop=FALSE]
				colnames(zz) = c("Series")
				rownames(zz) = paste("T+",1:NROW(object@forecast$seriesFor), sep="")
			}
			print(zz, digits = 4)
			cat("\n\n")
		})
setMethod("show",
		signature(object = "STARroll"),
		function(object){
			if(!is.null(object@model$noncidx)){
				cat("\nObject contains non-converged estimation windows. Use resume method to re-estimate.\n")
				invisible(object)
			} else{
				cat(paste("\n*-------------------------------------*", sep = ""))
				cat(paste("\n*              STAR Roll              *", sep = ""))
				cat(paste("\n*-------------------------------------*", sep = ""))
				N = object@model$n.refits
				model = object@model$spec@model
				modelinc = model$modelinc
				vmodel = model$modeldesc$vmodel
				dist = c("norm", "snorm", "std", "sstd","ged", "sged", "nig", "ghyp", "jsu", "ghst")[modelinc[44]]
				cat("\nNo.Refits\t\t:", N)
				cat("\nRefit Horizon\t:", object@model$refit.every)
				cat("\nNo.Forecasts\t:", NROW(object@forecast$density))
				cat(paste("\nstates      : ", modelinc[46], sep = ""))
				cat(paste("\nstatevar     : ", c("y","s")[modelinc[49]], sep = ""))
				cat(paste("\nstatear      : ", as.logical(modelinc[21]), sep = ""))
				vty = c("dynamic","mixture","static")[ifelse(modelinc[50]>0, ifelse(modelinc[50]==4,2,1),3)]
				cat(paste("\nvariance     : ", vty, sep=""))
				cat(paste("\ndistribution : ", dist, sep = ""))
				cat("\n")
				cat("\nForecast Density\n")
				print(round(head(object@forecast$density),4))
				cat("\n..........................\n")
				print(round(tail(object@forecast$density),4))
				cat("\nElapsed:", format(object@model$elapsed))
				cat("\n")				
				invisible(object)
			}
		})

setMethod("show",
		signature(object = "STARsim"),
		function(object){
			model = object@model
			cat(paste("\n*--------------------------------------*", sep = ""))
			cat(paste("\n*          STAR Simulation             *", sep = ""))
			cat(paste("\n*--------------------------------------*", sep = ""))
			sim = object@simulation
			sigma = sim$sigmaSim
			series = sim$seriesSim
			resids = sim$residSim
			m = dim(series)[2]
			N = dim(series)[1]
			cat(paste("\nHorizon: ",N))
			cat(paste("\nSimulations: ",m,"\n",sep=""))
			cat(paste("\nSummary[m.sim=1]:\n"))
			Z1 = cbind(fitted(object), states(object))
			colnames(Z1)[1] = "series"
			Z2 = cbind(fitted(object), states(object))
			colnames(Z2) = NULL
			print(head(Z1))
			cat("\n.......\n")
			print(tail(Z2))
			cat("\n\n")
		})

setMethod("show",
		signature(object = "STARpath"),
		function(object){
			model = object@model
			cat(paste("\n*--------------------------------------*", sep = ""))
			cat(paste("\n*          STAR Path Simulation        *", sep = ""))
			cat(paste("\n*--------------------------------------*", sep = ""))
			sim = object@path
			sigma = sim$sigmaSim
			series = sim$seriesSim
			resids = sim$residSim
			m = dim(series)[2]
			N = dim(series)[1]
			cat(paste("\nHorizon: ",N))
			cat(paste("\nSimulations: ",m,"\n",sep=""))
			cat(paste("\nSummary[m.sim=1]:\n"))
			Z1 = cbind(fitted(object), states(object))
			colnames(Z1)[1] = "series"
			Z2 = cbind(fitted(object), states(object))
			colnames(Z2) = NULL
			print(head(Z1))
			cat("\n.......\n")
			print(tail(Z2))
			cat("\n\n")
		})
#----------------------------------------------------------------------------------
# plot methods
#----------------------------------------------------------------------------------
setMethod(f = "plot", signature(x = "STARfit", y = "missing"), plot.starfit)
setMethod(f = "plot", signature(x = "STARforecast", y = "missing"), plot.starforecast)


#----------------------------------------------------------------------------------
# Tests
#----------------------------------------------------------------------------------

nonlinearTest = function(object, data, robust = FALSE, sig.level = 0.05, ...)
{
	UseMethod("linearTest")
}

setMethod(f = "nonlinearTest", signature(object = "STARfit",  data = "missing"), nonlinearTest.fit)
setMethod(f = "nonlinearTest", signature(object = "STARspec", data = "xts"), nonlinearTest.spec)