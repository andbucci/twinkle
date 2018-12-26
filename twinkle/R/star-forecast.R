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

# To avoid confusion and problems, the user is required to supply
# data for external forecast variables from (T-out.sample+1):(T-out.sample+n.ahead+n.roll)
# even if the spec object contains T data for these regressors.
.forcregressors = function(model, xregfor, vregfor, sfor, probfor, n.ahead, N, out.sample, n.roll)
{	
	# N is the original length
	treq = n.ahead + n.roll
	modelinc = model$modelinc
	mxn = modelinc[3]
	vxn = modelinc[39]
	sxn = modelinc[20]
	if(mxn>0){
		if(!is.null(xregfor)){
			nmex = NROW(as.matrix(xregfor))
			mmex = NCOL(as.matrix(xregfor))
		} else{
			nmex = 0
			mmex = 0
		}
		if(!is.null(xregfor) && mmex != mxn)
		{
			cat("\nstarforecast-->error: Column dimension of external mean forecast matrix is wrong.")
			cat(paste("\nModel has ", mxn, " external regressors but forecast matrix has ", mmex, sep = ""))
			stop("\n...exiting\n")
		}
		if(!is.null(xregfor) && nmex < treq)
		{
			cat(paste("\nstarforecast-->error: You requested ", treq ," actual forecasts (including the rolling periods) but external mean forecasts provided have only ", nmex, " rows",sep=""))
			cat(paste("\nA minimum of ", treq," rows are required", sep=""))
			stop("\n...exiting")
		}
		if(is.null(xregfor)){
			# if NULL, augment with what is available from specification object
			if(out.sample>=treq){
				mxf = as.matrix(model$modeldata$mexdata)[1:(N-out.sample+treq), ,drop = FALSE]
			} else{
				mxf = rbind(as.matrix(model$modeldata$mexdata)[1:N, ,drop = FALSE], matrix(0, ncol = mxn, nrow = treq-out.sample))
			}
		} else {
			mxf = rbind(as.matrix(model$modeldata$mexdata)[1:(N-out.sample), ,drop = FALSE], as.matrix(xregfor)[1:treq,,drop = FALSE])
		}
	} else{
		mxf = NULL
	}
	if(vxn>0){
		if(!is.null(vregfor)){
			nvex = dim(as.matrix(vregfor))[1]
			mvex = dim(as.matrix(vregfor))[2]
		} else{
			nvex = 0
			mvex = 0
		}
		if(!is.null(vregfor) && mvex != vxn)
		{
			cat("\nstarforecast-->error: Column dimension of external variance forecast matrix is wrong.")
			cat(paste("\nModel has ",vxn," external regressors but forecast matrix has", mvex, sep = ""))
			stop("\n...exiting\n")
		}
		# N is the original length
		if(!is.null(vregfor) && nvex < treq)
		{
			cat(paste("\nstarforecast-->error: You requested ", treq ," actual forecasts (including the rolling periods) but external variance forecasts provided have only ",
							nvex, " rows",sep=""))
			cat(paste("\nA minimum of ", treq," rows are required", sep=""))
			stop("\n...exiting")
		}
		if(is.null(vregfor)){
			if(out.sample>=treq){
				vxf = as.matrix(model$modeldata$vexdata)[1:(N-out.sample+treq), ,drop = FALSE]
			} else{
				vxf = rbind(as.matrix(model$modeldata$vexdata)[1:N, ,drop = FALSE], matrix(0, ncol = vxn, nrow = treq-out.sample))
			}			
		} else {
			vxf = rbind(as.matrix(model$modeldata$vexdata)[1:(N-out.sample), ,drop = FALSE], as.matrix(vregfor)[1:treq,,drop = FALSE])
		}
	} else{
		vxf = NULL
	}
	if(modelinc[49]==2 && sxn>0){
		if(!is.null(sfor)){
			nxex = dim(as.matrix(sfor))[1]
			mxex = dim(as.matrix(sfor))[2]
		} else{
			nxex = 0
			mxex = 0
		}
		if(!is.null(sfor) && mxex != sxn)
		{
			cat("\nstarforecast-->error: Column dimension of 's' (probability dynamics) forecast matrix is wrong.")
			cat(paste("\nModel has ",sxn," regressors but forecast matrix provided has", mxex, sep = ""))
			stop("\n...exiting\n")
		}
		if(!is.null(sfor) && nxex < treq)
		{
			cat(paste("\nstarforecast-->error: You requested ", treq ," actual forecasts (including the rolling periods) but 's' (probability dynamics) forecasts provided have only ",nxex, " rows",sep=""))
			cat(paste("\nA minimum of ", treq," rows are required\n", sep=""))
			stop("...exiting", call. = FALSE)
		}
		if(is.null(sfor)){
			if(out.sample>=treq){
				sxf = as.matrix(model$modeldata$s)[1:(N-out.sample+treq), ,drop = FALSE]
			} else{
				stop("\nstarforecast-->error: if (n.roll+n.ahead)>out.sample, sfor cannot be NULL!\n.")
			}
		} else {
			sxf = rbind(coredata(model$modeldata$s)[1:(N-out.sample), ,drop = FALSE], as.matrix(sfor)[1:treq,,drop = FALSE])
		}
	} else{
		sxf = NULL
	}
	
	if(modelinc[47]>0){
		# flag as error since it means user probably needs to rethink what they are doing
		if(is.null(probfor)) stop("\nstarforecast-->error: probfor cannot be NULL in a system estimated with fixed state probabilities!\n.")
		npex = nrow(as.matrix(probfor))[1]
		mpex = ncol(as.matrix(probfor))[2]
		# columns==states
		if(mpex != modelinc[46])
		{
			cat("\nstarforecast-->error: Column dimension of 'probfor' matrix is wrong.")
			cat(paste("\nModel has ", modelinc[46]," states but probfor has", mpex, sep = ""))
			stop("\n...exiting\n")
		}
		# N is the original length
		# need to account for the lags in each
		if(npex < treq)
		{
			cat(paste("\nstarforecast-->error: You requested ", treq ," actual forecasts (including the rolling periods) but probfor provided have only ",
							npex, " rows",sep=""))
			cat(paste("\nA minimum of ", treq," rows are required", sep=""))
			stop("\n...exiting")
		}
		pxf = rbind(coredata(model$modeldata$fixed.prob)[1:(N-out.sample), ,drop = FALSE], as.matrix(probfor)[1:treq,,drop = FALSE])
	} else{
		pxf = NULL
	}
	return(list(mxf = mxf, vxf = vxf, sxf = sxf, pxf = pxf))
}

.nsgarchforecast = function(arglist)
{
	modelinc = arglist$modelinc
	n.ahead = arglist$n.ahead
	N = arglist$np
	idx = arglist$idx
	if(arglist$modelinc[39]>0){
		omega = arglist$omega + arglist$vxfi%*%t(matrix(arglist$ipars[idx["vxreg",1]:idx["vxreg",2],1], ncol = modelinc[39]))
	} else{
		omega = arglist$omega
	}
	for(i in 1:n.ahead){
		if(modelinc[33]>0){
			arglist$h[N+i] = omega[N+i] + sum(arglist$ipars[idx["beta",1]:idx["beta",2],1]*arglist$h[N+i-(1:modelinc[33])]^2)
		} else{
			arglist$h[N+i] = omega[N+i]
		}
		if(modelinc[32]>0){
			for (j in 1:modelinc[32]){
				if (i-j > 0){
					s = arglist$h[N + i - j]^2
				} else{ 
					s = arglist$epsx[N + i - j]^2
				}
				arglist$h[N+i] = arglist$h[N+i] + arglist$ipars[idx["alpha",1]+j-1,1] * s
			}
		}
		arglist$h[N+i] = sqrt(arglist$h[N+i])
	}
	arglist = starf(arglist)
	return(arglist)
}
.ngjrgarchforecast = function(arglist)
{
	modelinc = arglist$modelinc
	N = arglist$N
	n.ahead = arglist$n.ahead
	idx = arglist$idx
	if(modelinc[39]>0){
		omega = arglist$omega + arglist$vxfi%*%t(matrix(arglist$ipars[idx["vxreg",1]:idx["vxreg",2],1], ncol = modelinc[39]))
	} else{
		omega = arglist$omega
	}
	for(i in 1:n.ahead){
		if(modelinc[33]>0){
			arglist$h[N+i] = omega[N+i] + sum(arglist$ipars[idx["beta",1]:idx["beta",2],1]*arglist$h[N+i-(1:modelinc[33])]^2)
		} else{
			arglist$h[N+i] = omega[N+i]
		}
		if(modelinc[32]>0){
			for (j in 1:modelinc[32]){
				if (i-j > 0){				
					s1 = (arglist$h[N + i - j]^2)
					s2 = arglist$kappa*arglist$ipars[idx["gamma",1]+j-1, 1] * (arglist$h[N + i - j]^2)
				} else{
					s1 = arglist$epsx[N + i - j]^2
					s2 = arglist$ipars[idx["gamma",1]+j-1, 1]*(arglist$epsx[N + i - j]^2)*as.integer(arglist$epsx[N + i - j]<0)
				}
				arglist$h[N+i] = arglist$h[N+i] + arglist$ipars[idx["alpha",1]+j-1, 1] * s1 + s2
			}
		}
		arglist$h[N+i] = sqrt(arglist$h[N+i])
	}
	arglist = starf(arglist)
	return(arglist)
}

.negarchforecast = function(arglist)
{
	modelinc = arglist$modelinc
	N = arglist$N
	n.ahead = arglist$n.ahead
	idx = arglist$idx	
	if(modelinc[39]>0){
		omega = arglist$omega + arglist$vxfi%*%t(matrix(arglist$ipars[idx["vxreg",1]:idx["vxreg",2],1], ncol = modelinc[39]))
	} else{
		omega = arglist$omega
	}
	for(i in 1:n.ahead){
		if(modelinc[33]>0){
			arglist$h[N+i] = omega[N+i] + sum(arglist$ipars[idx["beta",1]:idx["beta",2],1]*log(arglist$h[N+i-(1:modelinc[33])]^2))
		} else{
			arglist$h[N+i] = omega[N+i]
		}
		if(modelinc[32]>0){
			for(j in 1:modelinc[32]){
				if (i-j > 0){
					s = 0 #unconditional
				} else{
					# first step ahead z is active else z=0
					s = arglist$ipars[idx["alpha",1]+j-1,1] * arglist$z[N+i-j] + arglist$ipars[idx["gamma",1]+j-1,1]*(abs(arglist$z[N+i-j]) - arglist$kappa)
				}
				arglist$h[N+i] = arglist$h[N+i] + s
			}
		}
		arglist$h[N+i] = sqrt(exp(arglist$h[N+i]))
	}
	# forecast STAR process
	arglist = starf(arglist)
	return(arglist)
}

# 1. probability forecast (XL)
# 2. external regressors (xreg) + mu = constm
# 3. y forecast

# need to account for MA forecast
starf = function(arglist)
{
	modelinc = arglist$modelinc
	# np = N + i(roll) -1
	N = arglist$np
	n.ahead = arglist$n.ahead
	idx = arglist$idx
	constm = rbind(arglist$constm, matrix(0, ncol = modelinc[46], nrow = n.ahead))
	if(modelinc[3]>0){
		for(i in 1:modelinc[46]){
			constm[(N+1):(N+n.ahead),i] = arglist$mxfi[(N+1):(N+n.ahead),]%*%t(matrix(arglist$ipars[idx[paste("s",i,".xi",sep=""),1]:idx[paste("s",i,".xi",sep=""),2],1], ncol = modelinc[3]))
		}
	} else{
		for(i in 1:modelinc[46]){
			constm[(N+1):(N+n.ahead),i] = arglist$ipars[idx[paste("s",i,".phi0",sep=""),1],1]
		}
	}
	dynamic = as.logical(modelinc[50]>0)
	arglist$constm = constm
	# n.ahead-1 since we are already passing the n.ahead=1 values
	arglist$pmu = rbind(arglist$pmu, matrix(0, ncol = max(1, modelinc[46]-1), n.ahead-1))
	arglist$probs = rbind(arglist$probs, matrix(0, ncol = modelinc[46], n.ahead-1))
	if(is.null(arglist$mc.sims)) mc.sims = 12*N else mc.sims = arglist$mc.sims
	ydist = matrix(NA, nrow = mc.sims, ncol = n.ahead)
	esim = matrix(NA, nrow = mc.sims, ncol = n.ahead)
	for(j in 2:n.ahead){
		arglist$j = j
		if(modelinc[47]==1){
			arglist$XL = NULL
			# probs fixed: populate forecast probs with pxfi
			arglist$probs[N+j,] = arglist$pxfi[N+j,]
		} else{
			if(modelinc[48]==1) ytmp = arglist$yfun(as.numeric(arglist$y[1:(N+j-1)])) else ytmp = as.numeric(arglist$y[1:(N+j-1)])
			arglist$XL = build.lagmatrix(modelinc, s = arglist$sxfi[1:(N+j),], y = c(ytmp, NA), arglist$ylags)
			ptmp = switch(modelinc[46], dstar1f(arglist), dstar2f(arglist), dstar3f(arglist), dstar4f(arglist))
			arglist$probs = ptmp$probs
			arglist$pmu = ptmp$pmu
		}
		if(arglist$method=="mc.empirical"){
			arglist$m = N+j
			if(dynamic){
				esim[,j] = sample(as.numeric(arglist$zresiduals[1:N]), mc.sims, replace=TRUE)*arglist$h[N+j]
			} else{
				esim[,j] = sample(as.numeric(arglist$residuals[1:N]), mc.sims, replace=TRUE)
			}
			tmp = starforc.mc(esim[,j], arglist)
			arglist$y[N+j] = mean(tmp)
			ydist[,j] = tmp
		} else if(arglist$method=="mc.kernel"){
			arglist$m = N+j
			if(dynamic) fit = .kfit(arglist$zresiduals[1:N]) else fit = .kfit(arglist$residuals[1:N])
			if(dynamic){
				esim[,j] = rkde(mc.sims, fit)*arglist$h[N+j]
			} else{
				esim[,j] = rkde(mc.sims, fit)
			}
			tmp = starforc.mc(esim[,j], arglist)
			arglist$y[N+j] = mean(tmp)
			ydist[,j] = tmp
		} else if(arglist$method=="mc.parametric"){
			arglist$m = N+j
			if(dynamic) arglist$sigma = arglist$h[N+j] else arglist$sigma = arglist$ipars["sigma",1]
			esim[,j] = rdist(arglist$distribution, mc.sims, 0, arglist$sigma, arglist$ipars["skew",1], arglist$ipars["shape",1],
					arglist$ipars["ghlambda",1])
			tmp = starforc.mc(esim[,j], arglist)
			arglist$y[N+j] = mean(tmp)
			ydist[,j] = tmp
		} else if(arglist$method=="an.kernel"){
			arglist$m = N+j
			if(dynamic) arglist$sigma = arglist$h[N+j] else arglist$sigma = arglist$ipars["sigma",1]
			if(dynamic) fit = .kfit(arglist$zresiduals[1:N]) else fit = .kfit(arglist$residuals[1:N])
			arglist$fit = fit
			if(dynamic) tmp = integrate(starforc.an3, -Inf, Inf, arglist = arglist) else tmp = integrate(starforc.an2, -Inf, Inf, arglist = arglist)
			arglist$y[N+j] = tmp$value
		} else{
			arglist$m = N+j
			if(dynamic) arglist$sigma = arglist$h[N+j] else arglist$sigma = arglist$ipars["sigma",1]
			arglist$skew = arglist$ipars["skew",1]
			arglist$shape = arglist$ipars["shape",1]
			arglist$ghlambda = arglist$ipars["ghlambda",1]
			arglist$epsilon = NULL
			tmp = integrate(starforc.an1, -Inf, Inf, arglist = arglist)
			arglist$y[N+j] = tmp$value
		}
	}
	arglist$ydist = ydist
	arglist$esim = esim
	return(arglist)
}
# QUESTION: Should we be simulating from the residuals of each state (depending on prob[N+j]) or just the residuals?
starfmix = function(arglist)
{
	# only 2-states
	modelinc = arglist$modelinc
	# np = N + i(roll) -1
	N = arglist$np
	n.ahead = arglist$n.ahead
	idx = arglist$idx
	constm = rbind(arglist$constm, matrix(0, ncol = modelinc[46], nrow = n.ahead))
	if(modelinc[3]>0){
		for(i in 1:modelinc[46]){
			constm[(N+1):(N+n.ahead),i] = arglist$mxfi[(N+1):(N+n.ahead),]%*%t(matrix(arglist$ipars[idx[paste("s",i,".xi",sep=""),1]:idx[paste("s",i,".xi",sep=""),2],1], ncol = modelinc[3]))
		}
	} else{
		for(i in 1:modelinc[46]){
			constm[(N+1):(N+n.ahead),i] = arglist$ipars[idx[paste("s",i,".phi0",sep=""),1],1]
		}
	}
	arglist$constm = constm
	# n.ahead-1 since we are already passing the n.ahead=1 values
	arglist$pmu = rbind(arglist$pmu, matrix(0, ncol = max(1, modelinc[46]-1), n.ahead-1))
	arglist$probs = rbind(arglist$probs, matrix(0, ncol = modelinc[46], n.ahead-1))
	if(is.null(arglist$mc.sims)) mc.sims = 12*N else mc.sims = arglist$mc.sims
	ydist = matrix(NA, nrow = mc.sims, ncol = n.ahead)
	esim = matrix(NA, nrow = mc.sims, ncol = n.ahead)
	arglist$dynamic = TRUE	
	for(j in 2:n.ahead){
		arglist$j = j
		if(modelinc[47]==1){
			arglist$XL = NULL
			# probs fixed: populate forecast probs with pxfi
			arglist$probs[N+j,] = arglist$pxfi[N+j,]
		} else{
			if(modelinc[48]==1) ytmp = arglist$yfun(as.numeric(arglist$y[1:(N+j-1)])) else ytmp = as.numeric(arglist$y[1:(N+j-1)]) 
			arglist$XL = build.lagmatrix(modelinc, s = arglist$sxfi[1:(N+j),], y = c(ytmp, NA), arglist$ylags)
			ptmp = switch(modelinc[46],
					dstar2f(arglist),
					dstar2f(arglist),
					dstar3f(arglist),
					dstar4f(arglist))
			arglist$probs = ptmp$probs
			arglist$pmu = ptmp$pmu
		}
		msig = switch(modelinc[46],
				arglist$probs[N+j,1]*arglist$ipars["s1.sigma",1]+arglist$probs[N+j,2]*arglist$ipars["s2.sigma",1],
				arglist$probs[N+j,1]*arglist$ipars["s1.sigma",1]+arglist$probs[N+j,2]*arglist$ipars["s2.sigma",1],
				arglist$probs[N+j,1]*arglist$ipars["s1.sigma",1]+arglist$probs[N+j,2]*arglist$ipars["s2.sigma",1]+arglist$probs[N+j,3]*arglist$ipars["s3.sigma",1],
				arglist$probs[N+j,1]*arglist$ipars["s1.sigma",1]+arglist$probs[N+j,2]*arglist$ipars["s2.sigma",1]+arglist$probs[N+j,3]*arglist$ipars["s3.sigma",1]+arglist$probs[N+j,4]*arglist$ipars["s4.sigma",1]
				)
		arglist$h[N+j] = arglist$sigma = msig
		if(arglist$method=="mc.empirical"){
			arglist$m = N+j
			esim[,j] = sample(as.numeric(arglist$zresiduals[1:N]), mc.sims, replace=TRUE)*arglist$h[N+j]
			tmp = starforc.mc(esim[,j], arglist)
			arglist$y[N+j] = mean(tmp)
			ydist[,j] = tmp
		} else if(arglist$method=="mc.kernel"){
			arglist$m = N+j
			fit = .kfit(arglist$zresiduals[1:N])
			esim[,j] = rkde(mc.sims, fit)*arglist$h[N+j]
			tmp = starforc.mc(esim[,j], arglist)
			arglist$y[N+j] = mean(tmp)
			ydist[,j] = tmp
		} else if(arglist$method=="mc.parametric"){
			arglist$m = N+j			
			esim[,j] = rdist(arglist$distribution, mc.sims, 0, arglist$h[N+j], arglist$ipars["skew",1], arglist$ipars["shape",1],
					arglist$ipars["ghlambda",1])
			tmp = starforc.mc(esim[,j], arglist)
			arglist$y[N+j] = mean(tmp)
			ydist[,j] = tmp
		} else if(arglist$method=="an.kernel"){
			arglist$sigma = arglist$h[N+j]
			arglist$m = N+j
			fit = .kfit(arglist$zresiduals[1:N])
			arglist$fit = fit
			tmp = integrate(starforc.an3, -Inf, Inf, arglist = arglist)
			arglist$y[N+j] = tmp$value
		} else{
			arglist$m = N+j
			arglist$skew = arglist$ipars["skew",1]
			arglist$shape = arglist$ipars["shape",1]
			arglist$ghlambda = arglist$ipars["ghlambda",1]
			arglist$epsilon = NULL
			tmp = integrate(starforc.an1, -Inf, Inf, arglist = arglist)
			arglist$y[N+j] = tmp$value
		}
	}
	arglist$ydist = ydist
	arglist$esim = esim
	return(arglist)
}
starforc.an1 = function(e, arglist)
{
	n = length(e)
	T = arglist$m
	# Analytical forecast
	# use the starxsim C code
	ans = sapply(1:n, function(i){
				# starxsim expects a vector...
				# should supply c(arglist$epsx[1:(T-1)],e[i])
				# resx = as.double(rep(e[i],T))
				resx = c(arglist$epsx[1:(T-1)],e[i])
				tmp = switch(arglist$modelinc[46],
						try(.C("starxsim1", model = as.integer(arglist$modelinc), pars = as.double(arglist$ipars[,1]), 
										idx = as.integer(arglist$idx[,1]-1), x = as.double(arglist$y[1:T]), 
										s = double(arglist$modelinc[46]*T), res = resx, 
										prob = as.double(arglist$probs[1:T,,drop=FALSE]), 
										constm = as.double(arglist$constm[1:T,,drop=FALSE]), m = as.integer(T-1),
										T = as.integer(T), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim2", model = as.integer(arglist$modelinc), pars = as.double(arglist$ipars[,1]), 
										idx = as.integer(arglist$idx[,1]-1), x = as.double(arglist$y[1:T]), 
										s = double(arglist$modelinc[46]*T), res = resx, prob = as.double(arglist$probs[1:T,,drop=FALSE]), 
										constm = as.double(arglist$constm[1:T,,drop=FALSE]), m = as.integer(T-1),
										T = as.integer(T), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim3", model = as.integer(arglist$modelinc), pars = as.double(arglist$ipars[,1]), 
										idx = as.integer(arglist$idx[,1]-1), x = as.double(arglist$y[1:T]), 
										s = double(arglist$modelinc[46]*T), res = resx, prob = as.double(arglist$probs[1:T,,drop=FALSE]), 
										constm = as.double(arglist$constm[1:T,,drop=FALSE]), m = as.integer(T-1),
										T = as.integer(T), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim4", model = as.integer(arglist$modelinc), pars = as.double(arglist$ipars[,1]), 
										idx = as.integer(arglist$idx[,1]-1), x = as.double(arglist$y[1:T]), 
										s = double(arglist$modelinc[46]*T), res = resx, prob = as.double(arglist$probs[1:T,,drop=FALSE]), 
										constm = as.double(arglist$constm[1:T,,drop=FALSE]), m = as.integer(T-1),
										T = as.integer(T), PACKAGE="twinkle"), silent = TRUE)
				)
				return(tmp$x[T])
			})
	ans * ddist(arglist$distribution, e, 0, arglist$sigma, arglist$skew, arglist$shape, arglist$ghlambda)
}

starforc.an2 = function(e, arglist)
{
	n = length(e)	
	T = arglist$m
	# Analytical forecast
	# use the starxsim C code
	ans = sapply(1:n, function(i){
				resx = c(arglist$epsx[1:(T-1)],e[i])
				tmp = switch(arglist$modelinc[46],
						try(.C("starxsim1", model = as.integer(arglist$modelinc), pars = as.double(arglist$ipars[,1]), 
										idx = as.integer(arglist$idx[,1]-1), x = as.double(arglist$y[1:T]), 
										s = double(arglist$modelinc[46]*T), res = resx, prob = as.double(arglist$probs[1:T,,drop=FALSE]), 
										constm = as.double(arglist$constm[1:T,,drop=FALSE]), m = as.integer(T-1),
										T = as.integer(T), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim2",model = as.integer(arglist$modelinc), pars = as.double(arglist$ipars[,1]), 
										idx = as.integer(arglist$idx[,1]-1), x = as.double(arglist$y[1:T]), 
										s = double(arglist$modelinc[46]*T), res = resx, prob = as.double(arglist$probs[1:T,,drop=FALSE]), 
										constm = as.double(arglist$constm[1:T,,drop=FALSE]), m = as.integer(T-1),
										T = as.integer(T), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim3", model = as.integer(arglist$modelinc), pars = as.double(arglist$ipars[,1]), 
										idx = as.integer(arglist$idx[,1]-1), x = as.double(arglist$y[1:T]), 
										s = double(arglist$modelinc[46]*T), res = resx, prob = as.double(arglist$probs[1:T,,drop=FALSE]), 
										constm = as.double(arglist$constm[1:T,,drop=FALSE]), m = as.integer(T-1),
										T = as.integer(T), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim4", model = as.integer(arglist$modelinc), pars = as.double(arglist$ipars[,1]), 
										idx = as.integer(arglist$idx[,1]-1), x = as.double(arglist$y[1:T]), 
										s = double(arglist$modelinc[46]*T), res = resx, prob = as.double(arglist$probs[1:T,,drop=FALSE]), 
										constm = as.double(arglist$constm[1:T,,drop=FALSE]), m = as.integer(T-1),
										T = as.integer(T), PACKAGE="twinkle"), silent = TRUE)
				)
				return(tmp$x[T])
			})
	ans * dkde(e, arglist$fit)
}

starforc.an3 = function(e, arglist)
{
	n = length(e)	
	T = arglist$m
	# Analytical forecast
	# use the starxsim C code
	ans = sapply(1:n, function(i){
				# starxsim expects a vector...
				resx = c(arglist$epsx[1:(T-1)],e[i])
				tmp = switch(arglist$modelinc[46],
						try(.C("starxsim1", model = as.integer(arglist$modelinc), pars = as.double(arglist$ipars[,1]), 
										idx = as.integer(arglist$idx[,1]-1), x = as.double(arglist$y[1:T]), 
										s = double(arglist$modelinc[46]*T), res = resx, prob = as.double(arglist$probs[1:T,,drop=FALSE]), 
										constm = as.double(arglist$constm[1:T,,drop=FALSE]), m = as.integer(T-1),
										T = as.integer(T), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim2",model = as.integer(arglist$modelinc), pars = as.double(arglist$ipars[,1]), 
										idx = as.integer(arglist$idx[,1]-1), x = as.double(arglist$y[1:T]), 
										s = double(arglist$modelinc[46]*T), res = resx, prob = as.double(arglist$probs[1:T,,drop=FALSE]), 
										constm = as.double(arglist$constm[1:T,,drop=FALSE]), m = as.integer(T-1),
										T = as.integer(T), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim3", model = as.integer(arglist$modelinc), pars = as.double(arglist$ipars[,1]), 
										idx = as.integer(arglist$idx[,1]-1), x = as.double(arglist$y[1:T]), 
										s = double(arglist$modelinc[46]*T), res = resx, prob = as.double(arglist$probs[1:T,,drop=FALSE]), 
										constm = as.double(arglist$constm[1:T,,drop=FALSE]), m = as.integer(T-1),
										T = as.integer(T), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim4", model = as.integer(arglist$modelinc), pars = as.double(arglist$ipars[,1]), 
										idx = as.integer(arglist$idx[,1]-1), x = as.double(arglist$y[1:T]), 
										s = double(arglist$modelinc[46]*T), res = resx, prob = as.double(arglist$probs[1:T,,drop=FALSE]), 
										constm = as.double(arglist$constm[1:T,,drop=FALSE]), m = as.integer(T-1),
										T = as.integer(T), PACKAGE="twinkle"), silent = TRUE)
				)
				return(tmp$x[T])
			})
	ans * dkde(e, arglist$fit)*arglist$sigma
}

starforc.mc = function(e, arglist)
{
	n = length(e)	
	T = arglist$m
	# Analytical forecast
	# use the starxsim C code
	ans = sapply(1:n, function(i){
				# starxsim expects a vector...
				resx = c(arglist$epsx[1:(T-1)],e[i])
				tmp = switch(arglist$modelinc[46],
						try(.C("starxsim1", model = as.integer(arglist$modelinc), pars = as.double(arglist$ipars[,1]), 
										idx = as.integer(arglist$idx[,1]-1), x = as.double(arglist$y[1:T]), 
										s = double(arglist$modelinc[46]*T), res = resx, prob = as.double(arglist$probs[1:T,,drop=FALSE]), 
										constm = as.double(arglist$constm[1:T,,drop=FALSE]), m = as.integer(T-1),
										T = as.integer(T), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim2",model = as.integer(arglist$modelinc), pars = as.double(arglist$ipars[,1]), 
										idx = as.integer(arglist$idx[,1]-1), x = as.double(arglist$y[1:T]), 
										s = double(arglist$modelinc[46]*T), res = resx, prob = as.double(arglist$probs[1:T,,drop=FALSE]), 
										constm = as.double(arglist$constm[1:T,,drop=FALSE]), m = as.integer(T-1),
										T = as.integer(T), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim3", model = as.integer(arglist$modelinc), pars = as.double(arglist$ipars[,1]), 
										idx = as.integer(arglist$idx[,1]-1), x = as.double(arglist$y[1:T]), 
										s = double(arglist$modelinc[46]*T), res = resx, prob = as.double(arglist$probs[1:T,,drop=FALSE]), 
										constm = as.double(arglist$constm[1:T,,drop=FALSE]), m = as.integer(T-1),
										T = as.integer(T), PACKAGE="twinkle"), silent = TRUE),
						try(.C("starxsim4", model = as.integer(arglist$modelinc), pars = as.double(arglist$ipars[,1]), 
										idx = as.integer(arglist$idx[,1]-1), x = as.double(arglist$y[1:T]), 
										s = double(arglist$modelinc[46]*T), res = resx, prob = as.double(arglist$probs[1:T,,drop=FALSE]), 
										constm = as.double(arglist$constm[1:T,,drop=FALSE]), m = as.integer(T-1),
										T = as.integer(T), PACKAGE="twinkle"), silent = TRUE)
				)
				return(tmp$x[T])
			})
	ans
}


build.lagmatrix = function(modelinc, s, y, ylags)
{
	if(modelinc[47]==0){
		if(modelinc[49]==2){
			s = as.matrix(s)
			XL = s
		} else{
			y = as.matrix(y)
			n = NROW(y)
			XL = matrix(NA, ncol = length(ylags), nrow = n)
			for(i in 1:length(ylags)){
				XL[(ylags[i]+1):n,i] = y[1:(n-ylags[i]),i]
			}
			XL[is.na(XL)]=0
		}
	} else{
		XL = NULL
	}
	return(XL)
}
.kfit = function(Z){
	H.pi = hpi(Z)
	fit = kde(Z, h=H.pi)
	return(fit)
}

dstar1f = function(arglist)
{
	j = arglist$j
	# unless the user passes his own fixed probabilities (i.e. weights in the 1-state case)
	# we set the value equal to the previous value.
	N = arglist$np
	arglist$probs[N+j,1] = arglist$probs[N+j-1,1]
	return(list(probs = arglist$probs, pmu = arglist$pmu))
}

# CONTINUE HERE
dstar2f = function(arglist)
{
	j = arglist$j
	ipars = arglist$ipars
	probs = arglist$probs
	estidx = arglist$estidx
	idx = arglist$idx
	XL = arglist$XL
	N = arglist$np
	modelinc = arglist$modelinc
	pmu = arglist$pmu
	if(modelinc[21]>0){
		beta = ipars[idx["s1.beta",1],1]
	} else{
		beta = 0
	}
	alpha = ipars[idx["s1.alpha",1]:idx["s1.alpha",2],1]
	cnst = ipars[idx["s1.c",1],1]
	gamma = ipars[idx["s1.gamma",1],1]
	pmu[N+j,1] = gamma*(-cnst + as.numeric(XL[N+j,,drop=FALSE]%*%alpha))
	if(modelinc[21]>0) pmu[N+j,1] = pmu[N+j,1] + beta*pmu[N+j-1,1]
	probs[N+j,1] = 1/(1+exp(-pmu[N+j,1]))
	probs[N+j,2] = 1 - probs[N+j,1]
	return(list(probs = probs, pmu = pmu))
}


dstar3f = function(arglist)
{
	j = arglist$j
	ipars = arglist$ipars
	probs = arglist$probs
	estidx = arglist$estidx
	idx = arglist$idx
	XL = arglist$XL
	N = arglist$np
	modelinc = arglist$modelinc
	pmu = arglist$pmu
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
	pmu[N+j,1] = gamma1*(-cnst1 + as.numeric(XL[N+j,,drop=FALSE]%*%alpha1))
	pmu[N+j,2] = gamma2*(-cnst2 + as.numeric(XL[N+j,,drop=FALSE]%*%alpha2))
	if(modelinc[21]>0) pmu[N+j,1] = pmu[N+j,1] + beta1*pmu[N+j-1,1]
	if(modelinc[25]>0) pmu[N+j,2] = pmu[N+j,2] + beta2*pmu[N+j-1,2]
	probs[N+j,1] = 1/(1+exp(-pmu[N+j,1]))
	probs[N+j,2] = 1 - probs[N+j,1]
	
	p1 = 1/(1+exp(pmu[N+j,1]))
	p2 = 1/(1+exp(pmu[N+j,2]))
	p12 = p1+p2
	p3 = 1/(1+p12)
	p1 = p1/(1+p12)
	p2 = p2/(1+p12)
	probs[N+j,1] = p1
	probs[N+j,2] = p2
	probs[N+j,3] = p3
	return(list(probs = probs, pmu = pmu))
}

dstar4f = function(arglist)
{
	j = arglist$j
	ipars = arglist$ipars
	probs = arglist$probs
	estidx = arglist$estidx
	idx = arglist$idx
	XL = arglist$XL
	N = arglist$np
	modelinc = arglist$modelinc
	pmu = arglist$pmu
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
	
	pmu[N+j,1] = gamma1*(-cnst1 + as.numeric(XL[N+j,,drop=FALSE]%*%alpha1))
	pmu[N+j,2] = gamma2*(-cnst2 + as.numeric(XL[N+j,,drop=FALSE]%*%alpha2))
	pmu[N+j,3] = gamma3*(-cnst3 + as.numeric(XL[N+j,,drop=FALSE]%*%alpha3))
	if(modelinc[21]>0) pmu[N+j,1] = pmu[N+j,1] + beta1*pmu[N+j-1,1]
	if(modelinc[25]>0) pmu[N+j,2] = pmu[N+j,2] + beta2*pmu[N+j-1,2]
	if(modelinc[29]>0) pmu[N+j,3] = pmu[N+j,3] + beta3*pmu[N+j-1,3]	
	p1 = 1/(1+exp(pmu[N+j,1]))
	p2 = 1/(1+exp(pmu[N+j,2]))
	p3 = 1/(1+exp(pmu[N+j,3]))
	p123 = p1+p2+p3
	p4 = 1/(1+p123)
	p1 = p1/(1+p123)
	p2 = p2/(1+p123)
	p3 = p3/(1+p123)
	probs[N+j,1] = p1
	probs[N+j,2] = p2
	probs[N+j,3] = p3
	probs[N+j,4] = p4
	return(list(probs = probs, pmu = pmu))
}
