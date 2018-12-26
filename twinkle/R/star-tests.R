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
.information.test = function(LLH, nObs, nPars)
{
	AIC  = (-2*LLH)/nObs + 2 * nPars/nObs
	BIC  = (-2*LLH)/nObs + nPars * log(nObs)/nObs
	SIC  = (-2*LLH)/nObs + log((nObs+2*nPars)/nObs)
	HQIC = (-2*LLH)/nObs + (2*nPars*log(log(nObs)))/nObs
	informationTests = list(AIC = AIC, BIC = BIC, SIC = SIC, HQIC = HQIC)
	return(informationTests)
}

# McLeod and Li Test of Heteroscedasticity:
# asymptotically chisquared distribution, assuming that 
# lags/n is small and lags is moderately large

mclitest = function(x, lag.max = 10, p = 0, q = 0)
{
	n = length(x)
	y = acf(x^2, lag.max = lag.max, plot  = FALSE)$acf^2
	STATISTIC = n*(n+2)*cumsum( (1/(n-(1:lag.max))) * y )
	pvalues = 1-pchisq(STATISTIC, 1:lag.max - p - q)
	return(list(STATISTIC = STATISTIC, pvalues = pvalues))
}

bgstartest = function(object, method = "Chisq", lag.max = 10)
{
	z = score(object)
	r = as.numeric(residuals(object))
	n = nrow(z)
	e = matrix(0, ncol = lag.max, nrow = n)
	for(i in 1:lag.max){e[(i+1):n,i] = r[1:(n-i)]}
	z = cbind(z, e)
	z = z[-c(1:lag.max),,drop=FALSE]
	r = r[-c(1:lag.max)]
	Ru = residuals(lm(r~z-1))
	nq = ncol(z)
	T = nrow(z)
	if(method=="F"){
		# F-Statistic ~ F(df1=lag.max, df2=T-nq-lag.max) 
		stat = ((T-nq-lag.max)/lag.max)*(sum(r^2) - sum(Ru^2))/(sum(Ru^2))
		dof = c(lag.max, T-nq-lag.max)
		pvalue = pf(stat, df1 = dof[1], df2 = dof[2], lower.tail=FALSE)
		names(dof)=c("df1","df2")
	} else if(method=="Chisq"){
		# LM Statistic ~ chisq(lag.max)
		stat = T*(sum(r^2) - sum(Ru^2))/(sum(Ru^2))
		dof = lag.max
		pvalue = pchisq(stat, dof, lower.tail = FALSE)
		names(dof)="df"
	} else{
		stop("\nunrecognized method.\n")
	}
	names(stat)<-"LM test"
	names(pvalue)<-"p.value"
	test = list(statistic = stat, parameter = dof,
			method = paste("LM Test (",method,") for residual serial correlation in STAR model",sep=""), 
			p.value = pvalue,
			data.name = "STARfit object", coeffficients = NULL)
	class(test)<-"htest"
	return(test)
}


# LM linearity testing against 2 regime STAR
#   Performs an 3rd order Taylor expansion LM test
# In the case that s_t has more than 1 variable, then
# we can only use the fitted object since we need the coefficients
# to calculate the linear combination of s_t \alpha'
nonlinearTest.fit <- function(object, data, robust = FALSE, sig.level = 0.05)
{
	model = object@model
	# we subtract the max AR lags
	p = max(model$modelinc[c(2,6)])
	# check for MA terms (not supported)
	data = matrix(coredata(object@model$modeldata$data[1:object@model$modeldata$T]), ncol=1)
	MA = sum(model$modelinc[c(4,8,12,16,17)])
	if(model$modelinc[46]!=2) stop("\ntwinkle-->error: nonlinearTest only applies to 2-state model\n")
	if(MA>0) stop("\ntwinkle-->error: nonlinearTest does not support MA terms in STAR model\n")
	N = NROW(data) - p
	if(p>0){
		X = coredata(modelmatrix(object, linear=TRUE)[-c(1:p),])
	} else{
		X = coredata(modelmatrix(object, linear=TRUE))
	}
	if(p>0){
		S = coredata(modelmatrix(object, linear=FALSE)[-c(1:p),])
	} else{
		S = coredata(modelmatrix(object, linear=FALSE))
	}
	if(ncol(S)>1){
		warning("\ntwinkle-->warning: more than 1 threshold variable in equation...using state coefficients to collapse to vector.")
		alpha = object@fit$ipars[model$pidx["s1.alpha",1]:model$pidx["s1.alpha",2],1]
		S = S %*% alpha
	}
	if(p>0) Y = data[-c(1:p),,drop=FALSE] else Y = data
	# Regressors under the null
	# We consider the complete set (i.e. max(state1,state2) in terms of AR and Intercept parameters)
	Nn = substr(colnames(X),4,20)
	if(any(duplicated(Nn))) X0 = X[,-which(duplicated(Nn)),drop=FALSE]
	Nm = substr(colnames(X0),4,20)
	colnames(X0) = Nm
	if(any(Nm=="phi0")){
		idxc = which(Nm=="phi0")
		hasintercept = TRUE
	} else{
		hasintercept = FALSE
	}
	# Linear Model (Null)
	b0 = solve(t(X0) %*% X0)%*%t(X0)%*%Y
	e0 = Y - X0 %*% b0
	SSE0 = sum(e0^2)
	
	# Alternative Regressors under the alternative
	if(model$modelinc[49]==1){
		if(hasintercept){
			stmp = repmat(S, 1, ncol(X0)-1)
			X1 = cbind(X0[,-idxc]*stmp, X0[,-idxc]*stmp^2, X0[,-idxc]*stmp^3)
		} else{
			stmp = repmat(S, 1, ncol(X0))
			X1 = cbind(X0*stmp, X0*stmp^2, X0*stmp^3)
		}
	} else{
		stmp = repmat(S, 1, ncol(X0))
		X1 = cbind(X0*stmp, X0*stmp^2, X0*stmp^3)
	}
	Z = cbind(X0, X1)
	# standardize Z
	if(hasintercept){
		xs = apply(Z[,-1], 2, "sd")
		Z[,-1] = scale(Z[,-1], center=FALSE, scale=xs)
	} else{
		xs = apply(Z, 2, "sd")
		Z = scale(Z, center=FALSE, scale=xs)
	}
	n0 = NCOL(X0)
	n1 = NCOL(X1)
	if(!robust){	
		# Alternative regression
		b1   = ginv(t(Z)%*%Z, tol=.Machine$double.eps)%*%t(Z)%*%e0
		e1   = e0 - Z%*%b1
		SSE1 = sum(e1^2)
		
		# Statistics
		F.statistic = ((SSE0-SSE1)/n1)/(SSE1/(N-n0-n1))
		F.pvalue = pf(F.statistic, n1, N - n1 - n0, lower.tail = FALSE)
		F.H0 = qf(sig.level, n1, N-n1-n0, lower.tail = FALSE)>=F.statistic
		chisq.statistic = N*(SSE0 - SSE1)/SSE0
		chisq.pvalue = pchisq(chisq.statistic, n1, lower.tail=FALSE)
		chisq.H0 = chisq.statistic<=qchisq(sig.level,n1, lower.tail = FALSE)
		
		return(list(hypothesis = "H0: Linear Model",
						F.statistic = F.statistic, F.pvalue = F.pvalue, F.dof = c(n1, N-n0-n1), 
						F.decision = ifelse(F.H0, "Fail to Reject", "Reject"),
						chisq.statistic = chisq.statistic, chisq.pvalue = chisq.pvalue, 
						chisq.dof = n1, chisq.decision = ifelse(chisq.H0, "Fail to Reject", "Reject")))
		
	} else{
		b1  = solve(t(X0)%*%X0)%*%t(X0)%*%X1
		R  = X1 - X0%*%b1
		one = matrix(1,N,1)
		Re0  = repmat(e0,1,n1)*R
		b2 = solve(t(Re0)%*%Re0)%*%t(Re0)%*%one
		e1 = one - Re0%*%b2
		SSE1 = sum(e1[,1]^2)
		
		chisq.statistic = N - SSE1
		chisq.pvalue = pchisq(chisq.statistic, n1, lower.tail=FALSE)
		chisq.H0 = chisq.statistic<=qchisq(sig.level, n1, lower.tail = FALSE)
		
		F.statistic = ((N - SSE1)/n1)/(SSE1/(N-n0-n1))
		F.pvalue = pf(F.statistic, n1, N - n1 - n0, lower.tail = FALSE)
		F.H0 = qf(sig.level, n1, N-n1-n0, lower.tail = FALSE)>=F.statistic
		
		return(list(hypothesis = "H0: Linear Model",
						F.statistic = F.statistic, F.pvalue = F.pvalue, F.dof = c(n1, N-n0-n1), 
						F.decision = ifelse(F.H0, "Fail to Reject", "Reject"),
						chisq.statistic = chisq.statistic, chisq.pvalue = chisq.pvalue, 
						chisq.dof = n1, chisq.decision = ifelse(chisq.H0, "Fail to Reject", "Reject")))
	}
}

nonlinearTest.spec <- function(object, data, robust = FALSE, sig.level = 0.05)
{
	if(missing(data) | is.null(data)) stop("\ntwinkle-->error: nonlinearTest with STARspec object requires a valid xts data object\n")
	if(!is.xts(data)) stop("\ntwinkle-->error: nonlinearTest with STARspec object requires a valid xts data object\n")
	model = object@model
	# we subtract the max AR lags
	p = max(model$modelinc[c(2,6)])
	# check for MA terms (not supported)
	MA = sum(model$modelinc[c(4,8,12,16,17)])
	if(model$modelinc[46]!=2) stop("\ntwinkle-->error: nonlinearTest only applies to 2-state model\n")
	if(MA>0) stop("\ntwinkle-->error: nonlinearTest does not support MA terms in STAR model\n")
	N = nrow(data) - p
	if(p>0){
		X = coredata(modelmatrix(object, data = data, linear=TRUE)[-c(1:p),])
	} else{
		X = coredata(modelmatrix(object, data = data, linear=TRUE))
	}
	if(p>0){
		S = coredata(modelmatrix(object, data = data, linear=FALSE)[-c(1:p),])
	} else{
		S = coredata(modelmatrix(object, data = data, linear=FALSE))
	}
	if(ncol(S)>1){
		warning("\ntwinkle-->warning: more than 1 threshold variable in equation...using state coefficients to collapse to vector.")
		alpha = object@model$ipars[model$pidx["s1.alpha",1]:model$pidx["s1.alpha",2],1]
		S = S %*% alpha
	}
	if(p>0) Y = matrix(coredata(data)[-c(1:p)], ncol=1) else Y = matrix(coredata(data), ncol=1)
	# Regressors under the null
	# We consider the complete set (i.e. max(state1,state2) in terms of AR and Intercept parameters)
	Nn = substr(colnames(X),4,20)
	if(any(duplicated(Nn))) X0 = X[,-which(duplicated(Nn)),drop=FALSE]
	Nm = substr(colnames(X0),4,20)
	colnames(X0) = Nm
	if(any(Nm=="phi0")){
		idxc = which(Nm=="phi0")
		hasintercept = TRUE
	} else{
		hasintercept = FALSE
	}
	# Linear Model (Null)
	b0 = solve(t(X0) %*% X0)%*%t(X0)%*%Y
	e0 = Y - X0 %*% b0
	SSE0 = sum(e0^2)
	
	# Alternative Regressors under the alternative
	if(model$modelinc[49]==1){
		if(hasintercept){
			stmp = repmat(S, 1, ncol(X0)-1)
			X1 = cbind(X0[,-idxc]*stmp, X0[,-idxc]*stmp^2, X0[,-idxc]*stmp^3)
		} else{
			stmp = repmat(S, 1, ncol(X0))
			X1 = cbind(X0*stmp, X0*stmp^2, X0*stmp^3)
		}
	} else{
		stmp = repmat(S, 1, ncol(X0))
		X1 = cbind(X0*stmp, X0*stmp^2, X0*stmp^3)
	}
	Z = cbind(X0, X1)
	# standardize Z
	if(hasintercept){
		xs = apply(Z[,-1], 2, "sd")
		Z[,-1] = scale(Z[,-1], center=FALSE, scale=xs)
	} else{
		xs = apply(Z, 2, "sd")
		Z = scale(Z, center=FALSE, scale=xs)
	}
	n0 = NCOL(X0)
	n1 = NCOL(X1)
	if(!robust){
		# Alternative regression
		b1   = ginv(t(Z)%*%Z, tol=.Machine$double.eps)%*%t(Z)%*%e0
		e1   = e0 - Z%*%b1
		SSE1 = sum(e1^2)

		# Statistics
		F.statistic = ((SSE0-SSE1)/n1)/(SSE1/(N-n0-n1))
		F.pvalue = pf(F.statistic, n1, N - n1 - n0, lower.tail = FALSE)
		F.H0 = qf(sig.level, n1, N-n1-n0, lower.tail = FALSE)>=F.statistic
		chisq.statistic = N*(SSE0 - SSE1)/SSE0
		chisq.pvalue = pchisq(chisq.statistic, n1, lower.tail=FALSE)
		chisq.H0 = chisq.statistic<=qchisq(sig.level,n1, lower.tail = FALSE)
		
		return(list(hypothesis = "H0: Linear Model",
						F.statistic = F.statistic, F.pvalue = F.pvalue, F.dof = c(n1, N-n0-n1), 
						F.decision = ifelse(F.H0, "Fail to Reject", "Reject"),
						chisq.statistic = chisq.statistic, chisq.pvalue = chisq.pvalue, 
						chisq.dof = n1, chisq.decision = ifelse(chisq.H0, "Fail to Reject", "Reject")))
	} else{
		b1  = solve(t(X0)%*%X0)%*%t(X0)%*%X1
		R  = X1 - X0%*%b1
		one = matrix(1,N,1)
		Re0  = repmat(e0,1,n1)*R
		b2 = solve(t(Re0)%*%Re0)%*%t(Re0)%*%one
		e1 = one - Re0%*%b2
		SSE1 = sum(e1[,1]^2)
		
		chisq.statistic = N - SSE1
		chisq.pvalue = pchisq(chisq.statistic, n1, lower.tail=FALSE)
		chisq.H0 = chisq.statistic<=qchisq(sig.level,n1, lower.tail = FALSE)
		
		F.statistic = ((N - SSE1)/n1)/(SSE1/(N-n0-n1))
		F.pvalue = pf(F.statistic, n1, N - n1 - n0, lower.tail = FALSE)
		F.H0 = qf(sig.level, n1, N-n1-n0, lower.tail = FALSE)>=F.statistic
		
		return(list(hypothesis = "H0: Linear Model",
						F.statistic = F.statistic, F.pvalue = F.pvalue, F.dof = c(n1, N-n0-n1), 
						F.decision = ifelse(F.H0, "Fail to Reject", "Reject"),
						chisq.statistic = chisq.statistic, chisq.pvalue = chisq.pvalue, 
						chisq.dof = n1, chisq.decision = ifelse(chisq.H0, "Fail to Reject", "Reject")))
	}
}
# GIRF\left(k,\delta,\Omega_{t-1}\right) = E\left[y_{t+k}\left|{\varepsilon_t = \delta,\Omega_{t-1} \right.}\right] - E\left[y_{t+k}\left|{\Omega_{t-1}\right.}\right]
# k=0,1,\ldots,
# shock (delta) = standardized shock
# start (k)
# end (simulation periods after shock)
# sampling (type of simulation)
# ssim = list of matrices (length = end-start) for the "s" variables
star.girf = function(object, shock = 1, start=1, end = start+10, sampling = c("parametric","empirical","kernel","spd"), ssim = NULL)
{
	# choose the shock (assuming independence from other variables)
	# simulate based on shock
	for(i in 1:start){
		
	}
}