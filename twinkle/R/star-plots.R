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

#------------------------------------------------------------------------------------
# Fitted Plots
#------------------------------------------------------------------------------------

plot.starfit = function(x, ...)
{
	modelinc = x@model$modelinc
	def.par <- par(no.readonly = TRUE)
	T = x@model$modeldata$T
	D = as.Date(x@model$modeldata$index[1:T])
	if(modelinc[46]==4){
		mat = matrix(c(1,1,1,2,3,4,5,5,5),3,3,byrow=TRUE)
		nf=layout(mat)
		plot(D, as.numeric(states(x)[,1]), ylim=c(0,1), main = "Prob[state={1,2,3,4}]", type="l", ylab="Probability", xlab="Time")
		lines(D, as.numeric(states(x)[,2]), col=2, lty=2)
		lines(D, as.numeric(states(x)[,3]), col=3, lty=2)
		lines(D, as.numeric(states(x)[,4]), col=4, lty=2)
		grid()
		legend("topleft", c("[state=1]", "[state=2]", "[state=3]", "[state=4]"), col=1:4, lty=c(1,2,2,2), bty="n")
		if(modelinc[47]==0){
			plot(D, as.numeric(states(x, type="pmu")[,1]), main = "State Dynamics[1]", type="l", ylab=expression(pi['1,t']), xlab="Time")
			grid()
			plot(D, as.numeric(states(x, type="pmu")[,2]), main = "State Dynamics[2]", type="l", ylab=expression(pi['2,t']), xlab="Time")
			grid()
			plot(D, as.numeric(states(x, type="pmu")[,3]), main = "State Dynamics[3]", type="l", ylab=expression(pi['3,t']), xlab="Time")
			grid()
		}
		plot(D, as.numeric(x@model$modeldata$data[1:x@model$modeldata$T]), main = "Fitted vs Actual", type="l", ylab="Value", xlab="Time")
		lines(D, as.numeric(fitted(x)), col="steelblue", lwd=2)
		grid()
	} else if(modelinc[46]==3){
		mat = matrix(c(1,1,2,3,4,4),3,2,byrow=TRUE)
		nf=layout(mat)
		plot(D, as.numeric(states(x)[,1]), ylim=c(0,1), main = "Prob[state={1,2,3}]", type="l", ylab="Probability", xlab="Time")
		lines(D, as.numeric(states(x)[,2]), col=2, lty=2)
		lines(D, as.numeric(states(x)[,3]), col=3, lty=2)
		grid()
		legend("topleft", c("[state=1]", "[state=2]", "[state=3]"), col=1:3, lty=c(1,2,2), bty="n")
		if(modelinc[47]==0){
			plot(D, as.numeric(states(x, type="pmu")[,1]), main = "State Dynamics[1]", type="l", ylab=expression(pi['1,t']), xlab="Time")
			grid()
			plot(D, as.numeric(states(x, type="pmu")[,2]), main = "State Dynamics[2]", type="l", ylab=expression(pi['2,t']), xlab="Time")
			grid()
		}
		plot(D, as.numeric(x@model$modeldata$data[1:x@model$modeldata$T]), main = "Fitted vs Actual", type="l", ylab="Value", xlab="Time")
		lines(D, as.numeric(fitted(x)), col="steelblue", lwd=2)
		grid()
	} else if(modelinc[46]==2){
		par(mfrow=c(3,1))
		plot(D, as.numeric(states(x)[,1]), main = "Prob[state=1]", type="l", ylab="Probability", xlab="Time")
		grid()
		if(modelinc[47]==0){
			plot(D, as.numeric(states(x, type="pmu")), main = "State Dynamics", type="l", ylab=expression(pi[t]), xlab="Time")
			grid()
		}
		plot(D, as.numeric(x@model$modeldata$data[1:x@model$modeldata$T]), main = "Fitted vs Actual", type="l", ylab="Value", xlab="Time")
		lines(D, as.numeric(fitted(x)), col="steelblue", lwd=2)
		grid()
	} else{
		plot(D, as.numeric(x@model$modeldata$data[1:x@model$modeldata$T]), main = "Fitted vs Actual", type="l", ylab="Value", xlab="Time")
		lines(D, as.numeric(fitted(x)), col="steelblue", lwd=2)
	}
	par(def.par)
	return(invisible())
}


####################################################################################
# Transition Function for z
.transfunall = function(object)
{
	s = as.numeric(states(object, "u"))
	p = as.numeric(states(object, "prob")[,1])
	plot(s, p, type="l")
}


# only for 2 states
trans2fun2d = function(object, colidx = 1, fixed.values = NULL, doplot = TRUE, ...)
{
	if(object@model$modelinc[46]!=2) stop("\ntwinkle-->error: function only works for the 2 states case at present.")
	model = object@model
	modelinc = model$modelinc
	T = model$modeldata$T
	if(modelinc[47]==0){
		if(modelinc[49]==2){
			XL = model$modeldata$s[1:T, , drop = FALSE]
		} else{
			if(modelinc[48]==1) ytmp = model$modeldata$fun(as.numeric(model$modeldata$data[1:T])) else ytmp = as.numeric(model$modeldata$data[1:T])
			XL = NULL
			for(i in 1:length(model$modeldata$ylags)){
				if(i==1) XL = lagf.numeric(ytmp, model$modeldata$ylags[i]) else XL = cbind(XL, lagf.numeric(ytmp, model$modeldata$ylags[i]))
			}
			XL = as.matrix(XL)
			XL[is.na(XL)]=0
		}
		XL = coredata(XL)
	} else{
		stop("\ntwinkle-->error: cannot evaluate transition function since fixed.probs were used\n")
	}
	if(ncol(XL)>1){
		if(is.null(fixed.values)){
			meanX = apply(XL, 2, "median")
			newX  = matrix(meanX, ncol = ncol(XL), nrow = nrow(XL), byrow=TRUE)
			newX[,colidx] = sort(XL[,colidx])
		} else{
			meanX = rep(NA, ncol(XL))
			meanX[-colidx] = fixed.values
			newX  = matrix(meanX, ncol = ncol(XL), nrow = nrow(XL), byrow=TRUE)
			newX[,colidx] = sort(XL[,colidx])
		}
		initp = as.numeric(tail(object@fit$pmu, 1))
	} else{
		newX = matrix(sort(XL[,1]), ncol=1)
		initp = as.numeric(tail(object@fit$pmu, 1))
	}
	sol = dstar2trans(object, newX, initp)
	if(doplot){
		if(modelinc[49]==1){
			px = paste("y['",colidx,",t-",model$modeldata$ylags[colidx],"']",sep="")
			expr2 =  parse(text=paste("G(",paste("y['",colidx,",t-",model$modeldata$ylags[colidx],"']",sep=""),", ...)",sep=""))
		} else{
			px = paste("s['",colidx,",t-d']",sep="")
			expr2 =  parse(text=paste("G(",paste("s['",colidx,",t-d']",sep=""),", ...)",sep=""))
		}
		expr1 = parse(text=px)
		plot(newX[,colidx], sol, type="l", xlab = expr1, ylab =expr2,main="Transition Function", lwd=2, ...)
		grid()
	}
	return(invisible(list(x = newX, y = sol)))
}

trans2fun3d = function(object, colidx = c(1,2), fixed.values = NULL, doplot = TRUE, ...)
{
	if(object@model$modelinc[46]!=2) stop("\ntwinkle-->error: function only works for the 2 states case at present.")
	model = object@model
	modelinc = model$modelinc
	T = model$modeldata$T
	if(modelinc[47]==0){
		if(modelinc[49]==2){
			XL = model$modeldata$s[1:T, , drop = FALSE]
		} else{
			XL = NULL
			if(length(model$modeldata$lags)<2) stop("\ntwinkle-->error: transfun3d requires transition variable to have at least 2 columns.")
			if(modelinc[48]==1) ytmp = model$modeldata$fun(as.numeric(model$modeldata$data[1:T])) else ytmp = as.numeric(model$modeldata$data[1:T])
			for(i in 1:length(model$modeldata$lags)){
				if(i==1) XL = lagf.numeric(ytmp, model$modeldata$ylags[i]) else XL = cbind(XL, lagf.numeric(ytmp, model$modeldata$ylags[i]))
			}
			XL = as.matrix(XL)
			XL[is.na(XL)]=0
		}
		XL = coredata(XL)
	} else{
		stop("\ntwinkle-->error: cannot evaluate transition function since fixed.probs were used\n")
	}
	if(ncol(XL)>2){
		if(is.null(fixed.values)){
			meanX = apply(XL, 2, "mean")
			augX  = matrix(meanX, ncol = ncol(XL), nrow = 1, byrow=TRUE)
			x1 = seq(min(XL[,colidx[1]]), max(XL[,colidx[1]]), length.out=100)
			x2 = seq(min(XL[,colidx[2]]), max(XL[,colidx[2]]), length.out=100)
			d2X   = cbind(x1, x2)
		} else{
			meanX = rep(NA, ncol(XL))
			meanX[-colidx] = fixed.values
			augX  = matrix(meanX, ncol = ncol(XL), nrow = 1, byrow=TRUE)
			x1 = seq(min(XL[,colidx[1]]), max(XL[,colidx[1]]), length.out=100)
			x2 = seq(min(XL[,colidx[2]]), max(XL[,colidx[2]]), length.out=100)
			d2X   = cbind(x1, x2)
		}
		initp = as.numeric(tail(object@fit$pmu, 1))
	} else{
		augX = NULL
		x1 = seq(min(XL[,colidx[1]]), max(XL[,colidx[1]]), length.out=100)
		x2 = seq(min(XL[,colidx[2]]), max(XL[,colidx[2]]), length.out=100)
		d2X   = cbind(x1, x2)
		initp = as.numeric(tail(object@fit$pmu, 1))
	}
	xp = c(model$pars[model$pidx["s1.c",1],1], initp, 
			model$pars[model$pidx["s1.alpha",1]:model$pidx["s1.alpha",2],1][colidx],
			model$pars[model$pidx["s1.beta",1],1],model$pars[model$pidx["s1.gamma",1],1])
	# st[2] is varing
	#1/(1+exp(-(xp[1]+xp[3]*d2X[,1] + xp[4]*d2X[2,2] + xp[5]*xp[2])))
	sol<-try(.C("star3dfun", n = as.integer(100), pars = as.double(xp), zmat=double(100*100),
			v = as.double(d2X), PACKAGE="twinkle"), silent=TRUE)	
	zmat = matrix(sol$zmat,100,100,byrow=TRUE)
	if(doplot){
		if(modelinc[49]==1){
			px1 = paste("y['",colidx[1],",t-",model$modeldata$ylags[colidx[1]],"']",sep="")
			px2 = paste("y['",colidx[2],",t-",model$modeldata$ylags[colidx[2]],"']",sep="")
			expr3 =  parse(text=paste("G(",paste("y[t]",sep=""),", ...)",sep=""))			
		} else{
			px1 = paste("s['",colidx[1],",t-d']",sep="")
			px2 = paste("s['",colidx[2],",t-d']",sep="")
			expr3 =  parse(text=paste("G(",paste("s[t]",sep=""),", ...)",sep=""))
		}
		expr1 = parse(text=px1)
		expr2 = parse(text=px2)
		obj = wireframe(x = zmat, xlab = expr2, ylab =expr1, zlab = expr3, main="Transition 3D Function",
				row.values = d2X[,1], column.values = d2X[,2], drape = TRUE, colorkey = TRUE,
				scales = list(arrows = FALSE))
		print(obj)
		
	}
	return(invisible(list(z = zmat, x = d2X[,1], y = d2X[,2])))
}

trans3fun2d = function(object, colidx = 1, fixed.values = NULL, doplot = TRUE, ...)
{
	if(object@model$modelinc[46]!=3) stop("\ntwinkle-->error: function only works for the 3 state case.")
	model = object@model
	modelinc = model$modelinc
	T = model$modeldata$T
	if(modelinc[47]==0){
		if(modelinc[49]==2){
			XL = model$modeldata$s[1:T, , drop = FALSE]
		} else{
			if(modelinc[48]==1) ytmp = model$modeldata$fun(as.numeric(model$modeldata$data[1:T])) else ytmp = as.numeric(model$modeldata$data[1:T])
			XL = NULL
			for(i in 1:length(model$modeldata$lags)){
				if(i==1) XL = lagf.numeric(ytmp, model$modeldata$ylags[i]) else XL = cbind(XL, lagf.numeric(ytmp, model$modeldata$ylags[i]))
			}
			XL = as.matrix(XL)
			XL[is.na(XL)]=0
		}
		XL = coredata(XL)
	} else{
		stop("\ntwinkle-->error: cannot evaluate transition function since fixed.probs were used\n")
	}
	if(ncol(XL)>1){
		if(is.null(fixed.values)){
			meanX = apply(XL, 2, "median")
			newX  = matrix(meanX, ncol = ncol(XL), nrow = nrow(XL), byrow=TRUE)
			newX[,colidx] = sort(XL[,colidx])
		} else{
			meanX = rep(NA, ncol(XL))
			meanX[-colidx] = fixed.values
			newX  = matrix(meanX, ncol = ncol(XL), nrow = nrow(XL), byrow=TRUE)
			newX[,colidx] = sort(XL[,colidx])
		}
		initp = as.numeric(tail(object@fit$pmu, 1))
	} else{
		newX = matrix(sort(XL[,1]), ncol=1)
		initp = as.numeric(tail(object@fit$pmu, 1))
	}
	sol = dstar3trans2d(object, newX, initp)
	if(doplot){
		par(mfrow=c(3,1))
		if(modelinc[49]==1){
			px1 = paste("y['",1,",",colidx,",t-",model$modeldata$ylags[colidx],"']",sep="")
			px2 = paste("y['",2,",",colidx,",t-",model$modeldata$ylags[colidx],"']",sep="")
			px3 = paste("y['",3,",",colidx,",t-",model$modeldata$ylags[colidx],"']",sep="")
			expr4 =  parse(text=paste("G(",paste("y['",1,",",colidx,",t-",model$modeldata$ylags[colidx],"']",sep=""),", ...)",sep=""))
			expr5 =  parse(text=paste("G(",paste("y['",2,",",colidx,",t-",model$modeldata$ylags[colidx],"']",sep=""),", ...)",sep=""))
			expr6 =  parse(text=paste("G(",paste("y['",3,",",colidx,",t-",model$modeldata$ylags[colidx],"']",sep=""),", ...)",sep=""))
		} else{
			px1 = paste("s['",1,",",colidx,",t-d']",sep="")
			px2 = paste("s['",2,",",colidx,",t-d']",sep="")
			px3 = paste("s['",3,",",colidx,",t-d']",sep="")
			expr4 =  parse(text=paste("G(",paste("s['",1,",",colidx,",t-d']",sep=""),", ...)",sep=""))
			expr5 =  parse(text=paste("G(",paste("s['",2,",",colidx,",t-d']",sep=""),", ...)",sep=""))
			expr6 =  parse(text=paste("G(",paste("s['",3,",",colidx,",t-d']",sep=""),", ...)",sep=""))
		}
		expr1 = parse(text=px1)
		expr2 = parse(text=px2)
		expr3 = parse(text=px3)
		plot(newX[,colidx], sol[,1], type="l", xlab = expr1, ylab =expr4,main="Transition Function (State 1)", lwd=2, ...)
		grid()
		plot(newX[,colidx], sol[,2], type="l", xlab = expr2, ylab =expr5,main="Transition Function (State 2)", lwd=2, ...)
		grid()
		plot(newX[,colidx], sol[,3], type="l", xlab = expr3, ylab =expr6,main="Transition Function (State 3)", lwd=2, ...)
		grid()
	}
	return(invisible(list(x = newX, y = sol)))
}


dstar2trans = function(object, XL, initp)
{
	model = object@model
	idx = model$pidx
	modelinc = model$modelinc
	ipars = object@fit$ipars
	if(modelinc[21]>0){
		beta = ipars[idx["s1.beta",1],1]
	} else{
		beta = 0
	}
	alpha = ipars[idx["s1.alpha",1]:idx["s1.alpha",2],1]
	cnst = ipars[idx["s1.c",1],1]
	gamma = ipars[idx["s1.gamma",1],1]
	
	pmu = apply(XL, 1, function(x){ gamma*(as.numeric(sum(x*alpha))-cnst) })
	if(modelinc[21]>0) pmu = sapply(pmu, function(x) initp*beta+x)
	ans = 1/(1+exp(-pmu))
	return(ans)
}

dstar3trans2d = function(object, XL, initp)
{
	model = object@model
	idx = model$pidx
	modelinc = model$modelinc
	ipars = object@fit$ipars
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
	
	
	pmu1 = apply(XL, 1, function(x){ gamma1*(-cnst1 + as.numeric(sum(x*alpha1))) })
	pmu2 = apply(XL, 1, function(x){ gamma2*(-cnst2 + as.numeric(sum(x*alpha2))) })
	
	if(modelinc[21]>0) pmu1 = sapply(pmu1, function(x) initp[1]*beta1+x)
	if(modelinc[25]>0) pmu2 = sapply(pmu2, function(x) initp[2]*beta2+x)
	
	probs = matrix(0, ncol=3, nrow=nrow(pmu1))
	expmu1 = pmin(exp(100), exp(pmu1))
	expmu2 = pmin(exp(100), exp(pmu2))
	probs[,1] = expmu1/(1+expmu1+expmu2)
	probs[,2] = expmu2/(1+expmu1+expmu2)
	probs[,3] = 1-probs[,1]-probs[,2]	
	return(probs)
}

dstar2trans3d = function(object, d2X, augX, initp, colidx)
{
	model = object@model
	idx = model$pidx
	modelinc = model$modelinc
	ipars = object@fit$ipars
	if(modelinc[21]>0){
		beta = ipars[idx["s1.beta",1],1]
	} else{
		beta = 0
	}
	n = nrow(d2X)
	d3Y = matrix(NA, n, n)
	alpha = ipars[idx["s1.alpha",1]:idx["s1.alpha",2],1]
	cnst = ipars[idx["s1.c",1],1]
	gamma = ipars[idx["s1.gamma",1],1]
	
	# move this to C++
	for(i in 1:n){
		for(j in 1:n){
			if(!is.null(augX)){
				tmp = augX[1,]
				tmp[colidx] = c(d2X[i,1], d2X[j,2])
				d3Y[i,j] = gamma*(-cnst + as.numeric(sum(tmp*alpha)))
				if(modelinc[21]>0) d3Y[i,j] = initp*beta+d3Y[i,j]
			} else{
				tmp = c(d2X[i,1], d2X[j,2])
				d3Y[i,j] = gamma*(-cnst + as.numeric(sum(tmp*alpha)))
				if(modelinc[21]>0) d3Y[i,j] = initp*beta+d3Y[i,j]
			}
		}
	}
	ans = 1/(1+exp(-d3Y))
	return(ans)
}
#------------------------------------------------------------------------------------
# Forecast Plots
# custom functions for plotting forecasts allow either 1-ahead all rolling forecasts else
# given a roll number the n-ahead forecasts
#------------------------------------------------------------------------------------
plot.starforecast = function(x, roll = 0, ...)
{
	if(x@forecast$n.ahead==1) roll = "all"
	if(is.character(roll) && roll=="all"){
		plot.1starforc(x, ...)
	} else{
		plot.nstarforc(x, roll = roll, ...)
	}
}
plot.nstarforc = function(x, roll = 0, ...)
{
	T = x@model$modeldata$T+roll
	n = x@forecast$n.ahead
	origdata = x@model$modeldata$data[1:T]
	forc = fitted(x)[,roll+1]
	fitx = origdata - x@model$modeldata$residuals[1:T]
	usen = min(3*n, T)
	D = c(as.POSIXct(x@model$modeldata$index[1:T]), generatefwd(as.POSIXct(x@model$modeldata$index[T]), length.out = n, by = x@model$modeldata$period))
	D = as.Date(tail(D, n+usen))
	if(substr(x@forecast$method,1,2)=="mc"){
		incd = TRUE
		d = x@forecast$yDist[[roll+1]]
		# no uncertainty in 1-ahead (parameter uncertainty not accounted for by ML method)
		q5 = c(forc[1], apply(d[,-1], 2, function(x) quantile(x, 0.05)))
		q25 = c(forc[1], apply(d[,-1], 2, function(x) quantile(x, 0.25)))
		q75 = c(forc[1], apply(d[,-1], 2, function(x) quantile(x, 0.75)))
		q95 = c(forc[1], apply(d[,-1], 2, function(x) quantile(x, 0.95)))
		ylim = c(min(c(q5, tail(origdata, usen), forc)), max(c(q95,tail(origdata, usen), forc)))
	} else{
		incd = FALSE
		ylim = c(min(c(tail(origdata, usen), forc)), max(c(tail(origdata, usen), forc)))
		
	}
	# Need to also print T:
	plot(D, c(tail(origdata, usen), forc), type="l", col="black", main=paste("STAR Model Forecast [",n,"-ahead]\n","T0:",as.character(x@model$modeldata$index[T]),sep=""), xlab = "Time", ylab="Value", ylim = ylim, ...)
	lines(D, c(tail(fitx, usen), forc), type="l", col="steelblue")
	if(incd){
		polygon(x=c(D[(usen+1):(usen+n)],rev(D[(usen+2):(usen+n)])), 
				y = c(q95, rev(q5[-1])), border=NA, angle = 45, col = "aliceblue")
		lines(D, c(rep(NA, usen), q5), col = "tomato1", lty=2)
		lines(D, c(rep(NA, usen), q95), col = "tomato1", lty=2)
		lines(D, c(rep(NA, usen), q25), col = "tomato1", lty=2)
		lines(D, c(rep(NA, usen), q75), col = "tomato1", lty=2)
		lines(D, c(rep(NA, usen), forc), col = "orange", lwd=2)
		
	} else{
		polygon(x=c(tail(D,n), rev(tail(D,n))),  y = c(rep(ylim[2],n), rep(ylim[1], n)), border=NA, angle = 45, col = "aliceblue")
		lines(D, c(rep(NA, usen), forc), col = "orange", lwd=2)
	}
	grid()
	abline(v=D[usen+1], col="grey", lty=2)
	if(incd){
		legend("topleft", c("Actual", "Estimated","[5%,25%,75%,95%] MC Forecast Quantiles",paste("Forecast [",x@forecast$method,"]",sep="")), 
				col=c("black","steelblue","tomato1","orange"), bty="n",
				lty=c(1,1,2,1), lwd=c(1,1,1,2),cex=0.8)
	} else{
		legend("topleft", c("Actual", "Estimated",paste("Forecast [",x@forecast$method,"]",sep="")), 
				col=c("black","steelblue","orange"), bty="n",
				lty=c(1,1,1), lwd=c(1,1,2),cex=0.8)
	}
	return(invisible(x))
}

plot.1starforc = function(x, roll = 0, ...)
{
	T = x@model$modeldata$T
	nr = x@forecast$n.roll+1
	N = x@forecast$N
	origdata = x@model$modeldata$data[1:T]
	forc = fitted(x)[1,]
	D = as.POSIXct(names(forc))
	fitx = origdata - x@model$modeldata$residuals[1:T]
	usen = min(3*nr, T)
	D = c(as.POSIXct(x@model$modeldata$index[1:T]), D)
	D = as.Date(tail(D, nr+usen))
	ylim = c(min(c(tail(origdata, usen), forc), na.rm=TRUE), max(c(tail(origdata, usen), forc), na.rm=TRUE))
	plot(D, c(tail(origdata, usen), as.numeric(forc)), type="l", col="black", main="STAR Model Rolling Forecast", xlab = "Time", ylab="Value", ylim = ylim, ...)
	polygon(x=c(tail(D,nr), rev(tail(D,nr))),  y = c(rep(ylim[2],nr), rep(ylim[1], nr)), border=NA, angle = 45, col = "aliceblue")
	lines(D, c(rep(NA, usen), forc), col = "orange", lwd=2)
	lines(D, c(tail(fitx, usen), rep(NA, nr)), col = "steelblue")
	grid()
	abline(v=D[usen+1], col="grey", lty=2)
	legend("topleft", c("Actual", "Estimated", "Forecast"), col=c("black","steelblue","tomato1"), bty="n",
			lty=c(1,1,1), lwd = c(1,1,2))
	return(invisible(x))
}


# Plot Probabilities against returns (heatmap)