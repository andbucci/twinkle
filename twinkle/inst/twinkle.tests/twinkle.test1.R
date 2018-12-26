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

twinkle.test1a = function(cluster = NULL){
	tic = Sys.time()
	# Table 3.7 Parameter estimates for a STAR model for weekly
	# returns on the Dutch guilder exchange rate (with additional
	# results from their code output)
	require(quantmod)
	data(forex)
	# to replicate the results of van Dijk and Franses use:
	# next observation carried backward (i.e. fromLast=TRUE)
	# for missing values
	fx = na.locf(forex, fromLast = TRUE)
	fxw = fx[which(weekdays(index(forex))=="Wednesday"),4]
	dx = ROC(fxw, na.pad=FALSE)*100
	# threshold variable a running mean of the absolute returns
	dxv = runMean(abs(dx), n=4)
	dxv[1:3] = c(as.numeric(abs(dx[1,1])), as.numeric(mean(abs(dx[1:2,1]))), 
			as.numeric(mean(abs(dx[1:3,1]))))
	dxv = lag(dxv, 1)
	dxv[1] = 0
	# period considered
	idx = 1:521
	spec1 = starspec(
			mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(0,2), 
					statevar = c("y","s")[2], s = dxv[idx]))
	solver.control=list(maxit=150000, alpha=1, beta=0.4, gamma=1.2, reltol=1e-12, trace=1,method="BFGS",n.restarts=3)
	fit1 = starfit(spec1, data = dx[idx], out.sample = 0, solver = "msoptim", solver.control = solver.control)
	#######################################################################	
	vdfpars = c(-0.18019807, 0.05980984, 0.28705073, 0.21335697, 4.3161180, 1.3554227)
	vdfse   = c( 0.13140430, 0.10138187, 0.10719107, 0.10026741, 1.0771948, 0.1521352)
	names(vdfpars) = names(vdfse) = c("s1.phi0", "s2.phi0", "s2.phi1","s2.phi2", "gamma", "c")
	twinklepars = coef(fit1)
	# equivalent representation
	# twinklepars["s1.c"] = twinklepars["s1.c"]/twinklepars["s1.alpha1"]
	twinklese = sqrt(diag(vcov(fit1, robust=FALSE)))
	
	twinkleSSE = sum(residuals(fit1)^2)
	vdfSSE = 1233.321352

	options(width=150)
	zz <- file("test1a-1.txt", open="wt")
	sink(zz)
	print(data.frame(vdfpars = vdfpars, twinklepars = twinklepars[1:6]))
	print(data.frame(vdfse = vdfse, twinklese = twinklese[1:6]))
	
	print(c("twinkle.SSE"=twinkleSSE, "vdf.SSE" = vdfSSE))
	sink(type="message")
	sink()
	close(zz)
	
	spec2 = starspec(
			mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(1,1), 
					statevar = c("y","s")[2], s = dxv[idx]))
	solver.control=list(maxit=150000, alpha=1, beta=0.4, gamma=1.2, reltol=1e-12, trace=1,method="BFGS")
	fit2 = starfit(spec2, data = dx[idx], out.sample = 0, solver = "strategy", solver.control = solver.control, n=5)
	
	spec3 = starspec(
			mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(1,1), 
					maOrder=c(1,1), matype="state", statevar = c("y","s")[2], s = dxv[idx]))
	solver.control=list(maxit=150000, alpha=1, beta=0.4, gamma=1.2, reltol=1e-12, trace=1,method="BFGS")
	fit3 = starfit(spec3, data = dx[idx], out.sample = 0, solver = "strategy", solver.control = solver.control, n=5)
	
	spec4 = starspec(
			mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(0,0), 
					maOrder=c(1,1), matype="state", statevar = c("y","s")[2], s = dxv[idx]))
	solver.control=list(maxit=150000, alpha=1, beta=0.4, gamma=1.2, reltol=1e-12, trace=1,method="BFGS")
	fit4 = starfit(spec4, data = dx[idx], out.sample = 0, solver = "strategy", solver.control = solver.control, n=5)
	
	spec5 = starspec(
			mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(1,1), 
					maOrder=1, matype="linear", statevar = c("y","s")[2], s = dxv[idx]))
	solver.control=list(maxit=150000, alpha=1, beta=0.4, gamma=1.2, reltol=1e-12, trace=1,method="BFGS")
	fit5 = starfit(spec5, data = dx[idx], out.sample = 0, solver = "strategy", solver.control = solver.control, n=5)
	
	models = c("2-STAR(0,2)", "2-STAR(1,1)", "2-STARMA([1,1],[1,1])", "2-STMA(2,2)", "2-STARLMA([1,1],[1]")
	logl = c(likelihood(fit1), likelihood(fit2), likelihood(fit3), likelihood(fit4), likelihood(fit5))
	AIC = c(infocriteria(fit1)[1], infocriteria(fit2)[1], infocriteria(fit3)[1], infocriteria(fit4)[1], infocriteria(fit5)[1])
	BIC = c(infocriteria(fit1)[2], infocriteria(fit2)[2], infocriteria(fit3)[2], infocriteria(fit4)[2], infocriteria(fit5)[2])
	cfm = matrix(NA, ncol = 5, nrow = 11)
	colnames(cfm) = models
	rownames(cfm) = c("s1.phi0", "s1.phi1", "s1.psi1", "s2.phi0", "s2.phi1", "s2.phi2", "s2.psi1", "psi1", "s1.gamma","s1.c", "sigma")
	cfm[c(1,4,5,6,9:11),1] = coef(fit1) 
	cfm[c(1,2,4,5,9:11), 2] = coef(fit2)
	cfm[c(1,2,3,4,5,7,9:11), 3] = coef(fit3)
	cfm[c(1,3,4,7,9:11), 4] = coef(fit4)
	cfm[c(1,2,4,5,8:11), 5] = coef(fit5)
	
	zz <- file("test1a-2.txt", open="wt")
	sink(zz)
	print(round(rbind(cfm, logl, AIC, BIC), 4))
	sink(type="message")
	sink()
	close(zz)
	
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

# continue with models for variance
twinkle.test1b = function(cluster = NULL){
	require(quantmod)
	data(forex)
	# to replicate the results of van Dijk and Franses use:
	# next observation carried backward (i.e. fromLast=TRUE)
	# for missing values
	fx = na.locf(forex, fromLast = TRUE)
	fxw = fx[which(weekdays(index(forex))=="Wednesday"),4]
	dx = ROC(fxw, na.pad=FALSE)*100
	# threshold variable a running mean of the absolute returns
	dxv = runMean(abs(dx), n=4)
	dxv[1:3] = c(as.numeric(abs(dx[1,1])), as.numeric(mean(abs(dx[1:2,1]))), 
			as.numeric(mean(abs(dx[1:3,1]))))
	dxv = lag(dxv, 1)
	dxv[1] = 0
	# period considered
	idx = 1:521
	spec1 = starspec(
			mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(0,2), 
					statevar = c("y","s")[2], s = dxv[idx]), variance.model=list(dynamic=TRUE,
					model="mixture"))
	solver.control=list(maxit=150000, alpha=1, beta=0.4, gamma=1.2, reltol=1e-12, trace=1,method="BFGS")
	fit1 = starfit(spec1, data = dx[idx], out.sample = 0, solver = "strategy", solver.control = solver.control, n=4)
	
	spec2 = starspec(
			mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(1,1), 
					statevar = c("y","s")[2], s = dxv[idx]), variance.model=list(dynamic=TRUE,
					model="eGARCH"))
	solver.control=list(maxit=1500, alpha=1, beta=0.4, gamma=1.2, reltol=1e-12, trace=1,method="BFGS")
	fit2 = starfit(spec2, data = dx[idx], out.sample = 0, solver = "strategy", solver.control = solver.control, n=4)
	
	spec3 = starspec(
			mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(1,1), 
					maOrder=c(1,1), matype="state", statevar = c("y","s")[2], s = dxv[idx]), 
			variance.model=list(dynamic=TRUE, model="mixture"))
	solver.control=list(maxit=1500, alpha=1, beta=0.4, gamma=1.2, reltol=1e-12, trace=1,method="BFGS")
	fit3 = starfit(spec3, data = dx[idx], out.sample = 0, solver = "strategy", solver.control = solver.control, n=4)
	
	
	spec4 = starspec(
			mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(1,1), 
					maOrder=c(1,1), matype="state", statevar = c("y","s")[2], s = dxv[idx]), 
			variance.model=list(dynamic=TRUE, model="eGARCH"))
	solver.control=list(maxit=1500, alpha=1, beta=0.4, gamma=1.2, reltol=1e-12, trace=1,method="BFGS")
	fit4 = starfit(spec4, data = dx[idx], out.sample = 0, solver = "strategy", solver.control = solver.control, n=4)
	
	models = c("2-STAR(0,2)-NorMix", "2-STAR(1,1)-eGARCH", "2-STARMA([1,1],[1,1])-NorMix", "2-STARMA([1,1],[1,1])-eGARCH")
	logl = c(likelihood(fit1), likelihood(fit2), likelihood(fit3), likelihood(fit4))
	AIC = c(infocriteria(fit1)[1], infocriteria(fit2)[1], infocriteria(fit3)[1], infocriteria(fit4)[1])
	BIC = c(infocriteria(fit1)[2], infocriteria(fit2)[2], infocriteria(fit3)[2], infocriteria(fit4)[2])
	cfm = matrix(NA, ncol = 4, nrow = 16)
	colnames(cfm) = models
	rownames(cfm) = c("s1.phi0", "s1.phi1", "s1.psi1", "s2.phi0", "s2.phi1", "s2.phi2", "s2.psi1", "s1.gamma","s1.c", "sigma","s1.sigma","s2.sigma","omega","alpha1","beta1","gamma1")
	cfm[c(1,4,5,6,8,9,11,12),1] = coef(fit1) 
	cfm[c(1,2,4,5,8,9,13:16), 2] = coef(fit2)
	cfm[c(1,2,3,4,5,7,8,9,11,12), 3] = coef(fit3)
	cfm[c(1,2,3,4,5,7,8,9,13:16), 4] = coef(fit4)
	
	options(width=150)
	zz <- file("test1b-1.txt", open="wt")
	sink(zz)
	print(round(rbind(cfm, logl, AIC, BIC), 4))
	sink(type="message")
	sink()
	close(zz)
	
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

twinkle.test1c = function(cluster = NULL){
	tic = Sys.time()
	data(forex)
	require(quantmod)
	# equivalent ways of using type="y".
	# This method (fun on y) is the best since fun can then be used in the simulation (path depedence)
	fx = na.locf(forex, fromLast = TRUE)
	fxw = fx[which(weekdays(index(forex))=="Wednesday"),4]
	dx = ROC(fxw, na.pad=FALSE)*100
	
	fun = function(x){
		x = as.numeric(x)
		y = runMean(abs(x), n=4)
		y[1:3] = c(abs(x[1]), mean(abs(x[1:2])), mean(abs(x[1:3])))
		return(y)
	}
	spec1 = starspec(
			mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(2,2), 
					statevar = c("y","s")[1], ylags = 1, yfun = fun))
	solver.control=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, trace=1,method="BFGS")
	fit1 = starfit(spec1, data = dx[1:521], out.sample = 10, solver = "strategy", solver.control = solver.control, n=5)
	
	dxv = runMean(abs(dx), n=4)
	dxv[1:3] = c(as.numeric(abs(dx[1,1])), as.numeric(mean(abs(dx[1:2,1]))), 
			as.numeric(mean(abs(dx[1:3,1]))))
	dxv = lag(dxv, 1)
	dxv[1]=0
	spec2 = starspec(
			mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(2,2), 
					statevar = c("y","s")[2], s = dxv[1:521]))
	solver.control=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, trace=1,method="BFGS")
	fit2 = starfit(spec2, data = dx[1:521], out.sample = 10, solver = "strategy", solver.control = solver.control, n=5)
	
	options(width=150)
	zz <- file("test1c-1.txt", open="wt")
	sink(zz)	
	print(rbind(coef(fit1), coef(fit2)))
	print(c(likelihood(fit1), likelihood(fit2)))
	sink(type="message")
	sink()
	close(zz)
	
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

twinkle.test1d = function(cluster = NULL){
	# structural breaks
	tic = Sys.time()
	set.seed(100)
	
	# 1-break (2 regimes)
	b1 = xts(c(rnorm(1000, 0.1, 0.2), rnorm(500, -0.02, 0.2)), as.Date(1:1500))
	x1 = xts(seq(0, 1, length.out=1500), index(b1))
	spec0 = starspec(
			mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(0,0), 
					statevar = c("y","s")[2], s = x1),
			variance.model = list(dynamic = FALSE), distribution.model = "norm")
	solver.control=list(maxit=17000, alpha=1, beta=0.4, gamma=1.4,trace=1, method="BFGS")
	fit0 = starfit(spec0, data = b1, solver = "strategy", solver.control = solver.control, n=4)
	
	options(width=150)
	zz <- file("test1d-1.txt", open="wt")
	sink(zz)	
	cat("\nTrue states are N(0.1, 0.2) and N(-0.02, 0.2)\n")
	print(coef(fit0))
	print(fit0)
	sink(type="message")
	sink()
	close(zz)
	
	jpeg("plot_1d-1.jpeg", width = 800, height = 1200, quality=100)
	par(mfrow=c(3,1))
	plot(as.numeric(states(fit0)[,2]), main = "Prob[state=2]", type="l", ylab="Probability", xlab="Time")
	abline(v=1000, untf  = TRUE, col=2)
	grid()
	plot(as.numeric(states(fit0, type="pmu")), main = "State Dynamics", type="l", ylab="Time Trend Dynamics", xlab="Time")
	abline(v=1000, untf  = TRUE, col=2)
	grid()
	plot(as.numeric(b1), main = "Fitted vs Actual", type="l", ylab="Value", xlab="Time")
	lines(as.numeric(fitted(fit0)), col="steelblue", lwd=2)
	abline(v=1000, untf  = TRUE, col=2)
	grid()
	dev.off()
	
	
	# 1-break (2 regimes mixture)
	b1 = xts(c(rnorm(1000, 0.1, 0.2), rnorm(500, 0.1, 0.1)), as.Date(1:1500))
	x1 = xts(seq(0, 1, length.out=1500), index(b1))
	spec1 = starspec(
			mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(0,0), 
					statevar = c("y","s")[2], s = x1),
			variance.model = list(dynamic = TRUE, model="mixture"), distribution.model = "norm")	
	solver.control=list(maxit=17000, alpha=1, beta=0.4, gamma=1.4,trace=1, method="BFGS")
	fit1 = starfit(spec1, data = b1, solver = "strategy", solver.control = solver.control, n=6)
	
	zz <- file("test1d-1x.txt", open="wt")
	sink(zz)
	cat("\nTrue states are N(0.1, 0.2) and N(0.1, 0.1)\n")
	print(coef(fit1))
	sink(type="message")
	sink()
	close(zz)
	
	jpeg("plot_1d-1x.jpeg", width = 1000, height = 1200, quality=100)
	par(mfrow=c(3,1))
	plot(as.numeric(states(fit1)[,2]), main = "Prob[state=2]", type="l", ylab="Probability", xlab="Time")
	abline(v=1000, untf  = TRUE, col=2)
	grid()
	plot(as.numeric(states(fit1, type="pmu")), main = "State Dynamics", type="l", ylab="Time Trend Dynamics", xlab="Time")
	abline(v=1000, untf  = TRUE, col=2)
	grid()
	plot(as.numeric(b1), main = "Fitted vs Actual", type="l", ylab="Value", xlab="Time")
	lines(as.numeric(fitted(fit1)), col="steelblue", lwd=2)
	abline(v=1000, untf  = TRUE, col=2)
	grid()
	dev.off()
	
	
	# 2-breaks (3 regimes)
	set.seed(115)
	b2 = xts(c(rnorm(1500, 0.1, 0.1), rnorm(1500, -0.2, 0.1), rnorm(1500, 0.2, 0.1)), as.Date(1:4500))
	x2 = xts(seq(0, 10, length.out=4500), index(b2))
	spec2 = starspec(
			mean.model = list(states = 3, include.intercept = c(1,1,1), arOrder = c(0,0,0), 
					statevar = c("y","s")[2], s = x2))
	solver.control=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, trace=1,method="BFGS")
	fit2 = starfit(spec2, data = b2, solver = "strategy", solver.control = solver.control, n=6)
	
	zz <- file("test1d-2.txt", open="wt")
	sink(zz)
	cat("\nTrue states are N(0.1, 0.1), N(-0.2, 0.1) and N(0.2, 0.1)\n")
	print(coef(fit2))
	sink(type="message")
	sink()
	close(zz)
	
	idxx = coef(fit2)[1:3]
	idxx = match(c(0.1, -0.2, 0.2), round(idxx,2))
	jpeg("plot_1d-2.jpeg", width = 800, height = 1200, quality=100)
	mat = matrix(c(1,1,2,3,4,4),3,2,byrow=TRUE)
	nf=layout(mat)
	layout.show(nf)
	plot(as.numeric(states(fit2)[,idxx[1]]), ylim=c(0,1), main = "Prob[state={1,2,3}]", type="l", ylab="Probability", xlab="Time")
	lines(as.numeric(states(fit2)[,idxx[2]]), col=2, lty=2)
	lines(as.numeric(states(fit2)[,idxx[3]]), col=3, lty=3)
	abline(v=1500)
	abline(v=3000)
	grid()
	legend("topleft", c("[state=1]", "[state=2]", "[state=3]"), col=1:3, lty=c(1,2,3), bty="n")
	id1 = ifelse(idxx[1]>idxx[2], 1, 2)
	id2 = ifelse(idxx[1]>idxx[2], 2, 1)
	plot(as.numeric(states(fit2, type="pmu")[,id1]), main = "State Dynamics[1]", type="l", ylab="Time Trend Dynamics", xlab="Time")
	idx = min(which(as.numeric(states(fit2, type="pmu")[,id1])>0))
	abline(v=idx, h=0, col=2)
	grid()
	plot(as.numeric(states(fit2, type="pmu")[,id2]), main = "State Dynamics[2]", type="l", ylab="Time Trend Dynamics", xlab="Time")
	idx = min(which(as.numeric(states(fit2, type="pmu")[,id2])>0))
	abline(v=idx, h=0, col=2)
	grid()
	plot(as.numeric(b2), main = "Fitted vs Actual", type="l", ylab="Value", xlab="Time")
	lines(as.numeric(fitted(fit2)), col="steelblue", lwd=2)
	grid()
	dev.off()
	
	
	# 3-breaks (4 regimes)
	set.seed(125)
	b3 = xts(c(rnorm(1000, 0.1, 0.1), rnorm(1000, -0.2, 0.1), rnorm(1000, 0.2, 0.1), rnorm(1000, -0.1, 0.1)), as.Date(1:4000))
	x3 = xts(seq(0, 1, length.out=4000), index(b3))
	spec3 = starspec(
			mean.model = list(states = 4, include.intercept = c(1,1,1,1), arOrder = c(0,0,0,0), 
					statevar = c("y","s")[2], ylags = 1, s = exp(x3)))
	solver.control=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, trace=1,method="BFGS")
	fit3 = starfit(spec3, data = b3, solver = "strategy", solver.control = solver.control, n=6)
	
	
	options(width=150)
	zz <- file("test1d-3.txt", open="wt")
	sink(zz)
	cat("\nTrue states are N(0.1, 0.1), N(-0.2, 0.1), N(0.2, 0.1) and N(-0.1, 0.1)\n")
	print(round(coef(fit3),3))
	print(fit3)
	sink(type="message")
	sink()
	close(zz)
	
	
	idxx = coef(fit3)[1:4]
	idxx = match(c(0.1, -0.2, 0.2, -0.1), round(idxx,2))
	
	jpeg("plot_1d-3.jpeg", width = 1000, height = 1200, quality=100)
	mat = matrix(c(1,1,1,2,3,4,5,5,5),3,3,byrow=TRUE)
	nf=layout(mat)
	layout.show(nf)
	plot(as.numeric(states(fit3)[,idxx[1]]), ylim=c(0,1), main = "Prob[state={1,2,3,4}]", type="l", ylab="Probability", xlab="Time")
	lines(as.numeric(states(fit3)[,idxx[2]]), col=2, lty=2)
	lines(as.numeric(states(fit3)[,idxx[3]]), col=3, lty=2)
	lines(as.numeric(states(fit3)[,idxx[4]]), col=4, lty=2)
	abline(v=1000)
	abline(v=2000)
	abline(v=3000)
	grid()
	legend("topleft", c("[state=1]", "[state=2]", "[state=3]", "[state=4]"), col=c(3,4,1,2), lty=c(2,2,1,2), bty="n")
	plot(as.numeric(states(fit3, type="pmu")[,id1]), main = "State Dynamics[1]", type="l", ylab="Time Trend Dynamics", xlab="Time")
	idx = min(which(as.numeric(states(fit3, type="pmu")[,id1])>0))
	abline(v=idx, h=0, col=2)
	grid()
	plot(as.numeric(states(fit3, type="pmu")[,id2]), main = "State Dynamics[2]", type="l", ylab="Time Trend Dynamics", xlab="Time")
	idx = min(which(as.numeric(states(fit3, type="pmu")[,id2])>0))
	abline(v=idx, h=0, col=2)
	grid()
	plot(as.numeric(states(fit3, type="pmu")[,id3]), main = "State Dynamics[3]", type="l", ylab="Time Trend Dynamics", xlab="Time")
	idx = min(which(as.numeric(states(fit3, type="pmu")[,id3])>0))
	abline(v=idx, h=0, col=2)
	grid()
	plot(as.numeric(b3), main = "Fitted vs Actual", type="l", ylab="Value", xlab="Time")
	lines(as.numeric(fitted(fit3)), col="steelblue", lwd=2)
	grid()
	dev.off()
	
	
	# 2-breaks (3 regimes)
	set.seed(115)
	b2 = xts(c(rnorm(1500, 0.1, 0.2), rnorm(1500, -0.2, 0.1), rnorm(1500, 0.2, 0.25)), as.Date(1:4500))
	x2 = xts(seq(0, 10, length.out=4500), index(b2))
	spec4 = starspec(
			mean.model = list(states = 3, include.intercept = c(1,1,1), arOrder = c(0,0,0), 
					statevar = c("y","s")[2], s = x2), variance.model=list(dynamic=TRUE, model="mixture"))
	solver.control=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, trace=1,method="BFGS")
	fit4 = starfit(spec4, data = b2, solver = "strategy", solver.control = solver.control, n=6)
	
	zz <- file("test1d-4.txt", open="wt")
	sink(zz)
	cat("\nTrue states are N(0.1, 0.2), N(-0.2, 0.1) and N(0.2, 0.25)\n")
	print(coef(fit4))
	print(fit4)
	sink(type="message")
	sink()
	close(zz)
	
	idxx = coef(fit4)[1:3]
	idxx = match(c(0.1, -0.2, 0.2), round(idxx,2))
	jpeg("plot_1d-4.jpeg", width = 800, height = 1200, quality=100)
	mat = matrix(c(1,1,2,3,4,4),3,2,byrow=TRUE)
	nf=layout(mat)
	layout.show(nf)
	plot(as.numeric(states(fit2)[,idxx[1]]), ylim=c(0,1), main = "Prob[state={1,2,3}]", type="l", ylab="Probability", xlab="Time")
	lines(as.numeric(states(fit2)[,idxx[2]]), col=2, lty=2)
	lines(as.numeric(states(fit2)[,idxx[3]]), col=3, lty=3)
	abline(v=1500)
	abline(v=3000)
	grid()
	legend("topleft", c("[state=1]", "[state=2]", "[state=3]"), col=1:3, lty=c(1,2,3), bty="n")
	id1 = ifelse(idxx[1]>idxx[2], 1, 2)
	id2 = ifelse(idxx[1]>idxx[2], 2, 1)
	plot(as.numeric(states(fit2, type="pmu")[,id1]), main = "State Dynamics[1]", type="l", ylab="Time Trend Dynamics", xlab="Time")
	idx = min(which(as.numeric(states(fit2, type="pmu")[,id1])>0))
	abline(v=idx, h=0, col=2)
	grid()
	plot(as.numeric(states(fit2, type="pmu")[,id2]), main = "State Dynamics[2]", type="l", ylab="Time Trend Dynamics", xlab="Time")
	idx = min(which(as.numeric(states(fit2, type="pmu")[,id2])>0))
	abline(v=idx, h=0, col=2)
	grid()
	plot(as.numeric(b2), main = "Fitted vs Actual", type="l", ylab="Value", xlab="Time")
	lines(as.numeric(fitted(fit2)), col="steelblue", lwd=2)
	grid()
	dev.off()
	
	
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}



twinkle.test1e = function(cluster = NULL){
	# Zivot Book benchmark comparison
	tic = Sys.time()
	require(tsDyn)
	require(quantmod)
	data(ndx)
	# Translating part of the old S+ code in book example using quantmod
	ndx.ret2 = ROC(Cl(ndx), na.pad=FALSE)^2
	ndx.rvol = sqrt(apply.weekly(ndx.ret2, FUN="sum"))
	colnames(ndx.rvol) = "RVOL"

	jpeg("plot_1e-1.jpeg", width = 800, height = 1200, quality=100)
	par(mfrow=c(2,2))
	plot(ndx.rvol, main="RVOL")
	plot(log(ndx.rvol), main="Log RVOL")
	ndx.acf = acf(log(as.numeric(ndx.rvol)), 25, main="ACF log(ndx.rvol)")
	ndx.pacf = pacf(log(as.numeric(ndx.rvol)), 25, main="PACF log(ndx.rvol)")
	dev.off()
	
	# Luukkonen Test
	spec = starspec(mean.model=list(states=2,arOrder=c(2,2), statevar='y',ylags=1))
	tmp = nonlinearTest(spec, data=log(ndx.rvol))
	twinkle.ntest1 = c("statistic" = tmp$F.statistic, pvalue = tmp$F.pvalue)
	zivotbook.ntest1 = c("statistic" = 3.7068, "pvalue"=0.0014)
	
	spec = starspec(mean.model=list(states=2,arOrder=c(2,2), statevar='y',ylags=2))
	tmp = nonlinearTest(spec, data=log(ndx.rvol))
	twinkle.ntest2 = c("statistic" = tmp$F.statistic, pvalue = tmp$F.pvalue)
	zivotbook.ntest2 = c("statistic" = 2.3204, "pvalue" = 0.0333)
	
	# Luukkonen Test (Robust)
	
	zz <- file("test1e-1.txt", open="wt")
	sink(zz)
	print(data.frame(twinkle = round(twinkle.ntest1,4), benchmark = zivotbook.ntest1))
	print(data.frame(twinkle = round(twinkle.ntest2,4), benchmark = zivotbook.ntest2))
	sink(type="message")
	sink()
	close(zz)
	
	
	# Comparison between Zivot, tsDyn and twinkle
	spec = starspec(mean.model=list(states=2,arOrder=c(2,2), statevar='y',ylags=1))
	ctrl = list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, trace=1,method="BFGS",n.restarts=4)
	fit = starfit(spec, log(ndx.rvol), solver="msoptim", solver.control = ctrl)
	
	zivotbook.cf = c("s1.phi0" = -2.668, "s1.phi1" = -0.396, "s1.phi2" = 0.216,"s2.phi0" = -3.729, "s2.phi1" = -0.221, "s2.phi2" = 0.205, 
			"s1.gamma" = 1.608, "s1.c" = -2.845, "sigma" = 0.415, "logL" = -158.863)
	
	specf0 = spec
	setfixed(specf0)<-as.list(c("s1.phi0" = -2.668, "s1.phi1" = -0.396, "s1.phi2" = 0.216,"s2.phi0" = -3.729, "s2.phi1" = -0.221, "s2.phi2" = 0.205, 
					"s1.gamma" = 1.608, "s1.c" = -2.845, "sigma" = 0.415))
	tmp = starfilter(specf0, log(ndx.rvol))
	
	
	twinkle.cf = c(coef(fit), "logL" = likelihood(fit))
	
	
	mod.lstar <- lstar(as.numeric(log(ndx.rvol)), steps=1, d=1, mL = 2, mH=2, thDelay=0, control=list(maxit=3000), 
			starting.control=list(gammaInt=c(0,5)))
	tmp = coef(mod.lstar)
	tmp = unname(tmp)
	
	tsdyn.cf = c("s2.phi0" = tmp[1], "s2.phi1" = tmp[2], "s2.phi2" = tmp[3],"s1.phi0" = tmp[1]+tmp[4], "s1.phi1" = tmp[2]+tmp[5], 
			"s1.phi2" = tmp[3]+tmp[6], "s1.gamma" = tmp[7], "s1.c" = tmp[8], "sigma" = sd(residuals(mod.lstar), na.rm=TRUE), "logL" = NA)
	tsdyn.cf = tsdyn.cf[names(twinkle.cf)]
	specf1 = spec
	setfixed(specf1)<-as.list(tsdyn.cf[-10])
	tmp = starfilter(specf1, log(ndx.rvol))
	tsdyn.cf[10] = likelihood(tmp)
	
	

	zz <- file("test1e-2.txt", open="wt")
	sink(zz)
	print(data.frame(twinkle = round(twinkle.cf,4), benchmark = zivotbook.cf, tsDyn =  round(tsdyn.cf,4)))
	sink(type="message")
	sink()
	close(zz)
}

twinkle.test1f = function(cluster = NULL){
	# Test Forecasts
	tic = Sys.time()
	data(forex)
	fx = na.locf(forex, fromLast = TRUE)
	fxw = fx[which(weekdays(index(forex))=="Wednesday"),4]
	dx = ROC(fxw, na.pad=FALSE)*100
	fun = function(x){
		x = as.numeric(x)
		y = runMean(abs(x), n=4)
		y[1:3] = c(abs(x[1]), mean(abs(x[1:2])), mean(abs(x[1:3])))
		return(y)
	}
	spec1 = starspec(
			mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(0,2), 
					statevar = c("y","x")[1], ylags = 1, yfun = fun))
	solver.control=list(maxit=150000, alpha=1, beta=0.4, gamma=1.2, reltol=1e-12, trace=1,method="BFGS")
	fit1 = starfit(spec1, data = dx[1:521], out.sample = 25, solver = "strategy", solver.control = solver.control, n=5)
	
	f1x1 = starforecast(fit1, n.roll=2, n.ahead=20, method="an.parametric", mc.sims = 10000)
	f1x2 = starforecast(fit1, n.roll=2, n.ahead=20, method="an.kernel", mc.sims = 10000)
	f1x3 = starforecast(fit1, n.roll=2, n.ahead=20, method="mc.empirical", mc.sims = 10000)
	f1x4 = starforecast(fit1, n.roll=2, n.ahead=20, method="mc.parametric", mc.sims = 10000)
	f1x5 = starforecast(fit1, n.roll=2, n.ahead=20, method="mc.kernel", mc.sims = 10000)
	
	jpeg("plot_1f-2.jpeg", width = 800, height = 1200, quality=100)
	par(mfrow=c(2,4))
	plot(f1x1, roll = 0, cex.main=0.95)
	plot(f1x1, roll = 1, cex.main=0.95)
	plot(f1x2, cex.main=0.95)
	plot(f1x3, cex.main=0.95)
	plot(f1x4, cex.main=0.95)
	plot(f1x5, roll = 0, cex.main=0.95)
	plot(f1x5, roll = 1, cex.main=0.95)
	dev.off()
	
	spec2 = starspec(
			mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(0,2), 
					statevar = c("y","x")[1], ylags = 1, yfun = fun),
			variance.model = list(dynamic = TRUE, model = "sGARCH"), 
			distribution.model = "norm")

	solver.control=list(maxit=150000, alpha=1, beta=0.4, gamma=1.3, reltol=1e-12, trace=1,method="L-BFGS-B",n.restarts=5, rseed=110)
	fit2 = starfit(spec2, data = dx[1:521], out.sample = 25, solver = "msoptim", solver.control = solver.control)	
	
	f2x1 = starforecast(fit2, n.roll=2, n.ahead=20, method="an.parametric", mc.sims = 10000)
	f2x2 = starforecast(fit2, n.roll=2, n.ahead=20, method="an.kernel", mc.sims = 10000)
	f2x3 = starforecast(fit2, n.roll=2, n.ahead=20, method="mc.empirical", mc.sims = 10000)
	f2x4 = starforecast(fit2, n.roll=2, n.ahead=20, method="mc.parametric", mc.sims = 10000)
	f2x5 = starforecast(fit2, n.roll=2, n.ahead=20, method="mc.kernel", mc.sims = 10000)
	
	jpeg("plot_1f-2.jpeg", width = 800, height = 1200, quality=100)
	par(mfrow=c(2,4))
	plot(f2x1, roll = 0, cex.main=0.95)
	plot(f2x1, roll = 1, cex.main=0.95)
	plot(f2x2, cex.main=0.95)
	plot(f2x3, cex.main=0.95)
	plot(f2x4, cex.main=0.95)
	plot(f2x5, roll = 0, cex.main=0.95)
	plot(f2x5, roll = 1, cex.main=0.95)
	dev.off()
	
	spec3 = starspec(
			mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(0,2), 
					statevar = c("y","x")[1], ylags = 1, yfun = fun),
			variance.model = list(dynamic = TRUE, model = "mixture"), 
			distribution.model = "norm")
	solver.control=list(maxit=150000, alpha=1, beta=0.4, gamma=1.2, reltol=1e-12, trace=1,method="BFGS",n.restarts=10, rseed=10)
	fit3 = starfit(spec3, data = dx[1:521], out.sample = 25, solver = "msoptim", solver.control = solver.control)
	
	spec3x = spec3
	setfixed(spec3x)<-as.list(coef(fit3))
	filt3 = starfilter(spec3x, data = dx[1:521], n.old = 521-25)
	head(cbind(sigma(fit3), sigma(filt3)))
	
	f3x1 = starforecast(fit3, n.roll=2, n.ahead=20, method="an.parametric", mc.sims = 10000)
	f3x2 = starforecast(fit3, n.roll=2, n.ahead=20, method="an.kernel", mc.sims = 10000)
	f3x3 = starforecast(fit3, n.roll=2, n.ahead=20, method="mc.empirical", mc.sims = 10000)
	f3x4 = starforecast(fit3, n.roll=2, n.ahead=20, method="mc.parametric", mc.sims = 10000)
	f3x5 = starforecast(fit3, n.roll=2, n.ahead=20, method="mc.kernel", mc.sims = 10000)
	
	jpeg("plot_1f-3.jpeg", width = 800, height = 1200, quality=100)
	par(mfrow=c(2,4))
	plot(f3x1, roll = 0, cex.main=0.95)
	plot(f3x1, roll = 1, cex.main=0.95)
	plot(f3x2, cex.main=0.95)
	plot(f3x3, cex.main=0.95)
	plot(f3x4, cex.main=0.95)
	plot(f3x5, roll = 0, cex.main=0.95)
	plot(f3x5, roll = 1, cex.main=0.95)
	dev.off()
	
	# methods on forecasts:
	jpeg("plot_1f-4.jpeg", width = 800, height = 1200, quality=100)
	par(mfrow=c(2,1))
	plot(sigma(f2x1)[,1], type="l", ylim=c(1,2.4), main="10-ahead Sigma Forecast",
			ylab="sigma", xlab="Time")
	lines(sigma(f3x1)[,1], col=2)
	legend("topright", c("GARCH","NormMix(2)"), col=1:2, lty=c(1,1), bty="n")
	plot(states(f2x1)[,1,2], type="l", main = "10-ahead Prob[State=2] Forecast", ylab="Probability",
			xlab="Time")
	lines(states(f3x1)[,1,2], col=2)
	legend("bottomright", c("GARCH","NormMix(2)"), col=1:2, lty=c(1,1), bty="n")
	dev.off()
	
	
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}


twinkle.test1h = function(cluster = NULL){
	# Rolling Forecasts
	tic = Sys.time()
	load("./data/forex.rda")
	fx = na.locf(forex, fromLast = TRUE)
	fxw = fx[which(weekdays(index(forex))=="Wednesday"),4]
	dx = ROC(fxw, na.pad=FALSE)*100
	
	
	fun = function(x){
		x = as.numeric(x)
		y = runMean(abs(x), n=4)
		y[1:3] = c(abs(x[1]), mean(abs(x[1:2])), mean(abs(x[1:3])))
		return(y)
	}
	spec1 = starspec(
			mean.model = list(regimes = 2, include.intercept = c(1,1), arOrder = c(1,1), 
					type = c("y","x")[1], dstar = FALSE, lags = 1, 
					external.regressors = NULL, fun = fun),
			variance.model = list(dynamic = FALSE, model = "sGARCH", 
					garchOrder = c(1, 1), submodel = NULL, external.regressors = NULL, 
					variance.targeting = FALSE), 
			distribution.model = "norm")
	solver.control=list(maxit=5000, alpha=1, beta=0.4, gamma=2, reltol=1e-12, trace=1,method="BFGS",n.restarts=5)
	roll1 = rollstar(spec1, dx, n.ahead = 1, n.start = 500, refit.every = 25, refit.window = c("recursive", "moving")[1], 
		solver = "msoptim", fit.control = list(), solver.control = solver.control)
	d1 = as.data.frame(roll1, which = "VaR")
	
	VaRTest(0.05, actual = d1[,3], VaR = d1[,2])
	VaRTest(0.01, actual = d1[,3], VaR = d1[,1])

	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}



twinkle.test1i = function(cluster = NULL){
	# Deterministic Simulation for gauging the dynamics of the system (equilibrium, stability etc)
	tic = Sys.time()
	load("./data/forex.rda")
	fx = na.locf(forex, fromLast = TRUE)
	fxw = fx[which(weekdays(index(forex))=="Wednesday"),4]
	dx = ROC(fxw, na.pad=FALSE)*100
	fun = function(x){
		x = as.numeric(x)
		y = runMean(abs(x), n=4)
		y[1:3] = c(abs(x[1]), mean(abs(x[1:2])), mean(abs(x[1:3])))
		return(y)
	}
	spec = starspec(
			mean.model = list(regimes = 2, include.intercept = c(1,1), arOrder = c(2,2), 
					type = c("y","x")[1], dstar = FALSE, lags = 1, 
					external.regressors = NULL, fun = fun), 
			distribution.model = "norm")
	solver.control=list(maxit=5000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, trace=1,method="L-BFGS-B",n.restarts=5)
	fit = starfit(spec, data = dx[1:521], solver = "msoptim", solver.control = solver.control)
	
	plot(fit)
	transfun2d(fit)
	
	# Even though there are no shocks, so that we will be in the previous state
	# the state conditional 'mean' fluctuates in an oscillatory fashion around
	# an attractor
	sim1 = starsim(fit, n.sim=250, n.start=0, custom.dist = list(name = "sample", 
					distfit = matrix(0, ncol=1, nrow=250)))
	
	par(mfrow=c(2,2))
	plot(sim1@simulation$seriesSim[,1], type="l")
	plot(sim1@simulation$condmSim[[1]][,1], type="l")
	abline(h=coef(fit)["s1.phi0"], col=2)
	abline(h=mean(sim1@simulation$condmSim[[1]][,1]), col=3)
	plot(sim1@simulation$condmSim[[1]][,2], type="l")
	abline(h=coef(fit)["s2.phi0"], col=2)
	abline(h=mean(sim1@simulation$condmSim[[1]][,2]), col=3)
	plot(sim1@simulation$probSim[[1]][,2], type="l")
	
	sim2 = starsim(fit, n.sim=250, n.start=0, custom.dist = list(name = "sample", 
					distfit = matrix(0, ncol=1, nrow=250)), prereturns=c(0,0))
	
	par(mfrow=c(2,2))
	plot(sim2@simulation$seriesSim[,1], type="l")
	plot(sim2@simulation$condmSim[[1]][,1], type="l")
	abline(h=coef(fit)["s1.phi0"], col=2)
	abline(h=mean(sim2@simulation$condmSim[[1]][,1]), col=3)
	plot(sim2@simulation$condmSim[[1]][,2], type="l")
	abline(h=coef(fit)["s2.phi0"], col=2)
	abline(h=mean(sim2@simulation$condmSim[[1]][,2]), col=3)
	plot(sim2@simulation$probSim[[1]][,2], type="l", ylim=c(0,1))
	
	par(mfrow=c(2,1))
	plot(sim1@simulation$condmSim[[1]][,1], type="l")
	lines(sim2@simulation$condmSim[[1]][,1], col=2)
	plot(sim1@simulation$condmSim[[1]][,2], type="l")
	lines(sim2@simulation$condmSim[[1]][,2], col=2)
	
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}