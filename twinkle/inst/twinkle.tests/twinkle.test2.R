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

# Simulation Tests:
# -- path versus sim
# -- unconditional
# -- methods and plots
# -- GIRF

# static variance
twinkle.test2a = function(cluster = NULL){
	tic = Sys.time()
	require(twinkle)
	data(forex)	
	
	# 1. Equivalence of sim and path methods for self-exciting model with yfun
	
	# load("./data/forex.rda")
	fx = na.locf(forex, fromLast = TRUE)
	fxw = fx[which(weekdays(index(forex))=="Wednesday"),4]
	dx = ROC(fxw, na.pad=FALSE)*100
	fun = function(x){
		x = as.numeric(x)
		N = length(x)
		if(N<4){
			y = abs(x)	
		} else{
			y = runMean(abs(x), n=4)
			y[1:3] = c(abs(x[1]), mean(abs(x[1:2])), mean(abs(x[1:3])))
		}
		return(y)
	}
	
	arorder=list(c(0,0),c(0,1), c(1,0), c(1,1), c(2,1))
	int = list(c(0,0), c(1,0))
	test = matrix(NA, 2, 5)
	for(i in 1:2){
		for(j in 1:5){
			spec = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]], 
							statevar = c("y","s")[1], ylags = 1, yfun = fun))
			solver.control=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE)
			fit = starfit(spec, data = dx[1:521], solver = "strategy", solver.control = solver.control, n=5)
			sim  = starsim(fit, n.sim=5000)
			specx = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]],
							statevar = c("y","s")[1], ylags = 1, yfun = fun))
			setfixed(specx)<-as.list(coef(fit))
			path = starpath(specx, n.sim=5000, prereturns = tail(dx[1:521], 6), 
					custom.dist=list(name="sample", distfit=sim@simulation$residSim/coef(fit)["sigma"]))			
			test[i,j] = as.integer(all.equal(fitted(sim), fitted(path)))
			#print(test[i,j])
		}
	}
	zz <- file("test2a-1.txt", open="wt")
	sink(zz)
	print(test)
	sink(type="message")
	sink()
	close(zz)
	
	
	# 2. Equivalence of sim and path methods for self-exciting model without yfun
	arorder=list(c(0,0),c(0,1), c(1,0), c(1,1), c(2,1))
	int = list(c(0,0), c(1,0))
	test = matrix(NA, 2, 5)
	for(i in 1:2){
		for(j in 1:5){
			spec = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]], 
							statevar = c("y","s")[1], ylags = 1))
			solver.control=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE)
			fit = starfit(spec, data = dx[1:521], solver = "strategy", solver.control = solver.control, n=5)
			sim  = starsim(fit, n.sim=5000)
			specx = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]],
							statevar = c("y","s")[1], ylags = 1))
			setfixed(specx)<-as.list(coef(fit))
			path = starpath(specx, n.sim=5000, prereturns = tail(dx[1:521], 6), 
					custom.dist=list(name="sample", distfit=sim@simulation$residSim/coef(fit)["sigma"]))
			test[i,j] = as.integer(all.equal(fitted(sim), fitted(path)))
			#print(test[i,j])
		}
	}
	zz <- file("test2a-2.txt", open="wt")
	sink(zz)
	print(test)
	sink(type="message")
	sink()
	close(zz)
	
	
	# 3. Equivalence of sim and path methods for s-var model
	dxv = runMean(abs(dx), n=4)
	dxv[1:3] = c(as.numeric(abs(dx[1,1])), as.numeric(mean(abs(dx[1:2,1]))), 
			as.numeric(mean(abs(dx[1:3,1]))))
	dxv = lag(dxv, 1)
	dxv[1] = 0
	arorder=list(c(0,0),c(0,1), c(1,0), c(1,1), c(2,1))
	int = list(c(0,0), c(1,0))
	test = matrix(NA, 2, 5)
	
	for(i in 1:2){
		for(j in 1:5){
			spec = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]], 
							statevar = c("y","s")[2], s = dxv[1:521]))
			solver.control=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE)
			fit = starfit(spec, data = dx[1:521], solver = "strategy", solver.control = solver.control, n=5)
			ssim = sample(as.numeric(dxv), 5000, replace=TRUE)
			sim  = starsim(fit, n.sim=5000, ssim = list(matrix(ssim, ncol=1)))
			specx = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]], 
							statevar = c("y","s")[2], s = dxv[500:521]))
			setfixed(specx)<-as.list(coef(fit))
			path = starpath(specx, n.sim=5000, prereturns = tail(dx[1:521], 6), 
					custom.dist=list(name="sample", distfit=sim@simulation$residSim/coef(fit)["sigma"]),
					ssim= list(matrix(ssim, ncol=1)))
			test[i,j] = as.integer(all.equal(fitted(sim), fitted(path)))
		}
	}
	zz <- file("test2a-3.txt", open="wt")
	sink(zz)
	print(test)
	sink(type="message")
	sink()
	close(zz)
	
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

# mixture variance
twinkle.test2b = function(cluster = NULL){
	tic = Sys.time()
	require(twinkle)
	data(forex)	
	
	# 1. Equivalence of sim and path methods for self-exciting model with yfun
	
	# load("./data/forex.rda")
	fx = na.locf(forex, fromLast = TRUE)
	fxw = fx[which(weekdays(index(forex))=="Wednesday"),4]
	dx = ROC(fxw, na.pad=FALSE)*100
	fun = function(x){
		x = as.numeric(x)
		N = length(x)
		if(N<4){
			y = abs(x)	
		} else{
			y = runMean(abs(x), n=4)
			y[1:3] = c(abs(x[1]), mean(abs(x[1:2])), mean(abs(x[1:3])))
		}
		return(y)
	}
	
	arorder=list(c(0,0),c(0,1), c(1,0), c(1,1), c(2,1))
	int = list(c(0,0), c(1,0))
	test = matrix(NA, 2, 5)
	for(i in 1:2){
		for(j in 1:5){
			spec = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]], 
							statevar = c("y","s")[1], ylags = 1, yfun = fun),
					variance.model=list(dynamic=TRUE, model="mixture"))
			solver.control=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE)
			fit = starfit(spec, data = dx[1:521], solver = "strategy", solver.control = solver.control, n=5)
			sim  = starsim(fit, n.sim=5000)
			specx = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]],
							statevar = c("y","s")[1], ylags = 1, yfun = fun),
					variance.model=list(dynamic=TRUE, model="mixture"))
			setfixed(specx)<-as.list(coef(fit))
			path = starpath(specx, n.sim=5000, prereturns = tail(dx[1:521], 12), preresiduals = tail(residuals(fit),2),
					custom.dist=list(name="sample", distfit=sim@simulation$residSim/sigma(sim)))
			test[i,j] = as.integer(all.equal(fitted(sim), fitted(path)))
			print(test[i,j])
		}
	}
	zz <- file("test2b-1.txt", open="wt")
	sink(zz)
	print(test)
	sink(type="message")
	sink()
	close(zz)
	
	
	# 2. Equivalence of sim and path methods for self-exciting model without yfun
	arorder=list(c(0,0),c(0,1), c(1,0), c(1,1), c(2,1))
	int = list(c(0,0), c(1,0))
	test = matrix(NA, 2, 5)
	for(i in 1:2){
		for(j in 1:5){
			spec = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]], 
							statevar = c("y","s")[1], ylags = 1),
					variance.model=list(dynamic=TRUE, model="mixture"))
			solver.control=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE)
			fit = starfit(spec, data = dx[1:521], solver = "strategy", solver.control = solver.control, n=5)
			sim  = starsim(fit, n.sim=5000)
			specx = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]],
							statevar = c("y","s")[1], ylags = 1),
					variance.model=list(dynamic=TRUE, model="mixture"))
			setfixed(specx)<-as.list(coef(fit))
			path = starpath(specx, n.sim=5000, prereturns = tail(dx[1:521], 6), preresiduals = tail(residuals(fit),2),
					custom.dist=list(name="sample", distfit=sim@simulation$residSim/sigma(sim)))
			test[i,j] = as.integer(all.equal(fitted(sim), fitted(path)))
			print(test[i,j])
		}
	}
	zz <- file("test2b-2.txt", open="wt")
	sink(zz)
	print(test)
	sink(type="message")
	sink()
	close(zz)
	
	
	# 3. Equivalence of sim and path methods for s-var model
	dxv = runMean(abs(dx), n=4)
	dxv[1:3] = c(as.numeric(abs(dx[1,1])), as.numeric(mean(abs(dx[1:2,1]))), 
			as.numeric(mean(abs(dx[1:3,1]))))
	dxv = lag(dxv, 1)
	dxv[1] = 0
	arorder=list(c(0,0),c(0,1), c(1,0), c(1,1), c(2,1))
	int = list(c(0,0), c(1,0))
	test = matrix(NA, 2, 5)
	
	for(i in 1:2){
		for(j in 1:5){
			spec = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]], 
							statevar = c("y","s")[2], s = dxv[1:521]),
					variance.model=list(dynamic=TRUE, model="mixture"))
			solver.control=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE)
			fit = starfit(spec, data = dx[1:521], solver = "strategy", solver.control = solver.control, n=5)
			ssim = sample(as.numeric(dxv), 5000, replace=TRUE)
			sim  = starsim(fit, n.sim=5000, ssim = list(matrix(ssim, ncol=1)))
			specx = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]], 
							statevar = c("y","s")[2], s = dxv[500:521]),
					variance.model=list(dynamic=TRUE, model="mixture"))
			setfixed(specx)<-as.list(coef(fit))
			path = starpath(specx, n.sim=5000, prereturns = tail(dx[1:521], 6), 
					custom.dist=list(name="sample", distfit=sim@simulation$residSim/sigma(sim)),
					ssim=list(matrix(ssim, ncol=1)))
			test[i,j] = as.integer(all.equal(fitted(sim), fitted(path)))
		}
	}
	zz <- file("test2b-3.txt", open="wt")
	sink(zz)
	print(test)
	sink(type="message")
	sink()
	close(zz)
	
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

# GARCH variance
twinkle.test2c = function(cluster = NULL){
	tic = Sys.time()
	require(twinkle)
	data(forex)	
	
	# 1. Equivalence of sim and path methods for self-exciting model with yfun
	
	# load("./data/forex.rda")
	fx = na.locf(forex, fromLast = TRUE)
	fxw = fx[which(weekdays(index(forex))=="Wednesday"),4]
	dx = ROC(fxw, na.pad=FALSE)*100
	fun = function(x){
		x = as.numeric(x)
		N = length(x)
		if(N<4){
			y = abs(x)	
		} else{
			y = runMean(abs(x), n=4)
			y[1:3] = c(abs(x[1]), mean(abs(x[1:2])), mean(abs(x[1:3])))
		}
		return(y)
	}
	
	arorder=list(c(0,0),c(0,1), c(1,0), c(1,1), c(2,1))
	int = list(c(0,0), c(1,0))
	test = matrix(NA, 2, 5)
	for(i in 1:2){
		for(j in 1:5){
			spec = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]], 
							statevar = c("y","s")[1], ylags = 1, yfun = fun),
					variance.model=list(dynamic=TRUE, model="sGARCH"))
			solver.control=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, method="BFGS", parsearch=TRUE)
			fit = starfit(spec, data = dx[1:521], solver = "strategy", solver.control = solver.control, n=10)
			sim  = starsim(fit, n.sim=5000)
			specx = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]],
							statevar = c("y","s")[1], ylags = 1, yfun = fun),
					variance.model=list(dynamic=TRUE, model="sGARCH"))
			setfixed(specx)<-as.list(coef(fit))
			path = starpath(specx, n.sim=5000, prereturns = tail(dx[1:521], 12), preresiduals = tail(residuals(fit),2),
					presigma = tail(sigma(fit), 2),
					custom.dist=list(name="sample", distfit=sim@simulation$residSim/sigma(sim)))
			test[i,j] = as.integer(all.equal(fitted(sim), fitted(path)))
			print(test[i,j])
		}
	}
	zz <- file("test2c-1.txt", open="wt")
	sink(zz)
	print(test)
	sink(type="message")
	sink()
	close(zz)
	
	
	# 2. Equivalence of sim and path methods for self-exciting model without yfun
	arorder=list(c(0,0),c(0,1), c(1,0), c(1,1), c(2,1))
	int = list(c(0,0), c(1,0))
	test = matrix(NA, 2, 5)
	for(i in 1:2){
		for(j in 1:5){
			spec = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]], 
							statevar = c("y","s")[1], ylags = 1),
					variance.model=list(dynamic=TRUE, model="sGARCH"))
			# set reltol to a high number in order to converge quickly since
			solver.control=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-6, method="BFGS", parsearch=TRUE)
			fit = starfit(spec, data = dx[1:521], solver = "optim", solver.control = solver.control, n=5)
			sim  = starsim(fit, n.sim=5000)
			specx = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]],
							statevar = c("y","s")[1], ylags = 1),
					variance.model=list(dynamic=TRUE, model="sGARCH"))
			setfixed(specx)<-as.list(coef(fit))
			path = starpath(specx, n.sim=5000, prereturns = tail(dx[1:521], 6), preresiduals = tail(residuals(fit),2),
					presigma = tail(sigma(fit), 2),
					custom.dist=list(name="sample", distfit=sim@simulation$residSim/sigma(sim)))
			test[i,j] = as.integer(all.equal(fitted(sim), fitted(path)))
			print(test[i,j])
		}
	}
	zz <- file("test2c-2.txt", open="wt")
	sink(zz)
	print(test)
	sink(type="message")
	sink()
	close(zz)
	
	
	# 3. Equivalence of sim and path methods for s-var model
	dxv = runMean(abs(dx), n=4)
	dxv[1:3] = c(as.numeric(abs(dx[1,1])), as.numeric(mean(abs(dx[1:2,1]))), 
			as.numeric(mean(abs(dx[1:3,1]))))
	dxv = lag(dxv, 1)
	dxv[1] = 0
	arorder=list(c(0,0),c(0,1), c(1,0), c(1,1), c(2,1))
	int = list(c(0,0), c(1,0))
	test = matrix(NA, 2, 5)
	
	for(i in 1:2){
		for(j in 1:5){
			spec = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]], 
							statevar = c("y","s")[2], s = dxv[1:521]),
					variance.model=list(dynamic=TRUE, model="sGARCH"))
			solver.control=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-8, method="BFGS", parsearch=TRUE)
			fit = starfit(spec, data = dx[1:521], solver = "strategy", solver.control = solver.control, n=5)
			ssim = sample(as.numeric(dxv), 5000, replace=TRUE)
			sim  = starsim(fit, n.sim=5000, ssim = list(matrix(ssim, ncol=1)))
			specx = starspec(mean.model = list(states = 2, include.intercept = int[[i]], arOrder = arorder[[j]], 
							statevar = c("y","s")[2], s = dxv[500:521]),
					variance.model=list(dynamic=TRUE, model="sGARCH"))
			setfixed(specx)<-as.list(coef(fit))
			path = starpath(specx, n.sim=5000, prereturns = tail(dx[1:521], 6), 
					preresiduals = tail(residuals(fit),2),
					presigma = tail(sigma(fit), 2),
					custom.dist=list(name="sample", distfit=sim@simulation$residSim/sigma(sim)),
					ssim=list(matrix(ssim, ncol=1)))
			test[i,j] = as.integer(all.equal(fitted(sim), fitted(path)))
			print(test[i,j])
		}
	}
	zz <- file("test2c-3.txt", open="wt")
	sink(zz)
	print(test)
	sink(type="message")
	sink()
	close(zz)
	
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}