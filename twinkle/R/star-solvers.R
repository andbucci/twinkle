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

# need a more sophisticated parameter search method which looks at the only searching
# for the state variables
.scaleparsearch = function(fun, arglist)
{
	ipars = arglist$ipars
	arglist$scale = FALSE
	arglist$transform = TRUE
	arglist$returnType="llh"
	pars = ipars[which(ipars[,"Estimate"]==1),1]
	pars = runif(length(pars))
	Dup = Ddn = rep(0, length(pars))
	opt = optim(pars, fn = fun, arglist=arglist, method="BFGS", control=list(trace=0, maxit=20))
	pars = opt$par
	D0 = opt$value
	adup = pars +1
	addn = pars -1
	for(i in 1:length(pars)){
		xpars = pars
		xpars[i] = adup[i]
		Dup[i] = fun(xpars, arglist)
	}
	for(i in 1:length(pars)){
		xpars = pars
		xpars[i] = addn[i]
		Ddn[i] = fun(xpars, arglist)
	}
	parscale = apply(cbind(D0/Ddn, D0/Ddn), 1, "mean")
	return(1/parscale)
}

# 1. startpars uses upper and lower bounds so transform is FALSE.
#    Returned parameters should then be 'transformed' to the
#    untransformed domain (those with [,7]==1)
# 2. If parsearch=FALSE and the solver is unconstrained, the routine 
#    will transform those parameters [,7]==1 and belonging to the
#    estimation set [,'Estimate']==1.
.starparsearch = function(fun, arglist, solver, solver.control, cluster = NULL, solvertype="unconstrained")
{
	ipars = arglist$ipars
	if(is.null(solver.control$parsearch)) parsearch = TRUE else parsearch = as.logical(solver.control$parsearch)
	if(is.null(solver.control$searchtype)) searchtype = "fast" else searchtype = solver.control$searchtype
	if(is.null(solver.control$n.restarts)) n.restarts = 1 else n.restarts = as.integer(solver.control$n.restarts)
	if(!is.null(solver.control$rseed)) rseed = solver.control$rseed else rseed = NULL
	if(parsearch){
		if(searchtype=="likelihood"){
			if(is.null(solver.control$parsim)) parsim = 5000 else parsim = abs(as.integer(solver.control$parsim))
			arglist$transform = FALSE
			pars = ipars[which(ipars[,"Estimate"]==1),1]
			LB = ipars[which(ipars[,"Estimate"]==1),"LB"]
			UB = ipars[which(ipars[,"Estimate"]==1),"UB"]
			spars = startpars(pars  = pars, fun = fun, LB = LB, UB = UB, arglist = arglist, rseed = rseed,
					n.sim = parsim, bestN = n.restarts, cluster = cluster)
			if(n.restarts==1){
				spars = matrix(as.numeric(spars), nrow=1)
				spars = spars[,-length(spars),drop=FALSE]
			} else{
				spars = spars[,-length(spars),drop=FALSE]
			}
			if(solvertype=="unconstrained"){
				tr = ipars[which(ipars[,"Estimate"]==1),7]
				for(i in 1:ncol(spars)){
					if(tr[i]==1){
						spars[,i] = logtransform(spars[,i], lower = LB[i], upper = UB[i], inverse=TRUE)
						if(is.nan(spars[1,i]) | !is.finite(spars[1,i])) spars[,i] = logtransform(0.99*min(max(pars[i], LB[i]),UB[i]), lower = LB[i], upper = UB[i], inverse=TRUE)
					}
				}				
			}
		} else{
			pars = ipars[which(ipars[,"Estimate"]==1),1]
			if(!is.null(rseed)) set.seed(rseed)
			spars = matrix(runif(length(pars)*n.restarts), nrow = n.restarts, ncol = length(pars))
		}
	} else{
		pars = ipars[which(ipars[,"Estimate"]==1),1]
		tr = ipars[which(ipars[,"Estimate"]==1),7]
		spars = matrix(pars, nrow = n.restarts, ncol = length(pars), byrow=TRUE)
		if(solvertype=="unconstrained"){
			LB = ipars[which(ipars[,"Estimate"]==1),"LB"]
			UB = ipars[which(ipars[,"Estimate"]==1),"UB"]
			for(i in 1:ncol(spars)){
				if(tr[i]==1){
					spars[,i] = logtransform(pars[i], lower = LB[i], upper = UB[i], inverse=TRUE)
					if(is.nan(spars[1,i]) | !is.finite(spars[1,i])) spars[,i] = logtransform(0.99*min(max(pars[i], LB[i]),UB[i]), lower = LB[i], upper = UB[i], inverse=TRUE)
				}
			}
		} else{
			spars = matrix(ipars[which(ipars[,"Estimate"]==1),1], nrow = n.restarts, ncol = length(pars), byrow=TRUE)
		}
	}
	return(spars)
}

fun2dstar = function(pars, y, x)
{
	mu1 = mean(y[y<0])
	mu2 = mean(y[y>0])
	cx = pars[1]
	ax = pars[2:length(pars)]
	p = cx + ax %*% x
	sum( (y - (p*mu1 + (1-p)*mu2))^2 )
}


.starsolver = function(solver, pars, fun, Ifn, ILB, IUB, gr, hessian, parscale, control, arglist, cluster)
{
	if(is.null(control$usepars)){
		usepars=FALSE
	} else{
		usepars = control$usepars
		control$usepars = NULL
	}
	# ToDo: if it is one of the bounded solvers need to adjust LB and UB to be outside of the pars provided
	retval = switch(solver,
			optim = .optimsolver(pars = pars, fun = fun, gr = gr, control = control, arglist = arglist, cluster = cluster, usepars = usepars),
			msoptim = .msoptimsolver(pars = pars, fun = fun, gr = gr, control = control, arglist = arglist, cluster = cluster),
			nlminb = .nlminbsolver(pars = pars, fun = fun, gr = gr, control = control, arglist = arglist, cluster = cluster, usepars = usepars),
			solnp = .solnpsolver(pars = pars, fun = fun, Ifn = Ifn, ILB = ILB, IUB = IUB, control = control, arglist = arglist, usepars = usepars),
			cmaes = .cmaessolver(pars = pars, fun = fun, gr = gr, control = control, arglist = arglist, cluster = cluster, usepars = usepars),
			deoptim = .deoptimsolver(pars = pars, fun = fun, gr = gr, control = control, arglist = arglist, cluster = cluster, usepars = usepars))
	return(retval)
}

.solnpctrl = function(control){
	# parameters check is now case independent
	ans = list()
	params = unlist(control)
	if(is.null(params)) {
		ans$rho = 1
		ans$outer.iter = 50
		ans$inner.iter = 1800
		ans$delta = 1.0e-8
		ans$tol = 1.0e-8
		ans$trace = 0
	} else{
		npar = tolower(names(unlist(control)))
		names(params) = npar
		if(any(substr(npar, 1, 3) == "rho")) ans$rho = as.numeric(params["rho"]) else ans$rho = 1
		if(any(substr(npar, 1, 10) == "outer.iter")) ans$outer.iter = as.numeric(params["outer.iter"]) else ans$outer.iter = 50
		if(any(substr(npar, 1, 10) == "inner.iter")) ans$inner.iter = as.numeric(params["inner.iter"]) else ans$inner.iter = 1000
		if(any(substr(npar, 1, 5) == "delta")) ans$delta = as.numeric(params["delta"]) else ans$delta = 1.0e-8
		if(any(substr(npar, 1, 3) == "tol")) ans$tol = as.numeric(params["tol"]) else ans$tol = 1.0e-8
		if(any(substr(npar, 1, 5) == "trace")) ans$trace = as.numeric(params["trace"]) else ans$trace = 0
	}
	return(ans)
}
.solnpsolver = function(pars, fun, Ifn, ILB, IUB, control, LB, UB, arglist, usepars = FALSE, cluster = NULL){
	if(!usepars) pars = as.numeric( .starparsearch(fun, arglist, solver="optim", solver.control = control, cluster = cluster, solvertype="constrained") )
	if(!is.null(control$parsearch)) control$parsearch = NULL
	if(!is.null(control$searchtype)) control$searchtype = NULL
	
	if(!is.null(control$parsim)) control$parsim = NULL
	if(!is.null(control$rseed)) control$rseed = NULL
	lower = arglist$ipars[which(arglist$ipars[,"Estimate"]==1),"LB"]
	upper = arglist$ipars[which(arglist$ipars[,"Estimate"]==1),"UB"]
	arglist$transform = FALSE
	ans = try(solnp(pars, fun = fun, eqfun = NULL, eqB = NULL, ineqfun = Ifn, ineqLB = ILB, 
					ineqUB = IUB, LB = lower, UB = upper, control = control, arglist), silent = TRUE)
	if(inherits(ans, "try-error")){
		sol = list()
		sol$convergence = 1
		sol$message = "Error"
		sol$pars = rep(NA, length(pars))
		names(sol$pars) = names(pars)
		sol$lik = NA
		hess = NULL
	}
	else{
		arglist$returnType = "all"
		ipars = fun(ans$par, arglist)$ipars
		sol = ans
		sol$pars = ipars[which(ipars[,"Estimate"]==1),1]
		names(sol$pars) = names(pars)
		sol$par = NULL
		sol$convergence = ans$convergence
		hess = NULL
	}
	if(sol$convergence!=0) warning("\ntwinkle-->warning: no convergence...\n")
	return(list(sol = sol, hess = hess))
}


.optimsolver = function(pars, fun, gr = NULL, control = list(), arglist, cluster = NULL, usepars = FALSE){
	# optim will not check for default control parameters
	if(is.null(control$method)){
		solvertype="unconstrained"
	} else{
		if(control$method=="L-BFGS-B") solvertype="constrained" else solvertype="unconstrained"
	}
	if(!usepars) pars = as.numeric( .starparsearch(fun, arglist, solver="optim", solver.control = control, cluster = cluster,
						solvertype=solvertype) )
	if(!is.null(control$parsearch)) control$parsearch = NULL
	if(!is.null(control$searchtype)) control$searchtype = NULL
	if(!is.null(control$n.restarts)) control$n.restarts = NULL
	if(!is.null(control$parsim)) control$parsim = NULL
	if(!is.null(control$rseed)) control$rseed = NULL
	arglist$transform = TRUE
	if(is.null(control$method)){
		method = "BFGS"
		lower = -Inf
		upper = Inf
		#control$parscale = .scaleparsearch(fun, arglist)
	} else{
		method = control$method
		if(method=="L-BFGS-B"){
			lower = arglist$ipars[which(arglist$ipars[,"Estimate"]==1),"LB"]
			upper = arglist$ipars[which(arglist$ipars[,"Estimate"]==1),"UB"]
			arglist$transform = FALSE
		} else{
			#control$parscale = .scaleparsearch(fun, arglist)
			lower = -Inf
			upper = Inf
		}
	}
	control$method = NULL
	ans = optim(fn = fun, gr = NULL, par = pars, arglist = arglist, control = control, method = method,
			lower = lower, upper = upper)
	if (inherits(ans, "try-error")){
		sol = list()
		sol$convergence = 1
		sol$message = ans$message
		sol$pars = rep(NA, length(pars))
		names(sol$pars) = names(pars)
		sol$lik = NA
		hess = NULL
	}
	else{
		arglist$returnType = "all"
		ipars = fun(ans$par, arglist)$ipars
		sol = ans
		sol$pars = ipars[which(ipars[,"Estimate"]==1),1]
		names(sol$pars) = names(pars)
		sol$par = NULL
		sol$convergence = ans$convergence
		hess = NULL
	}
	if(sol$convergence!=0) warning("\ntwinkle-->warning: no convergence...\n")
	return(list(sol = sol, hess = hess))
}

.bfgssolver = function(pars, fun, gr = NULL, control = list(), arglist, cluster = NULL, usepars = FALSE, transform.result = FALSE){
	# optim will not check for default control parameters
	if(!usepars) pars = as.numeric( .starparsearch(fun, arglist, solver="optim", solver.control = control, cluster = cluster,
						solvertype = "unconstrained") )
	if(!is.null(control$parsearch)) control$parsearch = NULL
	if(!is.null(control$searchtype)) control$searchtype = NULL
	if(!is.null(control$n.restarts)) control$n.restarts = NULL
	if(!is.null(control$parsim)) control$parsim = NULL
	if(!is.null(control$rseed)) control$rseed = NULL
	arglist$transform = TRUE
	method = "BFGS"
	lower = -Inf
	upper = Inf
	control$method = NULL
	ans = optim(fn = fun, gr = NULL, par = pars, arglist = arglist, control = control, method = method, lower = lower, upper = upper)
	if (inherits(ans, "try-error")){
		sol = list()
		sol$convergence = 1
		sol$message = ans$message
		sol$pars = rep(NA, length(pars))
		names(sol$pars) = names(pars)
		sol$lik = NA
		hess = NULL
	}
	else{
		sol = ans
		if(transform.result){
			arglist$returnType = "all"
			ipars = fun(ans$par, arglist)$ipars
			sol = ans
			sol$pars = ipars[which(ipars[,"Estimate"]==1),1]
		} else{
			sol$pars = ans$par
		}
		names(sol$pars) = names(pars)
		sol$par = NULL
		sol$convergence = ans$convergence
		hess = NULL
	}
	if(sol$convergence!=0) warning("\ntwinkle-->warning: no convergence...\n")
	return(list(sol = sol, hess = hess))
}



.nlminbsolver = function(pars, fun, gr = NULL, control = list(), arglist, cluster = NULL, usepars = FALSE){
	# optim will not check for default control parameters
	if(!usepars) pars = as.numeric( .starparsearch(fun, arglist, solver="optim", solver.control = control, cluster = cluster,
						solvertype = "constrained") )
	if(!is.null(control$parsearch)) control$parsearch = NULL
	if(!is.null(control$searchtype)) control$searchtype = NULL
	if(!is.null(control$n.restarts)) control$n.restarts = NULL
	if(!is.null(control$parsim)) control$parsim = NULL
	if(!is.null(control$rseed)) control$rseed = NULL
	
	lower = arglist$ipars[which(arglist$ipars[,"Estimate"]==1),"LB"]
	upper = arglist$ipars[which(arglist$ipars[,"Estimate"]==1),"UB"]
	arglist$transform = FALSE
	ans = nlminb(start = pars, objective = fun, gradient = NULL, hessian = NULL, arglist = arglist, control = control,
			lower = lower, upper = upper)
	if (inherits(ans, "try-error")){
		sol = list()
		sol$convergence = 1
		sol$message = ans$message
		sol$pars = rep(NA, length(pars))
		names(sol$pars) = names(pars)
		sol$lik = NA
		hess = NULL
	}
	else{
		arglist$returnType = "all"
		ipars = fun(ans$par, arglist)$ipars
		sol = ans
		sol$pars = ipars[which(ipars[,"Estimate"]==1),1]
		names(sol$pars) = names(pars)
		sol$par = NULL
		sol$convergence = ans$convergence
		hess = NULL
	}
	if(sol$convergence!=0) warning("\ntwinkle-->warning: no convergence...\n")
	return(list(sol = sol, hess = hess))
}

.deoptimsolver = function(pars, fun, gr = NULL, control = list(), arglist, cluster = NULL, usepars = FALSE){
	if(is.null(control$NP)) control$NP = 10*length(pars)
	if(is.null(control$itermax)) control$itermax = 500
	if(is.null(control$F)) control$F = 0.5
	if(is.null(control$strategy)) control$strategy = 2
	control$n.restarts = control$NP
	if(is.null(control$parsearch)){
		control$parsearch = TRUE
		xpars = .starparsearch(fun, arglist, solver="optim", solver.control = control, cluster = cluster, solvertype="constrained")
	} else{
		if(control$parsearch){
			xpars = .starparsearch(fun, arglist, solver="optim", solver.control = control, cluster = cluster, solvertype="constrained")
		} else{
			xpars = matrix(pars, ncol = length(pars), nrow = control$NP, byrow=TRUE)
			xpars[-1,] = runif(length(pars)*(control$NP-1))
		}
	}
	control$initialpop = xpars
	if(!is.null(control$searchtype)) control$searchtype = NULL
	if(!is.null(control$parsearch)) control$parsearch = NULL
	if(!is.null(control$n.restarts)) control$n.restarts = NULL
	if(!is.null(control$parsim)) control$parsim = NULL
	if(!is.null(control$rseed)) control$rseed = NULL
	lower = arglist$ipars[which(arglist$ipars[,"Estimate"]==1),"LB"]
	upper = arglist$ipars[which(arglist$ipars[,"Estimate"]==1),"UB"]
	arglist$transform = FALSE
	ans = DEoptim(fn = fun, lower = lower, upper = upper, arglist = arglist, control = control)
	if (inherits(ans, "try-error")){
		sol = list()
		sol$convergence = 1
		sol$message = "no convergence"
		sol$pars = rep(NA, length(pars))
		names(sol$pars) = names(pars)
		sol$lik = NA
		hess = NULL
	}
	else{
		arglist$returnType = "all"
		ipars = fun(ans$optim$bestmem, arglist)$ipars
		sol = ans
		sol$pars = ipars[which(ipars[,"Estimate"]==1),1]
		names(sol$pars) = names(pars)
		sol$par = NULL
		sol$convergence = 0
		hess = NULL
	}
	if(sol$convergence!=0) warning("\ntwinkle-->warning: no convergence...\n")
	return(list(sol = sol, hess = hess))
}



.msoptimsolver = function(pars, fun, gr, control, arglist, cluster){
	if(is.null(control$method)){
		solvertype="unconstrained"
	} else{
		if(control$method=="L-BFGS-B") solvertype="constrained" else solvertype="unconstrained"
	}
	pars = .starparsearch(fun, arglist, solver="optim", solver.control = control, cluster = cluster, solvertype=solvertype)
	N = NROW(pars)
	xsol = vector(mode="list", length = N)
	if(!is.null(cluster)){
		clusterEvalQ(cluster, require(twinkle))
		clusterExport(cluster, c("fun", "control", "pars", "arglist"), envir = environment())
		xsol = parLapplyLB(cluster, 1:N, function(i){
					return( .optimsolver(pars[i,], fun, gr = NULL, control, arglist, cluster = NULL, usepars = TRUE))
				})
	} else{
		for(i in 1:N){
			xsol[[i]] = .optimsolver(pars[i,], fun, gr, control, arglist, cluster = NULL, usepars = TRUE)
		}
	}
	best = sapply(xsol, function(x) x$sol$value)
	best = which(best == min(best, na.rm=TRUE))[1]
	return(xsol[[best]])
}
.cmaes.ctrl = function(control){
	cmaes.control(options = control$options, CMA = control$CMA)
}

.cmaessolver = function(pars, fun, gr, control, arglist, cluster, usepars = FALSE)
{
	control1 = control
	control = .cmaes.ctrl(control)
	if(!usepars) pars = as.numeric( .starparsearch(fun, arglist, solver="cmaes", solver.control = control1, cluster = cluster,
						solvertype = "constrained") )
	arglist$transform = FALSE
	lower = arglist$ipars[which(arglist$ipars[,"Estimate"]==1),"LB"]
	upper = arglist$ipars[which(arglist$ipars[,"Estimate"]==1),"UB"]	
	ans = try(cmaes(pars, fun, lower = lower, upper = upper, insigma = 1, ctrl = control, arglist = arglist), silent=TRUE)
	if(inherits(ans, "try-error")){
		sol = list()
		sol$convergence = 1
		sol$message = ans$message
		sol$pars = rep(NA, length(pars))
		names(sol$par) = names(pars)
		sol$lik = NA
		hess = NULL
	} else{
			ansf = try(nlminb(start = ans$par, fun, lower = lower, upper = upper, arglist = arglist, 
						control = list(trace= ifelse(control$options$DispFinal, 1, 0), 
								eval.max = 2000, iter.max = 1000, step.min = 0.1)), silent = TRUE)
		if(inherits(ansf, "try-error") || ansf$convergence>0){
			# revert to cmaes solution
			arglist$returnType = "all"
			ipars = fun(ans$par, arglist)$ipars
			sol = ans
			sol$pars = ipars[which(ipars[,"Estimate"]==1),1]
			names(sol$pars) = names(pars)
			sol$par = NULL
			sol$convergence = 1
			hess = NULL
			sol$message = ans$stopflag
		}
		else{
			arglist$returnType = "all"
			ipars = fun(ansf$par, arglist)$ipars
			sol = ans
			sol$pars = ipars[which(ipars[,"Estimate"]==1),1]
			names(sol$pars) = names(pars)
			sol$par = NULL
			sol$convergence = 0
			hess = NULL
		}
	}
	if(sol$convergence!=0) warning("\ntwinkle-->warning: no convergence...\n")
	return(list(sol = sol, hess = hess))
}