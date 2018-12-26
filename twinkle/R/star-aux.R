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

#-----------------------------
# CONSTANTS used in package:
eps<-.Machine$double.eps
TinY = 1.0e-8
#-----------------------------
.makestarnames = function(modelinc)
{
	parnames = NULL
	if(modelinc[1]>0) parnames = c(parnames, paste('s',1,'.phi0',sep="",collapse="'"))
	if(modelinc[2]>0) parnames = c(parnames, strsplit(paste('s',1,'.phi',1:modelinc[2],sep="",collapse="'"),"'")[[1]])
	if(modelinc[3]>0) parnames = c(parnames, strsplit(paste('s',1,'.xi',1:modelinc[3],sep="",collapse="'"),"'")[[1]])
	if(modelinc[4]>0) parnames = c(parnames, strsplit(paste('s',1,'.psi',1:modelinc[4],sep="",collapse="'"),"'")[[1]])
	
	if(modelinc[5]>0) parnames = c(parnames, paste('s',2,'.phi0',sep="",collapse="'"))
	if(modelinc[6]>0) parnames = c(parnames, strsplit(paste('s',2,'.phi',1:modelinc[6],sep="",collapse="'"),"'")[[1]])
	if(modelinc[7]>0) parnames = c(parnames, strsplit(paste('s',2,'.xi',1:modelinc[7],sep="",collapse="'"),"'")[[1]])
	if(modelinc[8]>0) parnames = c(parnames, strsplit(paste('s',2,'.psi',1:modelinc[8],sep="",collapse="'"),"'")[[1]])
	
	if(modelinc[9]>0) parnames = c(parnames, paste('s',3,'.phi0',sep="",collapse="'"))
	if(modelinc[10]>0) parnames = c(parnames, strsplit(paste('s',3,'.phi',1:modelinc[10],sep="",collapse="'"),"'")[[1]])
	if(modelinc[11]>0) parnames = c(parnames, strsplit(paste('s',3,'.xi',1:modelinc[11],sep="",collapse="'"),"'")[[1]])
	if(modelinc[12]>0) parnames = c(parnames, strsplit(paste('s',3,'.psi',1:modelinc[12],sep="",collapse="'"),"'")[[1]])
	
	if(modelinc[13]>0) parnames = c(parnames, paste('s',4,'.phi0',sep="",collapse="'"))
	if(modelinc[14]>0) parnames = c(parnames, strsplit(paste('s',4,'.phi',1:modelinc[14],sep="",collapse="'"),"'")[[1]])
	if(modelinc[15]>0) parnames = c(parnames, strsplit(paste('s',4,'.xi',1:modelinc[15],sep="",collapse="'"),"'")[[1]])
	if(modelinc[16]>0) parnames = c(parnames, strsplit(paste('s',4,'.psi',1:modelinc[16],sep="",collapse="'"),"'")[[1]])
	
	if(modelinc[17]>0) parnames = c(parnames, strsplit(paste('psi',1:modelinc[17],sep="",collapse="'"),"'")[[1]])
	
	if(modelinc[18]>0) parnames = c(parnames, "s1.gamma")
	if(modelinc[19]>0) parnames = c(parnames, "s1.c")
	if(modelinc[20]>0) parnames = c(parnames, paste("s1.alpha",1:modelinc[20],sep=""))
	if(modelinc[21]>0) parnames = c(parnames, "s1.beta")
	
	if(modelinc[22]>0) parnames = c(parnames, "s2.gamma")
	if(modelinc[23]>0) parnames = c(parnames, "s2.c")
	if(modelinc[24]>0) parnames = c(parnames, paste("s2.alpha",1:modelinc[24],sep=""))
	if(modelinc[25]>0) parnames = c(parnames, "s2.beta")
	
	if(modelinc[26]>0) parnames = c(parnames, "s3.gamma")
	if(modelinc[27]>0) parnames = c(parnames, "s3.c")
	if(modelinc[28]>0) parnames = c(parnames, paste("s3.alpha",1:modelinc[28],sep=""))
	if(modelinc[29]>0) parnames = c(parnames, "s3.beta")
	
	if(modelinc[50]==4){
		if(modelinc[30]>0) parnames = c(parnames, "s1.sigma")
		if(modelinc[31]>0) parnames = c(parnames, "s2.sigma")
		if(modelinc[32]>0) parnames = c(parnames, "s3.sigma")
		if(modelinc[33]>0) parnames = c(parnames, "s4.sigma")
	} else{
		if(modelinc[30]>0) parnames = c(parnames, "sigma")
		if(modelinc[31]>0) parnames = c(parnames, "omega")
		if(modelinc[32]>0) parnames = c(parnames, paste("alpha",1:modelinc[32],sep=""))
		if(modelinc[33]>0) parnames = c(parnames, paste("beta",1:modelinc[33],sep=""))
	}
	if(modelinc[34]>0) parnames = c(parnames, paste("gamma",1:modelinc[34],sep=""))
	if(modelinc[35]>0) parnames = c(parnames, paste("eta1",1:modelinc[35],sep=""))
	if(modelinc[36]>0) parnames = c(parnames, paste("eta2",1:modelinc[36],sep=""))
	if(modelinc[37]>0) parnames = c(parnames, paste("delta",sep=""))
	if(modelinc[38]>0) parnames = c(parnames, paste("lambda",sep=""))
	if(modelinc[39]>0) parnames = c(parnames, paste("vxreg",1:modelinc[39],sep=""))
	if(modelinc[40]>0) parnames = c(parnames, paste("xi",sep=""))
	if(modelinc[41]>0) parnames = c(parnames, paste("skew",sep=""))
	if(modelinc[42]>0) parnames = c(parnames, paste("shape",sep=""))
	if(modelinc[43]>0) parnames = c(parnames, paste("ghlambda",sep=""))
	return(parnames)
}


logtransform = function(x, lower, upper, inverse = FALSE) 
{
	if(!inverse) {
		ans = lower + (upper - lower)/(1 + exp(-1 * x))
	}
	else {
		ans = -1 * log(-(upper - x)/(-x + lower))
	}
	return(ans)
}

minmaxtransform = function(x, inverse = FALSE, min.x = NULL, max.x = NULL)
{
	if(inverse){
		if(is.null(min.x)) stop("\ncannot inverse without min.x!")
		if(is.null(max.x)) stop("\ncannot inverse without max.x!")
		ans = (x * (max.x - min.x)) + min.x
	} else{
		min.x = min(x, na.rm = TRUE)
		max.x = max(x, na.rm = TRUE)
		ans = (x - min.x)/(max.x - min.x)
	}
	return(ans)
}

iqrtransform = function(x, type = 7, inverse = FALSE, median.x = FALSE, IQR.x = NULL)
{
	if(inverse){
		if(is.null(median.x)) stop("\ncannot inverse without median.x!")
		if(is.null(IQR.x)) stop("\ncannot inverse without IQR.x!")
		ans = (x * IQR.x) + median.x
	} else{
		ans = (x - median(x, na.rm = TRUE))/IQR(x, na.rm = TRUE, type = type)
	}
	return(ans)
}
check_fun = function(fun, n=1)
{
	if(length(formals(fun))!=1) stop("\nstarspec-->error: fun should only take one argument\n")
	if(n==1){
		test = do.call(fun, list(rnorm(100)))
		if(length(test)!=100) stop("\nstarspec-->error: fun should return same length as input\n")
	} else{
		test = do.call(fun, list(matrix(rnorm(100*n, 100, n))))
		if(!is.matrix(test)) stop(paste("\nstarspec-->error: fun should return a matrix when a matrix is input\n",sep=""))
		if(ncol(test)!=n) stop(paste("\nstarspec-->error: fun should return a matrix with ", n," columns (as inputed)\n",sep=""))
		if(nrow(test)!=100) stop(paste("\nstarspec-->error: fun should return a matrix with 25 rows (as inputed)\n",sep=""))
	}
	return(0)
}

.starstart = function(pars, arglist)
{
	.eps = .Machine$double.eps
	data = arglist$data
	model = arglist$model
	start.pars = model$start.pars
	start.names = names(start.pars)
	T = model$modeldata$T
	fixed.pars = model$fixed.pars
	fixed.names = names(fixed.pars)
	idx = model$pidx
	modelinc = model$modelinc
	mexdata = model$modeldata$mexdata[1:T, , drop=FALSE]
	if(!is.null(mexdata)) colnames(mexdata)  = NULL
	
	#----------------------------------
	# ARX dynamics
	#----------------------------------
	if(modelinc[1]>0){
		if(is.na(pars[idx["s1.phi0", 1]:idx["s1.phi0", 2], 5])) pars[idx["s1.phi0", 1]:idx["s1.phi0", 2], 5] = min(data)
		if(is.na(pars[idx["s1.phi0", 1]:idx["s1.phi0", 2], 6])) pars[idx["s1.phi0", 1]:idx["s1.phi0", 2], 6] = max(data)
		if(!is.null(start.pars$s1.phi0)) pars[idx["s1.phi0", 1]:idx["s1.phi0", 2], 1] = start.pars$s1.phi0[1] else pars[idx["s1.phi0", 1]:idx["s1.phi0", 2], 1] = mean(data)
		if(any(substr(fixed.names, 1, 7)=="s1.phi0")){
			pars[idx["s1.phi0", 1]:idx["s1.phi0", 2], 1] = as.numeric(fixed.pars$s1.phi0)
			pars[idx["s1.phi0", 1]:idx["s1.phi0", 2], 5] = fixed.pars$s1.phi0
			pars[idx["s1.phi0", 1]:idx["s1.phi0", 2], 6] = fixed.pars$s1.phi0
		}
	}
	if(modelinc[5]>0){
		if(is.na(pars[idx["s2.phi0", 1]:idx["s2.phi0", 2], 5])) pars[idx["s2.phi0", 1]:idx["s2.phi0", 2], 5] = min(data)
		if(is.na(pars[idx["s2.phi0", 1]:idx["s2.phi0", 2], 6])) pars[idx["s2.phi0", 1]:idx["s2.phi0", 2], 6] = max(data)
		if(!is.null(start.pars$s2.phi0)) pars[idx["s2.phi0", 1]:idx["s2.phi0", 2], 1] = start.pars$s2.phi0[1] else pars[idx["s2.phi0", 1]:idx["s2.phi0", 2], 1] = mean(data)
		if(any(substr(fixed.names, 1, 7)=="s2.phi0")){
			pars[idx["s2.phi0", 1]:idx["s2.phi0", 2], 1] = as.numeric(fixed.pars$s2.phi0)
			pars[idx["s2.phi0", 1]:idx["s2.phi0", 2], 5] = fixed.pars$s2.phi0
			pars[idx["s2.phi0", 1]:idx["s2.phi0", 2], 6] = fixed.pars$s2.phi0
		}
	}
	if(modelinc[9]>0){
		if(is.na(pars[idx["s3.phi0", 1]:idx["s3.phi0", 2], 5])) pars[idx["s3.phi0", 1]:idx["s3.phi0", 2], 5] = min(data)
		if(is.na(pars[idx["s3.phi0", 1]:idx["s3.phi0", 2], 6])) pars[idx["s3.phi0", 1]:idx["s3.phi0", 2], 6] = max(data)
		if(!is.null(start.pars$s3.phi0)) pars[idx["s3.phi0", 1]:idx["s3.phi0", 2], 1] = start.pars$s3.phi0[1] else pars[idx["s3.phi0", 1]:idx["s3.phi0", 2], 1] = mean(data)
		if(any(substr(fixed.names, 1, 7)=="s3.phi0")){
			pars[idx["s3.phi0", 1]:idx["s3.phi0", 2], 1] = as.numeric(fixed.pars$s3.phi0)
			pars[idx["s3.phi0", 1]:idx["s3.phi0", 2], 5] = fixed.pars$s3.phi0
			pars[idx["s3.phi0", 1]:idx["s3.phi0", 2], 6] = fixed.pars$s3.phi0
		}
	}
	if(modelinc[13]>0){
		if(is.na(pars[idx["s4.phi0", 1]:idx["s4.phi0", 2], 5])) pars[idx["s4.phi0", 1]:idx["s4.phi0", 2], 5] = min(data)
		if(is.na(pars[idx["s4.phi0", 1]:idx["s4.phi0", 2], 6])) pars[idx["s4.phi0", 1]:idx["s4.phi0", 2], 6] = max(data)
		if(!is.null(start.pars$s4.phi0)) pars[idx["s4.phi0", 1]:idx["s4.phi0", 2], 1] = start.pars$s4.phi0[1] else pars[idx["s4.phi0", 1]:idx["s4.phi0", 2], 1] = mean(data)
		if(any(substr(fixed.names, 1, 7)=="s4.phi0")){
			pars[idx["s4.phi0", 1]:idx["s4.phi0", 2], 1] = as.numeric(fixed.pars$s4.phi0)
			pars[idx["s4.phi0", 1]:idx["s4.phi0", 2], 5] = fixed.pars$s4.phi0
			pars[idx["s4.phi0", 1]:idx["s4.phi0", 2], 6] = fixed.pars$s4.phi0
		}
	}
	
	if(modelinc[2]>0){
		gpnames = paste("s1.phi",1:modelinc[2],sep="")
		pxd = which(is.na(pars[idx["s1.phi", 1]:idx["s1.phi", 2], 5]))
		if(length(pxd)>0) pars[(idx["s1.phi", 1]:idx["s1.phi", 2])[pxd], 5] = -2
		pxd = which(is.na(pars[idx["s1.phi", 1]:idx["s1.phi", 2], 6]))
		if(length(pxd)>0) pars[(idx["s1.phi", 1]:idx["s1.phi", 2])[pxd], 6] =  2-TinY
		pars[idx["s1.phi", 1]:idx["s1.phi", 2], 1] = 0
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
	}
	
	
	if(modelinc[6]>0){
		gpnames = paste("s2.phi",1:modelinc[6],sep="")
		pxd = which(is.na(pars[idx["s2.phi", 1]:idx["s2.phi", 2], 5]))
		if(length(pxd)>0) pars[(idx["s2.phi", 1]:idx["s2.phi", 2])[pxd], 5] = -2
		pxd = which(is.na(pars[idx["s2.phi", 1]:idx["s2.phi", 2], 6]))
		if(length(pxd)>0) pars[(idx["s2.phi", 1]:idx["s2.phi", 2])[pxd], 6] =  2-TinY
		pars[idx["s2.phi", 1]:idx["s2.phi", 2], 1] = 0
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
	}
	
	if(modelinc[10]>0){
		gpnames = paste("s3.phi",1:modelinc[10],sep="")
		pxd = which(is.na(pars[idx["s3.phi", 1]:idx["s3.phi", 2], 5]))
		if(length(pxd)>0) pars[(idx["s3.phi", 1]:idx["s3.phi", 2])[pxd], 5] = -2
		pxd = which(is.na(pars[idx["s3.phi", 1]:idx["s3.phi", 2], 6]))
		if(length(pxd)>0) pars[(idx["s3.phi", 1]:idx["s3.phi", 2])[pxd], 6] =  2-TinY
		pars[idx["s3.phi", 1]:idx["s3.phi", 2], 1] = 0
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
	}
	
	if(modelinc[14]>0){
		gpnames = paste("s4.phi",1:modelinc[14],sep="")
		pxd = which(is.na(pars[idx["s4.phi", 1]:idx["s4.phi", 2], 5]))
		if(length(pxd)>0) pars[(idx["s4.phi", 1]:idx["s4.phi", 2])[pxd], 5] = -2
		pxd = which(is.na(pars[idx["s4.phi", 1]:idx["s4.phi", 2], 6]))
		if(length(pxd)>0) pars[(idx["s4.phi", 1]:idx["s4.phi", 2])[pxd], 6] =  2-TinY
		pars[idx["s4.phi", 1]:idx["s4.phi", 2], 1] = 0
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
	}
	
	if(modelinc[3]>0){
		gpnames = paste("s1.xi",1:modelinc[3],sep="")
		pxd = which(is.na(pars[idx["s1.xi", 1]:idx["s1.xi", 2], 5]))
		if(length(pxd)>0) pars[(idx["s1.xi", 1]:idx["s1.xi", 2])[pxd], 5] = -2
		pxd = which(is.na(pars[idx["s1.xi", 1]:idx["s1.xi", 2], 6]))
		if(length(pxd)>0) pars[(idx["s1.xi", 1]:idx["s1.xi", 2])[pxd], 6] =  2-TinY
		pars[idx["s1.xi", 1]:idx["s1.xi", 2], 1] = 0
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
	}
	
	if(modelinc[7]>0){
		gpnames = paste("s2.xi",1:modelinc[7],sep="")
		pxd = which(is.na(pars[idx["s2.xi", 1]:idx["s2.xi", 2], 5]))
		if(length(pxd)>0) pars[(idx["s2.xi", 1]:idx["s2.xi", 2])[pxd], 5] = -2
		pxd = which(is.na(pars[idx["s2.xi", 1]:idx["s2.xi", 2], 6]))
		if(length(pxd)>0) pars[(idx["s2.xi", 1]:idx["s2.xi", 2])[pxd], 6] =  2-TinY
		pars[idx["s2.xi", 1]:idx["s2.xi", 2], 1] = 0
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
	}
	
	if(modelinc[11]>0){
		gpnames = paste("s3.xi",1:modelinc[11],sep="")
		pxd = which(is.na(pars[idx["s3.xi", 1]:idx["s3.xi", 2], 5]))
		if(length(pxd)>0) pars[(idx["s3.xi", 1]:idx["s3.xi", 2])[pxd], 5] = -2
		pxd = which(is.na(pars[idx["s3.xi", 1]:idx["s3.xi", 2], 6]))
		if(length(pxd)>0) pars[(idx["s3.xi", 1]:idx["s3.xi", 2])[pxd], 6] =  2-TinY
		pars[idx["s3.xi", 1]:idx["s3.xi", 2], 1] = 0
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
	}
	if(modelinc[15]>0){
		gpnames = paste("s4.xi",1:modelinc[15],sep="")
		pxd = which(is.na(pars[idx["s4.xi", 1]:idx["s4.xi", 2], 5]))
		if(length(pxd)>0) pars[(idx["s4.xi", 1]:idx["s4.xi", 2])[pxd], 5] = -2
		pxd = which(is.na(pars[idx["s4.xi", 1]:idx["s4.xi", 2], 6]))
		if(length(pxd)>0) pars[(idx["s4.xi", 1]:idx["s4.xi", 2])[pxd], 6] =  2-TinY
		pars[idx["s4.xi", 1]:idx["s4.xi", 2], 1] = 0
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
	}
	#----------------------------------
	# MA dynamics
	#----------------------------------
	if(modelinc[4]>0){
		gpnames = paste("s1.psi",1:modelinc[4],sep="")
		pxd = which(is.na(pars[idx["s1.psi", 1]:idx["s1.psi", 2], 5]))
		if(length(pxd)>0) pars[(idx["s1.psi", 1]:idx["s1.psi", 2])[pxd], 5] = -2
		pxd = which(is.na(pars[idx["s1.psi", 1]:idx["s1.psi", 2], 6]))
		if(length(pxd)>0) pars[(idx["s1.psi", 1]:idx["s1.psi", 2])[pxd], 6] =  2-TinY
		pars[idx["s1.psi", 1]:idx["s1.psi", 2], 1] = 0
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
	}
	if(modelinc[8]>0){
		gpnames = paste("s2.psi",1:modelinc[8],sep="")
		pxd = which(is.na(pars[idx["s2.psi", 1]:idx["s2.psi", 2], 5]))
		if(length(pxd)>0) pars[(idx["s2.psi", 1]:idx["s2.psi", 2])[pxd], 5] = -2
		pxd = which(is.na(pars[idx["s2.psi", 1]:idx["s2.psi", 2], 6]))
		if(length(pxd)>0) pars[(idx["s2.psi", 1]:idx["s2.psi", 2])[pxd], 6] =  2-TinY
		pars[idx["s2.psi", 1]:idx["s2.psi", 2], 1] = 0
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
	}
	if(modelinc[12]>0){
		gpnames = paste("s3.psi",1:modelinc[12],sep="")
		pxd = which(is.na(pars[idx["s3.psi", 1]:idx["s3.psi", 2], 5]))
		if(length(pxd)>0) pars[(idx["s3.psi", 1]:idx["s3.psi", 2])[pxd], 5] = -2
		pxd = which(is.na(pars[idx["s3.psi", 1]:idx["s3.psi", 2], 6]))
		if(length(pxd)>0) pars[(idx["s3.psi", 1]:idx["s3.psi", 2])[pxd], 6] =  2-TinY
		pars[idx["s3.psi", 1]:idx["s3.psi", 2], 1] = 0
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
	}
	if(modelinc[16]>0){
		gpnames = paste("s4.psi",1:modelinc[16],sep="")
		pxd = which(is.na(pars[idx["s4.psi", 1]:idx["s4.psi", 2], 5]))
		if(length(pxd)>0) pars[(idx["s4.psi", 1]:idx["s4.psi", 2])[pxd], 5] = -2
		pxd = which(is.na(pars[idx["s4.psi", 1]:idx["s4.psi", 2], 6]))
		if(length(pxd)>0) pars[(idx["s4.psi", 1]:idx["s4.psi", 2])[pxd], 6] =  2-TinY
		pars[idx["s4.psi", 1]:idx["s4.psi", 2], 1] = 0
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
	}
	
	if(modelinc[17]>0){
		gpnames = paste("psi",1:modelinc[17],sep="")
		pxd = which(is.na(pars[idx["psi", 1]:idx["psi", 2], 5]))
		if(length(pxd)>0) pars[(idx["psi", 1]:idx["psi", 2])[pxd], 5] = -2
		pxd = which(is.na(pars[idx["psi", 1]:idx["psi", 2], 6]))
		if(length(pxd)>0) pars[(idx["psi", 1]:idx["psi", 2])[pxd], 6] =  2-TinY
		pars[idx["psi", 1]:idx["psi", 2], 1] = 0
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
	}
	#----------------------------------
	# star probability dynamics
	#----------------------------------
	if(modelinc[18]>0){
		if(is.na(pars[idx["s1.gamma", 1]:idx["s1.gamma", 2], 5])) pars[idx["s1.gamma", 1]:idx["s1.gamma", 2], 5] = 1e-12
		if(is.na(pars[idx["s1.gamma", 1]:idx["s1.gamma", 2], 6])) pars[idx["s1.gamma", 1]:idx["s1.gamma", 2], 6] = 1300
		if(!is.null(start.pars$s1.gamma)) pars[idx["s1.gamma", 1]:idx["s1.gamma", 2], 1] = start.pars$s1.gamma[1] else pars[idx["s1.gamma", 1]:idx["s1.gamma", 2], 1] = 0
		if(any(substr(fixed.names, 1, 8)=="s1.gamma")){
			pars[idx["s1.gamma", 1]:idx["s1.gamma", 2], 1] = as.numeric(fixed.pars$s1.gamma)
			pars[idx["s1.gamma", 1]:idx["s1.gamma", 2], 5] = fixed.pars$s1.gamma
			pars[idx["s1.gamma", 1]:idx["s1.gamma", 2], 6] = fixed.pars$s1.gamma
		}
	}
	if(modelinc[19]>0){
		if(is.na(pars[idx["s1.c", 1]:idx["s1.c", 2], 5])) pars[idx["s1.c", 1]:idx["s1.c", 2], 5] = -1300
		if(is.na(pars[idx["s1.c", 1]:idx["s1.c", 2], 6])) pars[idx["s1.c", 1]:idx["s1.c", 2], 6] =  1300
		if(!is.null(start.pars$s1.c)) pars[idx["s1.c", 1]:idx["s1.c", 2], 1] = start.pars$s1.c[1] else pars[idx["s1.c", 1]:idx["s1.c", 2], 1] = 0
		if(any(substr(fixed.names, 1, 4)=="s1.c")){
			pars[idx["s1.c", 1]:idx["s1.c", 2], 1] = as.numeric(fixed.pars$s1.c)
			pars[idx["s1.c", 1]:idx["s1.c", 2], 5] = fixed.pars$s1.c
			pars[idx["s1.c", 1]:idx["s1.c", 2], 6] = fixed.pars$s1.c
		}
	}
	
	if(modelinc[22]>0){
		if(is.na(pars[idx["s2.gamma", 1]:idx["s2.gamma", 2], 5])) pars[idx["s2.gamma", 1]:idx["s2.gamma", 2], 5] = -1300
		if(is.na(pars[idx["s2.gamma", 1]:idx["s2.gamma", 2], 6])) pars[idx["s2.gamma", 1]:idx["s2.gamma", 2], 6] =  1300
		if(!is.null(start.pars$s2.gamma)) pars[idx["s2.gamma", 1]:idx["s2.gamma", 2], 1] = start.pars$s2.gamma[1] else pars[idx["s2.gamma", 1]:idx["s2.gamma", 2], 1] = 0
		if(any(substr(fixed.names, 1, 8)=="s2.gamma")){
			pars[idx["s2.gamma", 1]:idx["s2.gamma", 2], 1] = as.numeric(fixed.pars$s2.gamma)
			pars[idx["s2.gamma", 1]:idx["s2.gamma", 2], 5] = fixed.pars$s2.gamma
			pars[idx["s2.gamma", 1]:idx["s2.gamma", 2], 6] = fixed.pars$s2.gamma
		}
	}
	
	if(modelinc[23]>0){
		if(is.na(pars[idx["s2.c", 1]:idx["s2.c", 2], 5])) pars[idx["s2.c", 1]:idx["s2.c", 2], 5] = -1300
		if(is.na(pars[idx["s2.c", 1]:idx["s2.c", 2], 6])) pars[idx["s2.c", 1]:idx["s2.c", 2], 6] =  1300
		if(!is.null(start.pars$s2.c)) pars[idx["s2.c", 1]:idx["s2.c", 2], 1] = start.pars$s2.c[1] else pars[idx["s2.c", 1]:idx["s2.c", 2], 1] = 0
		if(any(substr(fixed.names, 1, 4)=="s2.c")){
			pars[idx["s2.c", 1]:idx["s2.c", 2], 1] = as.numeric(fixed.pars$s2.c)
			pars[idx["s2.c", 1]:idx["s2.c", 2], 5] = fixed.pars$s2.c
			pars[idx["s2.c", 1]:idx["s2.c", 2], 6] = fixed.pars$s2.c
		}
	}
	
	if(modelinc[26]>0){
		if(is.na(pars[idx["s3.gamma", 1]:idx["s3.gamma", 2], 5])) pars[idx["s3.gamma", 1]:idx["s3.gamma", 2], 5] = -1300
		if(is.na(pars[idx["s3.gamma", 1]:idx["s3.gamma", 2], 6])) pars[idx["s3.gamma", 1]:idx["s3.gamma", 2], 6] =  1300
		if(!is.null(start.pars$s3.gamma)) pars[idx["s3.gamma", 1]:idx["s3.gamma", 2], 1] = start.pars$s3.gamma[1] else pars[idx["s3.gamma", 1]:idx["s3.gamma", 2], 1] = 0
		if(any(substr(fixed.names, 1, 8)=="s3.gamma")){
			pars[idx["s3.gamma", 1]:idx["s3.gamma", 2], 1] = as.numeric(fixed.pars$s3.gamma)
			pars[idx["s3.gamma", 1]:idx["s3.gamma", 2], 5] = fixed.pars$s3.gamma
			pars[idx["s3.gamma", 1]:idx["s3.gamma", 2], 6] = fixed.pars$s3.gamma
		}
	}
	
	if(modelinc[27]>0){
		if(is.na(pars[idx["s3.c", 1]:idx["s3.c", 2], 5])) pars[idx["s3.c", 1]:idx["s3.c", 2], 5] = -1300
		if(is.na(pars[idx["s3.c", 1]:idx["s3.c", 2], 6])) pars[idx["s3.c", 1]:idx["s3.c", 2], 6] =  1300
		if(!is.null(start.pars$s3.c)) pars[idx["s3.c", 1]:idx["s3.c", 2], 1] = start.pars$s3.c[1] else pars[idx["s3.c", 1]:idx["s3.c", 2], 1] = 0
		if(any(substr(fixed.names, 1, 4)=="s3.c")){
			pars[idx["s3.c", 1]:idx["s3.c", 2], 1] = as.numeric(fixed.pars$s3.c)
			pars[idx["s3.c", 1]:idx["s3.c", 2], 5] = fixed.pars$s3.c
			pars[idx["s3.c", 1]:idx["s3.c", 2], 6] = fixed.pars$s3.c
		}
	}
	
	
	if(modelinc[20]>0){
		gpnames = paste("s1.alpha",1:modelinc[20],sep="")
		pxd = which(is.na(pars[idx["s1.alpha", 1]:idx["s1.alpha", 2], 5]))
		if(length(pxd)>0) pars[(idx["s1.alpha", 1]:idx["s1.alpha", 2])[pxd], 5] = -1300
		pxd = which(is.na(pars[idx["s1.alpha", 1]:idx["s1.alpha", 2], 6]))
		if(length(pxd)>0) pars[(idx["s1.alpha", 1]:idx["s1.alpha", 2])[pxd], 6] =  1300
		pars[idx["s1.alpha", 1]:idx["s1.alpha", 2], 1] = 1
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
		# identification restriction imposition
		pars[idx["s1.alpha", 1], 1:2] = 1
		pars[idx["s1.alpha", 1], 4] = 0
		pars[idx["s1.alpha", 1], 5] = 1
		pars[idx["s1.alpha", 1], 6] = 1
	}
	
	if(modelinc[24]>0){
		gpnames = paste("s2.alpha",1:modelinc[24],sep="")
		pxd = which(is.na(pars[idx["s2.alpha", 1]:idx["s2.alpha", 2], 5]))
		if(length(pxd)>0) pars[(idx["s2.alpha", 1]:idx["s2.alpha", 2])[pxd], 5] = -1300
		pxd = which(is.na(pars[idx["s2.alpha", 1]:idx["s2.alpha", 2], 6]))
		if(length(pxd)>0) pars[(idx["s2.alpha", 1]:idx["s2.alpha", 2])[pxd], 6] =  1300
		pars[idx["s2.alpha", 1]:idx["s2.alpha", 2], 1] = 1
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
		# identification restriction
		pars[idx["s2.alpha", 1], 1:2] = 1
		pars[idx["s2.alpha", 1], 5] = 1
		pars[idx["s2.alpha", 1], 6] = 1
		pars[idx["s2.alpha", 1], 4] = 0
		
	}
	
	if(modelinc[28]>0){
		gpnames = paste("s3.alpha",1:modelinc[28],sep="")
		pxd = which(is.na(pars[idx["s3.alpha", 1]:idx["s3.alpha", 2], 5]))
		if(length(pxd)>0) pars[(idx["s3.alpha", 1]:idx["s3.alpha", 2])[pxd], 5] = -1300
		pxd = which(is.na(pars[idx["s3.alpha", 1]:idx["s3.alpha", 2], 6]))
		if(length(pxd)>0) pars[(idx["s3.alpha", 1]:idx["s3.alpha", 2])[pxd], 6] =  1300
		pars[idx["s3.alpha", 1]:idx["s3.alpha", 2], 1] = 1
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
		# identification restriction
		pars[idx["s3.alpha", 1], 1:2] = 1
		pars[idx["s3.alpha", 1], 5] = 1
		pars[idx["s3.alpha", 1], 6] = 1
		pars[idx["s3.alpha", 1], 4] = 0
		
	}
	
	if(modelinc[21]>0){
		if(is.na(pars[idx["s1.beta", 1]:idx["s1.beta", 2], 5])) pars[idx["s1.beta", 1]:idx["s1.beta", 2], 5] = -0.999
		if(is.na(pars[idx["s1.beta", 1]:idx["s1.beta", 2], 6])) pars[idx["s1.beta", 1]:idx["s1.beta", 2], 6] =  0.999
		if(!is.null(start.pars$s1.beta)) pars[idx["s1.beta", 1]:idx["s1.beta", 2], 1] = start.pars$s1.beta[1] else pars[idx["s1.beta", 1]:idx["s1.beta", 2], 1] = 0.5
		if(any(substr(fixed.names, 1, 7)=="s1.beta")){
			pars[idx["s1.beta", 1]:idx["s1.beta", 2], 1] = as.numeric(fixed.pars$s1.beta)
			pars[idx["s1.beta", 1]:idx["s1.beta", 2], 5] = fixed.pars$s1.beta
			pars[idx["s1.beta", 1]:idx["s1.beta", 2], 6] = fixed.pars$s1.beta
		}
	}
	
	if(modelinc[25]>0){
		if(is.na(pars[idx["s2.beta", 1]:idx["s2.beta", 2], 5])) pars[idx["s2.beta", 1]:idx["s2.beta", 2], 5] = -0.999
		if(is.na(pars[idx["s2.beta", 1]:idx["s2.beta", 2], 6])) pars[idx["s2.beta", 1]:idx["s2.beta", 2], 6] =  0.999
		if(!is.null(start.pars$s2.beta)) pars[idx["s2.beta", 1]:idx["s2.beta", 2], 1] = start.pars$s2.beta[1] else pars[idx["s2.beta", 1]:idx["s2.beta", 2], 1] = 0.5
		if(any(substr(fixed.names, 1, 7)=="s2.beta")){
			pars[idx["s2.beta", 1]:idx["s2.beta", 2], 1] = as.numeric(fixed.pars$s2.beta)
			pars[idx["s2.beta", 1]:idx["s2.beta", 2], 5] = fixed.pars$s2.beta
			pars[idx["s2.beta", 1]:idx["s2.beta", 2], 6] = fixed.pars$s2.beta
		}
	}
	
	if(modelinc[29]>0){
		if(is.na(pars[idx["s3.beta", 1]:idx["s3.beta", 2], 5])) pars[idx["s3.beta", 1]:idx["s3.beta", 2], 5] = -0.999
		if(is.na(pars[idx["s3.beta", 1]:idx["s3.beta", 2], 6])) pars[idx["s3.beta", 1]:idx["s3.beta", 2], 6] =  0.999
		if(!is.null(start.pars$s3.beta)) pars[idx["s3.beta", 1]:idx["s3.beta", 2], 1] = start.pars$s3.beta[1] else pars[idx["s3.beta", 1]:idx["s3.beta", 2], 1] = 0.5
		if(any(substr(fixed.names, 1, 7)=="s3.beta")){
			pars[idx["s3.beta", 1]:idx["s3.beta", 2], 1] = as.numeric(fixed.pars$s3.beta)
			pars[idx["s3.beta", 1]:idx["s3.beta", 2], 5] = fixed.pars$s3.beta
			pars[idx["s3.beta", 1]:idx["s3.beta", 2], 6] = fixed.pars$s3.beta
		}
	}
	
	
	#----------------------------------
	# sigma dynamics
	#----------------------------------
	
	if(modelinc[50]==0){
		.sd = sd(data)
		pars[idx["sigma", 1]:idx["sigma", 2], 5] = 1e-12
		pars[idx["sigma", 1]:idx["sigma", 2], 6] = 100*.sd
		if(is.null(start.pars$sigma)) pars[idx["sigma", 1]:idx["sigma", 2], 1] = .sd else pars[idx["sigma", 1]:idx["sigma", 2], 1] = abs(start.pars$sigma[1])
		if(any(substr(fixed.names, 1, 5) == "sigma")){
			pars[idx["sigma", 1]:idx["sigma", 2], 1] = abs(as.numeric(fixed.pars$sigma))
			pars[idx["sigma", 1]:idx["sigma", 2], 5] = abs(fixed.pars$sigma)
			pars[idx["sigma", 1]:idx["sigma", 2], 6] = abs(fixed.pars$sigma)
		}
	} else{
		garchmodel = model$modeldesc$vmodel
		pars = switch(garchmodel, 
				sGARCH = .aparchstart(pars, arglist),
				eGARCH = .egarchstart(pars, arglist),
				gjrGARCH = .aparchstart(pars, arglist),
				apARCH = .aparchstart(pars, arglist),
				mixture = .mixturestart(pars, arglist))
	}
	#----------------------------------
	# distribution
	#----------------------------------
	dbounds = .DistributionBounds(distribution = model$modeldesc$distribution)
	if(modelinc[41]>0){
		pars[idx["skew", 1]:idx["skew", 2], 5] = dbounds$skew.LB
		pars[idx["skew", 1]:idx["skew", 2], 6] = dbounds$skew.UB
		if(is.null(start.pars$skew)) pars[idx["skew", 1]:idx["skew", 2], 1] = dbounds$skew else pars[idx["skew", 1]:idx["skew", 2], 1] = start.pars$skew[1]
		if(any(substr(fixed.names, 1, 4) == "skew")){
			pars[idx["skew", 1]:idx["skew", 2], 1] = as.numeric(fixed.pars$skew)
			pars[idx["skew", 1]:idx["skew", 2], 5] = as.numeric(fixed.pars$skew)
			pars[idx["skew", 1]:idx["skew", 2], 6] = as.numeric(fixed.pars$skew)
		}
	}
	if(modelinc[42]>0){
		pars[idx["shape", 1]:idx["shape", 2], 5] = dbounds$shape.LB
		pars[idx["shape", 1]:idx["shape", 2], 6] = dbounds$shape.UB
		if(is.null(start.pars$shape)) pars[idx["shape", 1]:idx["shape", 2], 1] = dbounds$shape else pars[idx["shape", 1]:idx["shape", 2], 1] = start.pars$shape[1]
		if(any(substr(fixed.names, 1, 5) == "shape")){
			pars[idx["shape", 1]:idx["shape", 2], 1] = as.numeric(fixed.pars$shape)
			pars[idx["shape", 1]:idx["shape", 2], 5] = as.numeric(fixed.pars$shape)
			pars[idx["shape", 1]:idx["shape", 2], 6] = as.numeric(fixed.pars$shape)
		}
	}
	if(modelinc[43]>0){
		pars[idx["ghlambda", 1]:idx["ghlambda", 2], 5] = dbounds$ghlambda.LB
		pars[idx["ghlambda", 1]:idx["ghlambda", 2], 6] = dbounds$ghlambda.UB
		if(is.null(start.pars$ghlambda)) pars[idx["ghlambda", 1]:idx["ghlambda", 2], 1] = dbounds$ghlambda else pars[idx["ghlambda", 1]:idx["ghlambda", 2], 1] = start.pars$ghlambda[1]
		if(any(substr(fixed.names, 1, 8) == "ghlambda")){
			pars[idx["ghlambda", 1]:idx["ghlambda", 2], 1] = as.numeric(fixed.pars$ghlambda)
			pars[idx["ghlambda", 1]:idx["ghlambda", 2], 5] = as.numeric(fixed.pars$ghlambda)
			pars[idx["ghlambda", 1]:idx["ghlambda", 2], 6] = as.numeric(fixed.pars$ghlambda)
		}
	}
	return( pars )
}

.mixturestart = function(pars, arglist)
{
	data = arglist$data
	model = arglist$model
	modelinc = model$modelinc
	start.pars = model$start.pars
	fixed.pars = model$fixed.pars
	idx = model$pidx
	fixed.names = names(fixed.pars)
	start.names = names(start.pars)
	
	if(modelinc[30]>0){
		.sd = sd(data)
		pars[idx["s1.sigma", 1]:idx["s1.sigma", 2], 5] = 1e-12
		pars[idx["s1.sigma", 1]:idx["s1.sigma", 2], 6] = 10*.sd
		if(is.null(start.pars$s1.sigma)) pars[idx["s1.sigma", 1]:idx["s1.sigma", 2], 1] = .sd else pars[idx["s1.sigma", 1]:idx["s1.sigma", 2], 1] = abs(start.pars$s1.sigma[1])
		if(any(substr(fixed.names, 1, 8) == "s1.sigma")){
			pars[idx["s1.sigma", 1]:idx["s1.sigma", 2], 1] = abs(as.numeric(fixed.pars$s1.sigma))
			pars[idx["s1.sigma", 1]:idx["s1.sigma", 2], 5] = abs(fixed.pars$s1.sigma)
			pars[idx["s1.sigma", 1]:idx["s1.sigma", 2], 6] = abs(fixed.pars$s1.sigma)
		}
	}
	if(modelinc[31]>0){
		.sd = sd(data)
		pars[idx["s2.sigma", 1]:idx["s2.sigma", 2], 5] = 1e-12
		pars[idx["s2.sigma", 1]:idx["s2.sigma", 2], 6] = 10*.sd
		if(is.null(start.pars$s2.sigma)) pars[idx["s2.sigma", 1]:idx["s2.sigma", 2], 1] = .sd else pars[idx["s2.sigma", 1]:idx["s2.sigma", 2], 1] = abs(start.pars$s2.sigma[1])
		if(any(substr(fixed.names, 1, 8) == "s2.sigma")){
			pars[idx["s2.sigma", 1]:idx["s2.sigma", 2], 1] = abs(as.numeric(fixed.pars$s2.sigma))
			pars[idx["s2.sigma", 1]:idx["s2.sigma", 2], 5] = abs(fixed.pars$s2.sigma)
			pars[idx["s2.sigma", 1]:idx["s2.sigma", 2], 6] = abs(fixed.pars$s2.sigma)
		}
	}
	if(modelinc[32]>0){
		.sd = sd(data)
		pars[idx["s3.sigma", 1]:idx["s3.sigma", 2], 5] = 1e-12
		pars[idx["s3.sigma", 1]:idx["s3.sigma", 2], 6] = 10*.sd
		if(is.null(start.pars$s3.sigma)) pars[idx["s3.sigma", 1]:idx["s3.sigma", 2], 1] = .sd else pars[idx["s3.sigma", 1]:idx["s3.sigma", 2], 1] = abs(start.pars$s3.sigma[1])
		if(any(substr(fixed.names, 1, 8) == "s3.sigma")){
			pars[idx["s3.sigma", 1]:idx["s3.sigma", 2], 1] = abs(as.numeric(fixed.pars$s3.sigma))
			pars[idx["s3.sigma", 1]:idx["s3.sigma", 2], 5] = abs(fixed.pars$s3.sigma)
			pars[idx["s3.sigma", 1]:idx["s3.sigma", 2], 6] = abs(fixed.pars$s3.sigma)
		}
	}
	
	if(modelinc[31]>0){
		.sd = sd(data)
		pars[idx["s4.sigma", 1]:idx["s4.sigma", 2], 5] = 1e-12
		pars[idx["s4.sigma", 1]:idx["s4.sigma", 2], 6] = 10*.sd
		if(is.null(start.pars$s4.sigma)) pars[idx["s4.sigma", 1]:idx["s4.sigma", 2], 1] = .sd else pars[idx["s4.sigma", 1]:idx["s4.sigma", 2], 1] = abs(start.pars$s4.sigma[1])
		if(any(substr(fixed.names, 1, 8) == "s4.sigma")){
			pars[idx["s4.sigma", 1]:idx["s4.sigma", 2], 1] = abs(as.numeric(fixed.pars$s4.sigma))
			pars[idx["s4.sigma", 1]:idx["s4.sigma", 2], 5] = abs(fixed.pars$s4.sigma)
			pars[idx["s4.sigma", 1]:idx["s4.sigma", 2], 6] = abs(fixed.pars$s4.sigma)
		}
	}
	
	
	return( pars )
}
# apARCH model start parameters
# index match to rugarch:
# modelinc+21
.aparchstart = function(pars, arglist)
{
	data = arglist$data
	model = arglist$model
	# little trick to copy and paste the indices which are used in rugarch
	modelinc = model$modelinc[-c(1:24)]
	start.pars = model$start.pars
	fixed.pars = model$fixed.pars
	idx = model$pidx
	fixed.names = names(fixed.pars)
	start.names = names(start.pars)
	if(modelinc[7]>0){
		if(is.na(pars[idx["omega", 1]:idx["omega", 2], 5])) pars[idx["omega", 1]:idx["omega", 2], 5] = eps
		if(is.na(pars[idx["omega", 1]:idx["omega", 2], 6])) pars[idx["omega", 1]:idx["omega", 2], 6] = var(data)*1000
		if(is.null(start.pars$omega)) pars[idx["omega", 1]:idx["omega", 2], 1] = (var(data, na.rm = TRUE))/1000 else 
			pars[idx["omega", 1]:idx["omega", 2], 1] = start.pars$omega[1]
		if(any(substr(fixed.names, 1, 5) == "omega")){
			pars[idx["omega", 1]:idx["omega", 2], 1] = as.numeric(fixed.pars$omega)
			pars[idx["omega", 1]:idx["omega", 2], 5] = fixed.pars$omega
			pars[idx["omega", 1]:idx["omega", 2], 6] = fixed.pars$omega
		}
	}
	if(modelinc[8]>0){
		gpnames = paste("alpha",1:modelinc[8],sep="")
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 5]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 6]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 6] =  1-TinY
		pars[idx["alpha", 1]:idx["alpha", 2], 1] = rep(0.05/modelinc[8], modelinc[8])
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
	}
	
	if(modelinc[10] > 0){
		gqnames = paste("gamma",1:modelinc[10],sep="")
		pxd = which(is.na(pars[idx["gamma", 1]:idx["gamma", 2], 5]))
		if(length(pxd)>0) pars[(idx["gamma", 1]:idx["gamma", 2])[pxd], 5] = -1+TinY
		pxd = which(is.na(pars[idx["gamma", 1]:idx["gamma", 2], 6]))
		if(length(pxd)>0) pars[(idx["gamma", 1]:idx["gamma", 2])[pxd], 6] =  1-TinY
		pars[idx["gamma", 1]:idx["gamma", 2], 1] = rep(0.05/modelinc[10],modelinc[10])
		sp = na.omit(match(start.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gqnames[sp[i]], 1] = as.numeric(start.pars[gqnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gqnames[sp[i]], 1] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 5] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 6] = as.numeric(fixed.pars[gqnames[sp[i]]])
			}
		}
	}
	
	if(modelinc[13]>0){
		if(is.na(pars[idx["delta", 1]:idx["delta", 2], 5])) pars[idx["delta", 1]:idx["delta", 2], 5] = 0.01
		if(is.na(pars[idx["delta", 1]:idx["delta", 2], 6])) pars[idx["delta", 1]:idx["delta", 2], 6] = 3.5
		if(is.null(start.pars$delta)) pars[idx["delta", 1]:idx["delta", 2], 1] = 2 else pars[idx["delta", 1]:idx["delta", 2], 1] = start.pars$delta[1]
		if(any(substr(fixed.names, 1, 5) == "delta")){
			pars[idx["delta", 1]:idx["delta", 2], 1] = as.numeric(fixed.pars$delta)
			pars[idx["delta", 1]:idx["delta", 2], 5] = as.numeric(fixed.pars$delta)
			pars[idx["delta", 1]:idx["delta", 2], 6] = as.numeric(fixed.pars$delta)
		}
	}
	
	if(modelinc[9] > 0){
		gqnames = paste("beta",1:modelinc[9],sep="")
		pxd = which(is.na(pars[idx["beta", 1]:idx["beta", 2], 5]))
		if(length(pxd)>0) pars[(idx["beta", 1]:idx["beta", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["beta", 1]:idx["beta", 2], 6]))
		if(length(pxd)>0) pars[(idx["beta", 1]:idx["beta", 2])[pxd], 6] =  1-TinY
		pars[idx["beta", 1]:idx["beta", 2], 1] = rep(0.9/modelinc[9], modelinc[9])
		sp = na.omit(match(start.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gqnames[sp[i]], 1] = as.numeric(start.pars[gqnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gqnames[sp[i]], 1] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 5] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 6] = as.numeric(fixed.pars[gqnames[sp[i]]])
			}
		}
	}
	if(modelinc[15]>0){
		vxnames = paste("vxreg",1:modelinc[15],sep="")
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 5]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 6]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 6] = 100
		pars[idx["vxreg", 1]:idx["vxreg", 2], 1] = rep(TinY, modelinc[15])
		sp = na.omit(match(start.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[vxnames[sp[i]], 1] = as.numeric(start.pars[vxnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[vxnames[sp[i]], 1] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 5] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 6] = as.numeric(fixed.pars[vxnames[sp[i]]])
			}
		}
	}	
	return( pars )
}


.egarchstart = function(pars, arglist)
{
	data = arglist$data
	model = arglist$model
	modelinc = model$modelinc[-c(1:24)]
	start.pars = model$start.pars
	fixed.pars = model$fixed.pars
	idx = model$pidx
	fixed.names = names(fixed.pars)
	start.names = names(start.pars)
	if(modelinc[7]>0){
		if(is.na(pars[idx["omega", 1]:idx["omega", 2], 5])) pars[idx["omega", 1]:idx["omega", 2], 5] = -10
		if(is.na(pars[idx["omega", 1]:idx["omega", 2], 6])) pars[idx["omega", 1]:idx["omega", 2], 6] = 10
		if(is.null(start.pars$omega)) pars[idx["omega", 1]:idx["omega", 2], 1] = (var(data, na.rm = TRUE)) else pars[idx["omega", 1]:idx["omega", 2], 1] = start.pars$omega[1]
		if(any(substr(fixed.names, 1, 5) == "omega")){
			pars[idx["omega", 1]:idx["omega", 2], 1] = as.numeric(fixed.pars$omega)
			pars[idx["omega", 1]:idx["omega", 2], 5] = fixed.pars$omega
			pars[idx["omega", 1]:idx["omega", 2], 6] = fixed.pars$omega
		}
	}
	if(modelinc[8]>0){
		gpnames = paste("alpha",1:modelinc[8],sep="")
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 5]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 5] = -10
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 6]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 6] =  10
		pars[idx["alpha", 1]:idx["alpha", 2], 1] = rep(0.05/modelinc[8], modelinc[8])
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
	}
	
	if(modelinc[10] > 0){
		gqnames = paste("gamma",1:modelinc[10],sep="")
		pxd = which(is.na(pars[idx["gamma", 1]:idx["gamma", 2], 5]))
		if(length(pxd)>0) pars[(idx["gamma", 1]:idx["gamma", 2])[pxd], 5] = -10
		pxd = which(is.na(pars[idx["gamma", 1]:idx["gamma", 2], 6]))
		if(length(pxd)>0) pars[(idx["gamma", 1]:idx["gamma", 2])[pxd], 6] =  10
		pars[idx["gamma", 1]:idx["gamma", 2], 1] = rep(0.1/modelinc[10],modelinc[10])
		sp = na.omit(match(start.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gqnames[sp[i]], 1] = as.numeric(start.pars[gqnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gqnames[sp[i]], 1] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 5] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 6] = as.numeric(fixed.pars[gqnames[sp[i]]])
			}
		}
	}
	
	if(modelinc[9] > 0){
		gqnames = paste("beta",1:modelinc[9],sep="")
		pxd = which(is.na(pars[idx["beta", 1]:idx["beta", 2], 5]))
		if(length(pxd)>0) pars[(idx["beta", 1]:idx["beta", 2])[pxd], 5] = -1+TinY
		pxd = which(is.na(pars[idx["beta", 1]:idx["beta", 2], 6]))
		if(length(pxd)>0) pars[(idx["beta", 1]:idx["beta", 2])[pxd], 6] =  1-TinY
		pars[idx["beta", 1]:idx["beta", 2], 1] = rep(0.9/modelinc[9], modelinc[9])
		sp = na.omit(match(start.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gqnames[sp[i]], 1] = as.numeric(start.pars[gqnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gqnames[sp[i]], 1] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 5] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 6] = as.numeric(fixed.pars[gqnames[sp[i]]])
			}
		}
	}
	if(modelinc[15]>0){
		vxnames = paste("vxreg",1:modelinc[15],sep="")
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 5]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 5] = -100
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 6]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 6] = 100		
		pars[idx["vxreg", 1]:idx["vxreg", 2], 1] = rep(TinY, modelinc[15])
		sp = na.omit(match(start.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[vxnames[sp[i]], 1] = as.numeric(start.pars[vxnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[vxnames[sp[i]], 1] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 5] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 6] = as.numeric(fixed.pars[vxnames[sp[i]]])
			}
		}
	}	
	return( pars )
}
##################################################################################
# post-estimation functions
.starmodelpost = function(f, T, timer, convergence, message, hess, arglist)
{
	ipars = arglist$ipars
	model = arglist$model
	modelinc = model$modelinc
	estidx = arglist$estidx
	idx = model$pidx
	data = arglist$data
	arglist$returnType = "llh"
	arglist$transform = FALSE
	fit = list()
	if(is.null(hess)){
		fit$hessian = hessian(func = f, x = ipars[estidx, 1], method="Richardson", method.args=list(zero.tol=.Machine$double.eps), arglist = arglist)
	} else{
		fit$hessian = hess
	}
	fit$cvar = try(solve(fit$hessian), silent = TRUE)
	# (might also fail in which case user will see an error about the inversion failure)
	if(inherits(fit$cvar, "try-error")){
		zz = try(solve(optimHess(ipars[estidx, 1], f, arglist = arglist)), silent=TRUE)
		if(inherits(zz, "try-error")) {
			fit$cvar = NULL
			warning("\ntwinkle-->warning: failed to invert hessian\n")
		} else{
			fit$cvar = zz
		}
	}
	arglist$returnType = "all"
	temp = f(pars = ipars[estidx, 1], arglist = arglist)
	fit$z = temp$z
	fit$LLH = -temp$llh
	fit$log.likelihoods = temp$LHT
	fit$residuals = temp$res
	fit$condm = temp$condm
	fit$constm = temp$constm
	fit$probability = temp$probs
	if(modelinc[47]==0) fit$pmu = temp$pmu
	if(sum(modelinc[31:40])>0){
		if(modelinc[50]==4) fit$sigma = temp$h else fit$sigma = sqrt(temp$h)
	}
	arglist$returnType = "LHT"
	if(sum(ipars[,2])>0){
		pall = ipars[estidx | as.logical(ipars[,2]==1), 1]
		fixed = match(rownames(ipars[ipars[,2]==1, , drop = FALSE]), names(pall))
		fixedn = length(fixed)
		fNA = rep(NA, fixedn)
		nfixedn = length(pall) - fixedn
		fit$coef = pall
		if(is.null(fit$cvar)){
			fit$se.coef = rep(NA, nfixedn)
			fit$tval = rep(NA, nfixedn)
			# change here
			fit$matcoef = matrix(NA, ncol = 4, nrow = length(pall))
			fit$matcoef[-fixed,] = cbind(ipars[estidx, 1], fit$se.coef, fit$tval, rep(NA, nfixedn))
			fit$matcoef[fixed,] = cbind(fit$coef[fixed], fNA, fNA, fNA)
			fit$robust.se.coef = rep(NA, nfixedn)
			fit$robust.tval = rep(NA, nfixedn)
			fit$robust.matcoef = matrix(NA, ncol = 4, nrow = length(fit$coef))
			fit$robust.matcoef[-fixed,] = cbind(fit$coef[-fixed], fit$robust.se.coef, fit$robust.tval, rep(NA, nfixedn))
			fit$robust.matcoef[fixed,] = cbind(fit$coef[fixed], fNA, fNA, fNA)
			fit$hessian.message = "failed to invert hessian"
		} else{
			arglist$returnType="LHT"
			tmp = robustvcv(fun = f, pars = ipars[estidx, 1], nlag = 0, hess = fit$hessian, n = T, arglist = arglist)
			fit$robust.cvar = tmp$vcv
			fit$scores = jacobian(func = f, x = ipars[estidx, 1], method="Richardson", method.args=list(zero.tol=.Machine$double.eps), arglist = arglist)
			colnames(fit$scores) = names(ipars[estidx, 1])
			fit$se.coef = sqrt(diag(abs(fit$cvar)))
			fit$tval = fit$coef[-fixed]/fit$se.coef
			# change here
			fit$matcoef = matrix(NA, ncol = 4, nrow = length(fit$coef))
			fit$matcoef[-fixed,] = cbind(fit$coef[-fixed], fit$se.coef, fit$tval, 2*(1-pnorm(abs(fit$tval))))
			fit$matcoef[fixed,] = cbind(fit$coef[fixed], fNA,fNA, fNA)
			fit$robust.se.coef = sqrt(diag(fit$robust.cvar))
			fit$robust.tval = fit$coef[-fixed]/fit$robust.se.coef
			# change here
			fit$robust.matcoef = matrix(NA, ncol = 4, nrow = length(fit$coef))
			fit$robust.matcoef[-fixed,] = cbind(fit$coef[-fixed], fit$robust.se.coef, fit$robust.tval, 2*(1-pnorm(abs(fit$robust.tval))))
			fit$robust.matcoef[fixed,] = cbind(fit$coef[fixed], fNA, fNA, fNA)
			fit$hessian.message = NULL
		}
	} else{
		fit$coef = ipars[estidx, 1]
		if(is.null(fit$cvar)){
			fit$se.coef = rep(NA,length(fit$coef))
			fit$tval = rep(NA,length(fit$coef))
			# change here
			fit$matcoef = cbind(fit$coef, fit$se.coef,fit$tval, rep(NA,length(fit$coef)))
			fit$robust.se.coef = rep(NA,length(fit$coef))
			fit$robust.tval = rep(NA,length(fit$coef))
			# change here
			fit$robust.matcoef = cbind(fit$coef, fit$robust.se.coef,fit$robust.tval, rep(NA,length(fit$coef)))
			fit$hessian.message = "failed to invert hessian"
		} else{
			nlag=min(floor(1.2*(T)^(1/3)),(T))
			arglist$returnType = "LHT"
			tmp = robustvcv(fun = f, pars = ipars[estidx,1], nlag = nlag, hess = fit$hessian, n = T, arglist = arglist)
			fit$robust.cvar = tmp$vcv
			fit$scores = jacobian(func = f, x = ipars[estidx, 1], method="Richardson", method.args=list(zero.tol=.Machine$double.eps), arglist = arglist) 
			colnames(fit$scores) = names(ipars[estidx, 1])
			fit$se.coef = sqrt(diag(abs(fit$cvar)))
			fit$tval = fit$coef/fit$se.coef
			# change here
			fit$matcoef = cbind(fit$coef, fit$se.coef,
					fit$tval, 2*(1-pnorm(abs(fit$tval))))
			fit$robust.se.coef = sqrt(diag(fit$robust.cvar))
			fit$robust.tval = fit$coef/fit$robust.se.coef
			# change here
			fit$robust.matcoef = cbind(fit$coef, fit$robust.se.coef,fit$robust.tval, 2*(1-pnorm(abs(fit$robust.tval))))
			fit$hessian.message = NULL
		}
	}
	# return the correct indicators to the ipars (changed in the case of fixed pars and fixed.se = TRUE)
	ipars[,2:4] = model$pars[,2:4]
	dimnames(fit$matcoef) = list(names(fit$coef), c(" Estimate",
					" Std. Error", " t value", "Pr(>|t|)"))
	dimnames(fit$robust.matcoef) = list(names(fit$coef), c(" Estimate",
					" Std. Error", " t value", "Pr(>|t|)"))
	fit$fitted.values = data-fit$residuals
	fit$convergence = convergence
	fit$message = message
	fit$kappa = temp$kappa
	fit$persistence = temp$persistence
	fit$timer = timer
	fit$ipars = ipars
	return(fit)
}


####################################################################################
# Imported from rugarch package
# probability that x<0
pneg<-function(ghlambda, shape, skew, cond.density,...)
{
	kappa = try(expr=integrate(.pnegfunE, lower = -Inf, upper = 0, ghlambda, shape, skew, cond.density,...)[[1]],
			silent=TRUE)
	if(inherits(kappa, "try-error")){
		kappa<-NA}
	else{kappa<-integrate(.pnegfunE, lower = -Inf, upper = 0, ghlambda, shape, skew, cond.density,...)[[1]]}
	
	kappa
}

.pnegfunE<-function(x, ghlambda, shape, skew, cond.density,...)
{   # A function implemented by Alexios Ghalanos
	# Compute Expectation Value
	cond.density = cond.density[1]
	if (cond.density == "norm"){
		fun = dnorm(x)
	}
	else if(cond.density == "ged") {
		fun = dged(x, shape = shape)
	}
	else if(cond.density == "std") {
		fun = dstd(x, shape = shape)
	}
	else if(cond.density == "snorm") {
		fun =dsnorm(x, skew = skew)
	}
	else if(cond.density == "sged") {
		fun = dsged(x, shape = shape, skew = skew)
	}
	else if(cond.density == "sstd") {
		fun = dsstd(x, shape = shape, skew = skew)
	}
	else if(cond.density == "nig") {
		fun = dsnig(x, shape = shape, skew = skew)
	}
	else if(cond.density == "ghyp") {
		fun = dsgh(x, shape = shape, skew = skew, lambda = ghlambda)
	}
	else if(cond.density == "jsu") {
		fun = djsu(x, skew = skew, shape = shape)
	}
	else if(cond.density == "ghst") {
		fun = dsghst(x, skew = skew, shape = shape)
	}
	else{
		temp<-paste("d",cond.density,sep="")
		.ddist<-eval(parse(text=paste(temp)))
		fun = .ddist(x,...)
	}
	# Return Value:
	fun
}

egarchKappa<-function(ghlambda, shape, skew, cond.density,...)
{
	kappa = try(expr=integrate(.efunE, lower = -Inf, upper = Inf, ghlambda, shape, skew, cond.density,...)[[1]],
			silent=TRUE)
	if(inherits(kappa, "try-error")){
		kappa<-NA}
	else{kappa<-integrate(.efunE, lower = -Inf, upper = Inf, ghlambda, shape, skew, cond.density,...)[[1]]}
	
	kappa
}

.efunE<-function(x, ghlambda, shape, skew, cond.density,...)
{   # A function implemented by Alexios Ghalanos
	# Compute Expectation Value
	cond.density = cond.density[1]
	if(cond.density == "norm"){
		fun = abs(x) * dnorm(x)
	}
	else if(cond.density == "ged") {
		fun = abs(x) * dged(x, shape = shape)
	}
	else if(cond.density == "std") {
		fun = abs(x) * dstd(x, shape = shape)
	}
	else if(cond.density == "snorm") {
		fun = abs(x) * dsnorm(x, skew = skew)
	}
	else if(cond.density == "sged") {
		fun = abs(x) * dsged(x, shape = shape, skew = skew)
	}
	else if(cond.density == "sstd") {
		fun = abs(x) * dsstd(x, shape = shape, skew = skew)
	}
	else if(cond.density == "nig") {
		fun = abs(x) * dsnig(x, shape = shape, skew = skew)
	}
	else if(cond.density == "ghyp") {
		fun = abs(x) * dsgh(x, shape = shape, skew = skew, lambda = ghlambda)
	}
	else if(cond.density == "jsu") {
		fun = abs(x) * djsu(x, skew = skew, shape = shape)
	}
	else if(cond.density == "ghst") {
		fun = abs(x) * dsghst(x, skew = skew, shape = shape)
	}
	else{
		temp<-paste("d",cond.density,sep="")
		.ddist<-eval(parse(text=paste(temp)))
		fun = abs(x) * .ddist(x,...)
	}
	# Return Value:
	fun
}

lagf.numeric = function(x, lags=1)
{
	n = length(x)
	y = rep(NA, n)
	y[(lags+1):n] = x[1:(n-lags)]
	return(y)
}

lagf.matrix = function(x, L)
{
	x = as.matrix(x)
	m = ncol(x)
	n = nrow(x)
	y = matrix(0, nrow=n, ncol = m)
	for(i in 1:m){
		y[(L[i]+1):n,i] = x[1:(n-L[i]),i]
	}
	return(y)
}

repmat = function(a, n, m)
{
	kronecker(matrix(1, n, m), a)
}
size = function(x, n = NULL)
{
	x = as.matrix(x)
	if(missing(n)) sol = c(n = dim(x)[1], m = dim(x)[2]) else sol = dim(x)[n]
	return(sol)
}

zeros = function(n = 1, m = 1)
{
	if(missing(m)) m = n
	sol = matrix(0, nrow = n, ncol = m)
	return(sol)
}

ones = function(n = 1, m = 1)
{
	if(missing(m)) m = n
	sol = matrix(1, nrow = n, ncol = m)
	return(sol)
}

newlagmatrix = function(x,nlags,xc)
{
	nlags = nlags+1
	xt = size(x, 1);
	newX = rbind(x, zeros(nlags, 1))
	lagmatrix = repmat(newX, nlags, 1)
	lagmatrix = matrix(lagmatrix[1:(size(lagmatrix,1)-nlags)], nrow = (xt+nlags-1), ncol = nlags)
	lagmatrix = lagmatrix[nlags:xt,]
	y = lagmatrix[,1]
	x = lagmatrix[,2:nlags]
	if(xc == 1) x = cbind(ones(size(x,1), 1), x)
	return(data.frame(y = y, x = x))
}
getFFdaily = function(enddate = Sys.Date()){
	# match:
	# include=c("FF","MOM","STR","LTR","IND10"), 
	# name:
	# FF factors
	ff.url <-"http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Research_Data_Factors_daily.zip"
	f <- tempfile()
	download.file(ff.url, f)
	file.list <- unzip(f, list=TRUE)
	ff_daily_factors <- read.fwf(unzip(f, files=as.character(file.list[1,1])), 
			widths=c(8,8,8,8,10), header=FALSE,stringsAsFactors=FALSE, skip=5)
	n = dim(ff_daily_factors)
	ff_daily_factors = ff_daily_factors[-c(n-1, n), ]
	for (i in 2:5) ff_daily_factors[,i] <- na.omit(as.numeric(ff_daily_factors[,i]))
	for (i in 2:5) ff_daily_factors[,i] <- ff_daily_factors[,i]/100
	names(ff_daily_factors) <- c("date", "mktrf", "smb", "hml", "rf")
	ff_daily_factors$date <- as.Date(ff_daily_factors$date, format="%Y%m%d")
	# FF Momentum
	ff.url <- "http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Momentum_Factor_daily.zip"
	f <- tempfile()
	download.file(ff.url, f)
	file.list <- unzip(f, list=TRUE)
	ff_mom_factor <- read.fwf(unzip(f, files=as.character(file.list[1,1])), 
			widths=c(8,8), header=FALSE, stringsAsFactors=FALSE, skip=14)
	n = dim(ff_mom_factor)
	ff_mom_factor = ff_mom_factor[-c(n-1, n), ]
	ff_mom_factor[,2] <- na.omit(as.numeric(ff_mom_factor[,2]))/100
	names(ff_mom_factor) <- c("date", "mom")
	ff_mom_factor$date <- as.Date(ff_mom_factor$date, format="%Y%m%d")
	# FF ST reversal
	ff.url <- "http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_ST_Reversal_Factor_daily.zip"
	f <- tempfile()
	download.file(ff.url, f)
	file.list <- unzip(f, list=TRUE)
	ff_str_factor <- read.fwf(unzip(f, files=as.character(file.list[1,1])), 
			widths=c(8,8), header=FALSE, stringsAsFactors=FALSE, skip=14)
	n = dim(ff_str_factor)
	ff_str_factor = ff_str_factor[-c(n-1, n), ]
	ff_str_factor[,2] <- na.omit(as.numeric(ff_str_factor[,2]))/100
	names(ff_str_factor) <- c("date", "str")
	ff_str_factor$date <- as.Date(ff_str_factor$date, format="%Y%m%d")
	# FF LT reversal
	ff.url <- "http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_LT_Reversal_Factor_daily.zip"
	f <- tempfile()
	download.file(ff.url, f)
	file.list <- unzip(f, list=TRUE)
	ff_ltr_factor <- read.fwf(unzip(f, files=as.character(file.list[1,1])), 
			widths=c(8,8), header=FALSE, stringsAsFactors=FALSE, skip=14)
	n = dim(ff_ltr_factor)
	ff_ltr_factor = ff_ltr_factor[-c(n-1, n), ]
	ff_ltr_factor[,2] <- na.omit(as.numeric(ff_ltr_factor[,2]))/100
	names(ff_ltr_factor) <- c("date", "ltr")
	ff_ltr_factor$date <- as.Date(ff_ltr_factor$date, format="%Y%m%d")
	ff_daily_factors <- merge(ff_daily_factors, ff_mom_factor, by="date", all.x=TRUE)
	ff_daily_factors <- merge(ff_daily_factors, ff_str_factor, by="date", all.x=TRUE)
	ff_daily_factors <- merge(ff_daily_factors, ff_ltr_factor, by="date", all.x=TRUE)	
	ff_daily_factors <- subset(ff_daily_factors, subset=!is.na(date))
	ff_daily_factors = as.xts(as.matrix(ff_daily_factors[,-1]), as.Date(ff_daily_factors[,1]))
	xt = time(ff_daily_factors)
	if(!is.finite(z<-min(which(xt>as.Date(enddate))))){
		dx = paste(as.character(tail(xt, 1)), enddate, "d", sep = "/")
		y  = timeBasedSeq(dx, retclass = "Date")
		y  = y[-1] 
		n = length(y)
		fx = matrix(tail(ff_daily_factors, 1), n, NCOL(ff_daily_factors), byrow = TRUE)
		fx[,c(1:3,5)] = 0
		fx = xts(fx, y)
		colnames(fx) = c( "mktrf", "smb", "hml", "rf", "mom","str","ltr")
		ff_daily_factors = rbind(ff_daily_factors, fx)
	} else{
		ff_daily_factors = ff_daily_factors[1:z,]
	}
	return(ff_daily_factors)
}

ar_root = function(x)
{
	model = x@model
	n = model$modelinc[46]
	maxar = max(idx<-model$modelinc[c("s1.phi", "s2.phi", "s3.phi", "s4.phi")])
	cf = x@model$pars[which(x@model$pars[,"Include"]==1),1]
	if(maxar==0) return(list(use=FALSE, rmat=NA))
	rmat  = matrix(NA, ncol = maxar, nrow=n)
	for(i in 1:n){
		if(idx[i]>0){
			arpars  = cf[model$pos.matrix[paste("s",i,".phi",sep=""),1]:model$pos.matrix[paste("s",i,".phi",sep=""),2]]
			rmat[i,] = Mod(polyroot(c(1,-arpars)))
		}
	}
	rownames(rmat) = paste("state_",1:n,sep="")
	colnames(rmat) = paste("Moduli",1:maxar,sep="")
	
	return(list(use=TRUE, rmat = rmat))
}

ma_root = function(x)
{
	model = x@model
	n = model$modelinc[46]
	maxma = max(idx<-model$modelinc[c("s1.psi", "s2.psi", "s3.psi", "s4.psi")])
	cf = x@model$pars[which(x@model$pars[,"Include"]==1),1]
	if(maxma==0) return(list(use=FALSE, rmat=NA))
	rmat  = matrix(NA, ncol = maxma, nrow=n)
	for(i in 1:n){
		if(idx[i]>0){
			mapars  = cf[model$pos.matrix[paste("s",i,".psi",sep=""),1]:model$pos.matrix[paste("s",i,".psi",sep=""),2]]
			rmat[i,] = Mod(polyroot(c(1,mapars)))
		}
	}
	rownames(rmat) = paste("state_",1:n,sep="")
	colnames(rmat) = paste("Moduli",1:maxma,sep="")
	
	return(list(use=TRUE, rmat = rmat))
}

.checkrec = function(init, T){
	if(is.null(init)){
		type = 1
		n = T
	} else{
		if(is.character(init)){
			if(init == "all"){
				type = 1
				n = T
			} else{
				warning("\nstarfit-->warning: unrecognized option in rec.init...using 'all' option instead.\n")
				type = 1
				n = T
			}
		} else{
			if(init>=1){
				type = 1
				if(init<=T){
					n = round(init, 0)
				} else{
					n = T
					warning("\nstarfit-->warning: rec.init value > n.obs (less out.sample)...setting value to n.obs instead.\n")
				}
			} else{
				if(init>0){
					type = 2
					n = init
				} else{
					warning("\nstarfit-->warning: unrecognized option in rec.init...using 'all' option instead \n")
					type = 1
					n = T
				}
			}
		}
	}
	return(list(type = type, n = n))
}

.extractdata = function(data, warn = FALSE)
{
	xdata = try(as.xts(data), silent = TRUE)
	if(inherits(xdata, "try-error")){
		if(warn) warning("\nrugarch-->warning: data indexing not recognized by xts...coercing to Date with origin 1970-01-01.")
		if(is.data.frame(data) | is.matrix(data)) data = as.numeric(data[,1]) else data = as.numeric(data)
		data = unname(data)
		xdata = xts(data, as.POSIXct(as.Date(seq_along(data), origin="1970-01-01")))
	}
	obj = list()
	obj$data = as.numeric(coredata(xdata))
	obj$index = index(xdata)
	obj$period = median(diff(index(xdata)))
	return(obj)
}

.numeric2xts = function(data){
	data = as.numeric(data)
	return(xts(data, as.Date(1:NROW(data), origin="1950-01-01")))
}

.genxts = function(index0, length.out = 10, period = "days"){
	Z = seq(index0, by = period, length.out=length.out)
	return(Z)
}

.makedate = function(x)
{
	# find the divisor: 4 cases "-", "/", ".", and no divisor
	allc = strsplit(x[1], "")
	
	if(any(allc[[1]] == "-")){
		dt = "-"
		dte = t(apply(as.data.frame(x), 1, FUN=function(z) as.numeric(strsplit(z, dt)[[1]]) ))
		ld = max(apply(dte, 1, FUN = function(x) max(nchar(as.character(x)))))+3		
	} else if(any(allc[[1]] == "/")){
		dt = "/"
		dte = t(apply(as.data.frame(x), 1, FUN=function(z) as.numeric(strsplit(z, dt)[[1]]) ))
		ld = max(apply(dte, 1, FUN = function(x) max(nchar(as.character(x)))))+3
	} else if(any(allc[[1]] == ".")){
		dt = "."
		dte = t(apply(as.data.frame(x), 1, FUN=function(z) as.numeric(strsplit(z, dt)[[1]]) ))
	} else{
		# this is a little more complicated
		ld = length(allc[[1]])
		if(ld==6){
			dte = t(apply(as.data.frame(x), 1, FUN=function(z) 
								as.numeric(c(substr(z, 1,2), substr(z, 3,4), substr(z, 5,6)))))
		} else if(ld==8){
			# 2 cases either the 4 digit year is at the beginning or else at the end
			dte.1 = as.vector(t(apply(as.data.frame(x), 1, FUN=function(z) 
										as.numeric(c(substr(z, 1,2))))))
			dte.2 = as.vector(t(apply(as.data.frame(x), 1, FUN=function(z) 
										as.numeric(c(substr(z, 5,6))))))
			if(all(dte.1>18)){
				dte = t(apply(as.data.frame(x), 1, FUN=function(z) 
									as.numeric(c(substr(z, 1,4), substr(z, 5,6), substr(z, 7,8)))))
			} else if(all(dte.2>18)){
				dte = t(apply(as.data.frame(x), 1, FUN=function(z) 
									as.numeric(c(substr(z, 1,2), substr(z, 3,4), substr(z, 5,8)))))
			} else{
				return(list(status=0))
			}
		} else{
			return(list(status=0))	
		}
	}
	m = 0
	for(i in 1:3){
		if(all(dte[,i]<=12)) m = i
	}
	if(m==0) return(list(status=0))
	sq = 1:3
	sq = sq[-m]
	y = 0 
	for(i in sq){
		if(any(dte[,i]>31)) y = i
	}
	if(y==0) return(list(status=0))
	d = sq[-which(sq==y)]
	dmatrix = cbind(dte[,d], dte[,m], dte[,y])
	colnames(dmatrix) = c("d","m","y")
	if(ld==6){
		ddates = as.Date(paste(dmatrix[,3], dmatrix[,2], dmatrix[,1], sep = "-"), format="%y-%m-%d")
		dformat = "%y-%m-%d"
	} else{
		ddates = as.Date(paste(dmatrix[,3], dmatrix[,2], dmatrix[,1], sep = "-"), format="%Y-%m-%d")
		dformat = "%Y-%m-%d"
	}
	
	
	return(list(datesmat = dmatrix, dates = ddates, dformat = dformat, status=1))
}

.checkallfixed = function( spec ){
	# check that a given spec with fixed parameters
	model = spec@model
	pars = model$pars
	pnames = rownames(pars)
	estpars = pnames[as.logical(pars[,2] * pars[,3] + pars[,3] * pars[,4])]
	return( estpars )
	
}

# Exponential smoothing backcast
backcastv = function(res, T, lambda, delta=2){
	s = mean(res^delta)
	v = (lambda^T)*s + (1-lambda)*sum(lambda^(0:(T-1))*(res^delta))
	return(v)
}

# The following functions are based on on Kevin Sheppard's MFE toolbox
#---------------------------------------------------------------------
neweywestcv = function(data, nlag = NULL, center = TRUE)
{
	# Long-run covariance estimation using Newey-West (Bartlett) weights
	#  if nlag empty=NULL NLAG=min(floor(1.2*T^(1/3)),T)
	N = dim(as.matrix(data))[1]
	if(is.null(nlag)) nlag=min(floor(1.2*N^(1/3)),N)
	if(center) data = apply(data, 2, FUN = function(x) scale(x, center = TRUE, scale = FALSE))
	# weights
	bw = (nlag+1-(seq(0,nlag,by=1)))/(nlag+1)
	cv = 1/N * t(data)%*%data
	for(i in 1:nlag){
		gmi = 1/N * (t(data[(i+1):N,])%*%data[1:(N-i),])
		gpp = gmi + t(gmi)
		cv = cv + bw[i+1]*gpp
	}
	return(cv)
}

robustvcv = function(fun, pars, nlag = 0, hess, n, ...)
{
	k = length(pars)
	#h = apply(as.data.frame(pars), 1, FUN = function(x) max(abs(x*eps^(1/3)), 1e-7))
	h = apply(as.data.frame(pars), 1,FUN = function(z) max(abs(z)*1e-4, 1e-9))
	hplus =  pars+h
	hminus = pars-h
	hparsminus = hparsplus = matrix(pars, ncol = k, nrow = k, byrow = T)
	diag(hparsplus) = hplus
	diag(hparsminus) = hminus
	
	likelihoodsminus = likelihoodsplus = zeros(n, k)
	likelihoodsplus =  apply(hparsplus, 1,  FUN = function(x) fun(x, ...))
	likelihoodsminus = apply(hparsminus, 1, FUN = function(x) fun(x, ...))
	
	scores = zeros(n,k)
	likpm = likelihoodsplus-likelihoodsminus
	scores = likpm/(2*repmat(t(h), n, 1))
	A = hess/n
	hess = A
	Ainv = try( solve(A), silent = TRUE )
	if( inherits(Ainv, "try-error")){
		info = 1
		vcv = NA
	} else{
		if(nlag>0){
			B = neweywestcv(scores, nlag)
			vcv = (Ainv%*%B%*%Ainv)/n
		} else{
			B = cov(scores)
			vcv = (Ainv%*%B%*%Ainv)/n
		}
		info = 0
	}
	return(list(vcv = vcv, scores = scores, info = info))
}


#----------------------------------------------------------------------------------
# From HyperbolicDist of D.Scott
#----------------------------------------------------------------------------------
### Test for whole number, with tolerance for representation
### From post by Tony Plate <tplate_at_acm.org>
is.wholenumber <- function(x, tolerance = .Machine$double.eps^0.5){
	if (!is.numeric(x)){
		return(FALSE)
	} else {
		return(isTRUE(all(abs(x - round(x)) < tolerance)))
	}
}

### Calculate recursion for any generalized hyperbolic distribution
### Christine Yang Dong and David Scott from code by Diethelm Wuertz
ghypMean = function(lambda, alpha, beta, delta, mu) 
{	
	gamma <- sqrt(alpha^2 - beta^2)
	mu + delta * beta * besselRatio(delta * gamma, lambda, 1)/gamma
}

momRecursion = function(order = 12, printMatrix = FALSE) {
	## Description:
	##   Computes the moment coefficients recursively
	## Setting Start Values:
	a <- matrix(rep(0, times = order*order), ncol = order)
	a[1, 1] <- 1
	if (order > 1) {
		a[2, 1] <- 1
	}
	## Compute all Cofficients by Recursion:
	if (order > 1) {
		for (d in 2:order) {
			for (l in 2:d) {
				a[d,l] <- a[d - 1,l - 1] + a[d - 1, l]*(2*l + 1 - d)
			}
		}
	}
	rownames(a) <- paste("order=", 1:order, sep = "")
	colnames(a) <- paste("l=", 1:order, sep = "")
	## Print the matrix:
	if (printMatrix) {
		cat("\n")
		print(a)
		cat("\n")
	}
	for (k in 1:order) {
		L <- trunc((k + 1)/2):k
		M <- 2*L - k
	}
	return(list(a = a[order, L], L = L, M = M,
					lmin = trunc((order + 1)/2)))
}

besselRatio <- function(x, nu, orderDiff, useExpScaled = 700){
	if (x > useExpScaled){
		besselK(x, nu + orderDiff, expon.scaled = TRUE)/
				besselK(x, nu, expon.scaled = TRUE)
	}else{
		besselK(x, nu + orderDiff)/besselK(x, nu)
	}
}

## transfer moments about different locations for any distributions
momChangeAbout = function(order = "all", oldMom, oldAbout, newAbout) {
	if (!is.vector(oldMom)){
		stop("A vector of moments must be supplied")
	}
	if (order == "all") {
		## Compute moment of up to length(oldMom) about location new
		mom <- rep(NA,length(oldMom))
		oldMoment <- c(1,oldMom)
		for (i in 1:length(oldMom)) {
			oldMom <- oldMoment[1:(i+1)]
			binomCoeff <- choose(i, 0:i)
			diffPower <- (oldAbout - newAbout)^(i:0)
			mom[i] <- sum(binomCoeff*diffPower*oldMom)
		}
	} else {
		## Check order is within in the right range
		if (length(oldMom) < order) {
			stop("The length of of the vector oldMom must not be less than the
							value of order")
		}
		if (!is.wholenumber(order)){
			stop("Order must be a whole number")
		}
		if ((order < 0)) {
			stop("Order must be positive")
		}
		## Compute moment of a specific order about location new
		oldMom <- c(1,oldMom)
		oldMom <- oldMom[1:(order+1)]
		binomCoeff <- choose(order, 0:order)
		diffPower <- (oldAbout - newAbout)^(order:0)
		mom <- sum(binomCoeff*diffPower*oldMom)
	}
	
	## Return moment
	return(mom)
}


ghypMom = function(order, lambda, alpha, beta, delta, mu, momType = "raw", about = 0) {
	
	## check order is whole number 
	if (!is.wholenumber(order)){
		stop("Order must be a whole number")
	}
	if ((order < 0)) {
		stop("Order must be positive")
	} 
	
	## check momType
	momType <- as.character(momType)
	momType <- tolower(momType)
	if (momType != "raw" & momType != "central" & momType != "mu") {
		stop ("Unrecognised moment type")
	} 
	
	## unpack parameters	
	gm <- sqrt(alpha^2 - beta^2)
	zeta <- delta*gm
	
	if (order == 0) {
		mom <- 1
	} else {         
		## calculate mu moments     
		muMom <- rep (NA,order)
		for (i in 1:order) {
			a <- momRecursion(order = i) 
			coeff <- a$a              
			betaPow <- a$M        
			deltaPow <- 2*a$L
			zetaPow <- a$L
			lengthZetaPow <- length(zetaPow)
			
			## calculate terms and sum
			muM <- coeff*(delta^deltaPow)*(beta^betaPow)*
					sapply(zetaPow, besselRatio, x = zeta, nu = lambda)/(zeta^zetaPow)
			muMom[i] <- sum(muM)   
		}
	}  
	
	if (about != 0) {                    
		mom <- momChangeAbout(order = order, oldMom = muMom, 
				oldAbout = mu, newAbout = about)
	} else {
		if (momType == "mu") {
			mom = muMom[order]
		} else if (momType == "raw") {
			about <- 0
			mom <- momChangeAbout(order = order, oldMom = muMom, 
					oldAbout = mu, newAbout = about)
		} else if (momType == "central") {
			about <- ghypMean(lambda, alpha, beta, delta, mu)
			mom <- momChangeAbout(order = order, oldMom = muMom, 
					oldAbout = mu, newAbout = about)
		}
	}  
	return(mom)
}