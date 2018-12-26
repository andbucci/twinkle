/*################################################################################
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
#################################################################################*/
# include <R.h>
# include <math.h>
# include "starfilters.h"
// probabilities are calculated in R
void starxfilter1(int* model, double *pars, int *idx, double *x, double *res,
		double *mexdata, double *prob, double *constm, double *condm, double *marx,
		int i, int T)
{
/* --------------------------------------------------------------------------------
 * STARX Process :
 * (1-L)(y[t] - mu[t]) = (phi(1-L)(y[t] - mu[t])).prob[t]
 * where mu[t] includes constant, external and arch-in-mean
 * L is the lag operator
 * --------------------------------------------------------------------------------
 * */
	/*0 constm, 1 condm, 2 res*/
	int j, k, ind;
	double im = (double) model[50];
	if(model[0]>0) constm[i] = pars[idx[0]];
	// Exogenous Regressor Initialization
	if(model[2]>0)
	{
		for(k=0;k<model[2];k++)
		{
			ind=i+(T*k);
			constm[i]+=pars[idx[2]+k]*mexdata[ind];
		}
	}
	condm[i]+=constm[i];
	//AR initialization
	if(model[1]>0)
	{
		for(j=0; j<model[1];j++)
		{
			if(i>=j){
				condm[i]+=pars[idx[1]+j]*(x[i-(j+1)]-im*constm[i-(j+1)]);
			}
		}
	}
	if(model[3]>0)
	{
		for(j=0; j<model[3];j++)
		{
			if(i>=j){
				condm[i]+=pars[idx[3]+j]*res[i-(j+1)];
			}
		}
	}
	if(model[16]>0)
	{
		for(j=0; j<model[16];j++)
		{
			if(i>=j){
				marx[i]+=pars[idx[16]+j]*res[i-(j+1)];
			}
		}
	}
	res[i]=x[i]-(prob[i]*condm[i]) - marx[i];
}
void starxfilter2(int* model, double *pars, int *idx, double *x, double *res,
		double *mexdata, double *prob, double *constm, double *condm, double *marx,
		int i, int T)
{
/* --------------------------------------------------------------------------------
 * STARX Process :
 * (1-L)(y[t] - mu[t]) = (phi(1-L)(y[t] - mu[t])).prob[t]
 * where mu[t] includes constant, external and arch-in-mean
 * L is the lag operator
 * --------------------------------------------------------------------------------
 * */
	/*0 constm, 1 condm, 2 res*/
	int j, k, ind;
	double im = (double) model[50];
	if(model[0]>0) constm[i] = pars[idx[0]];
	if(model[4]>0) constm[i+T] = pars[idx[4]];
	// Exogenous Regressor Initialization
	if(model[2]>0)
	{
		for(k=0;k<model[2];k++)
		{
			ind=i+(T*k);
			constm[i]+=pars[idx[2]+k]*mexdata[ind];
		}
	}
	if(model[6]>0)
	{
		for(k=0;k<model[6];k++)
		{
			ind=i+(T*k);
			constm[i+T]+=pars[idx[6]+k]*mexdata[ind];
		}
	}
	condm[i]+=constm[i];
	condm[i+T]+=constm[i+T];
	//AR initialization
	if(model[1]>0)
	{
		for(j=0; j<model[1];j++)
		{
			if(i>=j){
				condm[i]+=pars[idx[1]+j]*(x[i-(j+1)]-im*constm[i-(j+1)]);
			}
		}
	}
	if(model[5]>0)
	{
			for(j=0; j<model[5];j++)
			{
				if(i>=j){
					condm[i+T]+=pars[idx[5]+j]*(x[i-(j+1)]-im*constm[T+i-(j+1)]);
				}
			}
	}

	if(model[3]>0)
	{
		for(j=0; j<model[3];j++)
		{
			if(i>=j){
				condm[i]+=pars[idx[3]+j]*res[i-(j+1)];
			}
		}
	}
	if(model[7]>0)
	{
		for(j=0; j<model[7];j++)
		{
			if(i>=j){
				condm[i+T]+=pars[idx[7]+j]*res[i-(j+1)];
			}
		}
	}

	if(model[16]>0)
	{
		for(j=0; j<model[16];j++)
		{
			if(i>=j){
				marx[i]+=pars[idx[16]+j]*res[i-(j+1)];
			}
		}
	}
	res[i]=x[i]-(prob[i]*condm[i] + prob[i+T]*condm[i+T]) - marx[i];
}
void starxfilter3(int* model, double *pars, int *idx, double *x, double *res,
		double *mexdata, double *prob, double *constm, double *condm, double *marx,
		int i, int T)
{
/* --------------------------------------------------------------------------------
 * STARX Process :
 * (1-L)(y[t] - mu[t]) = (phi(1-L)(y[t] - mu[t])).prob[t]
 * where mu[t] includes constant, external and arch-in-mean
 * L is the lag operator
 * --------------------------------------------------------------------------------
 * */
	/*0 constm, 1 condm, 2 res*/
	int TT = 2*T;
	int j, k, ind;
	double im = (double) model[50];
	if(model[0]>0) constm[i] = pars[idx[0]];
	if(model[4]>0) constm[i+T] = pars[idx[4]];
	if(model[8]>0) constm[i+TT] = pars[idx[8]];
	// Exogenous Regressor Initialization
	if(model[2]>0)
	{
		for(k=0;k<model[2];k++)
		{
			ind=i+(T*k);
			constm[i]+=pars[idx[2]+k]*mexdata[ind];
		}
	}
	if(model[6]>0)
	{
		for(k=0;k<model[6];k++)
		{
			ind=i+(T*k);
			constm[i+T]+=pars[idx[6]+k]*mexdata[ind];
		}
	}
	if(model[10]>0)
	{
		for(k=0;k<model[10];k++)
		{
			ind=i+(T*k);
			constm[i+TT]+=pars[idx[10]+k]*mexdata[ind];
		}
	}
	condm[i]+=constm[i];
	condm[i+T]+=constm[i+T];
	condm[i+TT]+=constm[i+TT];
	//AR initialization
	if(model[1]>0)
	{
		for(j=0; j<model[1];j++)
		{
			if(i>=j){
				condm[i]+=pars[idx[1]+j]*(x[i-(j+1)]-im*constm[i-(j+1)]);
			}
		}
	}
	if(model[5]>0)
	{
		for(j=0; j<model[5];j++)
		{
			if(i>=j){
				condm[i+T]+=pars[idx[5]+j]*(x[i-(j+1)]-im*constm[T+i-(j+1)]);
			}
		}
	}
	if(model[9]>0)
	{
		for(j=0; j<model[9];j++)
		{
			if(i>=j){
				condm[i+TT]+=pars[idx[9]+j]*(x[i-(j+1)]-im*constm[TT+i-(j+1)]);
			}
		}
	}
	if(model[3]>0)
	{
		for(j=0; j<model[3];j++)
		{
			if(i>=j){
				condm[i]+=pars[idx[3]+j]*res[i-(j+1)];
			}
		}
	}
	if(model[7]>0)
	{
		for(j=0; j<model[7];j++)
		{
			if(i>=j){
				condm[i+T]+=pars[idx[7]+j]*res[i-(j+1)];
			}
		}
	}
	if(model[11]>0)
	{
		for(j=0; j<model[11];j++)
		{
			if(i>=j){
				condm[i+TT]+=pars[idx[11]+j]*res[i-(j+1)];
			}
		}
	}

	if(model[16]>0)
	{
		for(j=0; j<model[16];j++)
		{
			if(i>=j){
				marx[i]+=pars[idx[16]+j]*res[i-(j+1)];
			}
		}
	}

	res[i]=x[i]-(prob[i]*condm[i] + prob[i+T]*condm[i+T] + prob[i+TT]*condm[i+TT]) - marx[i];
}

void starxfilter4(int* model, double *pars, int *idx, double *x, double *res,
		double *mexdata, double *prob, double *constm, double *condm, double *marx,
		int i, int T)
{
/* --------------------------------------------------------------------------------
 * STARX Process :
 * (1-L)(y[t] - mu[t]) = (phi(1-L)(y[t] - mu[t])).prob[t]
 * where mu[t] includes constant, external and arch-in-mean
 * L is the lag operator
 * --------------------------------------------------------------------------------
 * */
	/*0 constm, 1 condm, 2 res*/
	int TT = 2*T;
	int TTT = 3*T;
	int j, k, ind;
	double im = (double) model[50];
	if(model[0]>0) constm[i] = pars[idx[0]];
	if(model[4]>0) constm[i+T] = pars[idx[4]];
	if(model[8]>0) constm[i+TT] = pars[idx[8]];
	if(model[12]>0) constm[i+TTT] = pars[idx[12]];
	// Exogenous Regressor Initialization
	if(model[2]>0)
	{
		for(k=0;k<model[2];k++)
		{
			ind=i+(T*k);
			constm[i]+=pars[idx[2]+k]*mexdata[ind];
		}
	}
	if(model[6]>0)
	{
		for(k=0;k<model[6];k++)
		{
			ind=i+(T*k);
			constm[i+T]+=pars[idx[6]+k]*mexdata[ind];
		}
	}
	if(model[10]>0)
	{
		for(k=0;k<model[10];k++)
		{
			ind=i+(T*k);
			constm[i+TT]+=pars[idx[10]+k]*mexdata[ind];
		}
	}
	if(model[14]>0)
	{
		for(k=0;k<model[14];k++)
		{
			ind=i+(T*k);
			constm[i+TTT]+=pars[idx[14]+k]*mexdata[ind];
		}
	}
	condm[i]+=constm[i];
	condm[i+T]+=constm[i+T];
	condm[i+TT]+=constm[i+TT];
	condm[i+TTT]+=constm[i+TTT];
	//AR initialization
	if(model[1]>0)
	{
		for(j=0; j<model[1];j++)
		{
			if(i>=j){
				condm[i]+=pars[idx[1]+j]*(x[i-(j+1)]-im*constm[i-(j+1)]);
			}
		}
	}
	if(model[5]>0)
	{
		for(j=0; j<model[5];j++)
		{
			if(i>=j){
				condm[i+T]+=pars[idx[5]+j]*(x[i-(j+1)]-im*constm[T+i-(j+1)]);
			}
		}
	}
	if(model[9]>0)
	{
		for(j=0; j<model[9];j++)
		{
			if(i>=j){
				condm[i+TT]+=pars[idx[9]+j]*(x[i-(j+1)]-im*constm[TT+i-(j+1)]);
			}
		}
	}
	if(model[13]>0)
	{
			for(j=0; j<model[13];j++)
			{
				if(i>=j){
					condm[i+TTT]+=pars[idx[13]+j]*(x[i-(j+1)]-im*constm[TTT+i-(j+1)]);
				}
			}
	}

	if(model[3]>0)
	{
		for(j=0; j<model[3];j++)
		{
			if(i>=j){
				condm[i]+=pars[idx[3]+j]*res[i-(j+1)];
			}
		}
	}
	if(model[7]>0)
	{
		for(j=0; j<model[7];j++)
		{
			if(i>=j){
				condm[i+T]+=pars[idx[7]+j]*res[i-(j+1)];
			}
		}
	}
	if(model[11]>0)
	{
		for(j=0; j<model[11];j++)
		{
			if(i>=j){
				condm[i+TT]+=pars[idx[11]+j]*res[i-(j+1)];
			}
		}
	}

	if(model[15]>0)
	{
		for(j=0; j<model[15];j++)
		{
			if(i>=j){
				condm[i+TTT]+=pars[idx[15]+j]*res[i-(j+1)];
			}
		}
	}

	if(model[16]>0)
	{
		for(j=0; j<model[16];j++)
		{
			if(i>=j){
				marx[i]+=pars[idx[16]+j]*res[i-(j+1)];
			}
		}
	}

	res[i]=x[i]-(prob[i]*condm[i] + prob[i+T]*condm[i+T] + prob[i+TT]*condm[i+TT] + prob[i+TTT]*condm[i+TTT]) - marx[i];
}


void sgarchfilter(int *model, double *pars, int *idx, double *vexdata, double *e, int T, int i, double *h)
{
	int j, ind;
	h[i] = h[i] + pars[idx[30]];
	if( model[38]>0 )
	{
		for( j=0; j<model[38]; j++ )
		{
			ind = i + ( T * j );
			h[i] = h[i] + pars[idx[38]+j]*vexdata[ind];
		}
	}
	for( j=0; j<model[31]; j++ )
	{
		h[i] = h[i] + pars[idx[31]+j]*e[i-(j+1)];
	}
	for( j=0; j<model[32]; j++ )
	{
		h[i] = h[i] + pars[idx[32]+j]*h[i-(j+1)];
	}
}

void gjrgarchfilter(int *model, double *pars, int *idx, double *vexdata, double *nres, double *e, int T, int i, double *h)
{
	int j, ind;
	h[i] = h[i] + pars[idx[30]];
	if( model[38]>0 )
	{
		for( j=0; j<model[38]; j++ )
		{
			ind = i + ( T * j );
			h[i] = h[i] + pars[idx[38]+j]*vexdata[ind];
		}
	}
	for( j=0; j<model[31]; j++ )
	{
		h[i] = h[i] + pars[idx[31]+j]*e[i-(j+1)]+pars[idx[33]+j]*nres[i-(j+1)];
	}
	for( j=0; j<model[32]; j++ )
	{
		h[i] = h[i] + pars[idx[32]+j]*h[i-(j+1)];
	}
}

void egarchfilter(int *model, double *pars, int *idx, double meanz, double *z, double *vexdata, int T, int i, double *h)
{
	int j, ind;
	h[i] = h[i] +  pars[idx[30]];
	if( model[38]>0 )
	{
		for( j=0; j<model[38]; j++ )
		{
			ind = i + ( T * j );
			h[i] = h[i] + pars[idx[38]+j]*vexdata[ind];
		}
	}
	for( j=0; j<model[31]; j++ )
	{
		h[i] = h[i] + pars[idx[31]+j]*z[i-(j+1)] + pars[idx[33]+j]*( fabs(z[i-(j+1)]) - meanz );
	}
	for( j=0; j<model[32]; j++ )
	{
		h[i] = h[i] +  pars[idx[32]+j]*log( h[i-(j+1)] );
	}
	h[i] = exp( h[i] );
}

void star3dfun(int *n, double *pars, double *zmat, double *v)
{
	// pars = {const,initp,alpha1,alpha2,beta,gamma}
	int i, j, k;
	for(i=0;i<*n;i++){
		k = i* (*n);
		for(j=0;j<*n;j++){
			zmat[k+j]+=pars[5]*(v[i]*pars[2]+v[*n+j]*pars[3]-pars[0]);
			zmat[k+j]+=pars[4]*pars[1];
			zmat[k+j]=1.0/(1.0+exp(-zmat[j+k]));
		}
	}
}

void rfilter(double *x, double *filter, double *out, int *n)
{
	// n[0] length(x) n[1] length(filter)
    int i, j;
    double sum, tmp;
    for(i = 0; i < n[0]; i++) {
    	sum = x[i];
    	for (j = 0; j < n[1]; j++) {
    		tmp = out[n[1] + i - j - 1];
    		sum += tmp * filter[j];
    	}
    	out[n[1] + i] = sum;
    }
}

