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
# include "twinkle.h"
# include "distributions.h"
// starxfilter(states: 1, 2, 3, 4)(static/dynamic: s,d)
void starxfilters1s(int* model, double *pars, int *idx, double *x, double *res, double *z,
		double *mexdata, double *prob, double *constm, double *condm, double *marx, int *T, double *LHT, double *llh)
{
	int i;
	double lk=0;
	for(i=0; i<*T; i++)
	{
		starxfilter1(model, pars, idx, x, res, mexdata, prob, constm, condm, marx, i, *T);
		z[i] = res[i]/pars[idx[29]];
		LHT[i] = log(garchdistribution(z[i], pars[idx[29]], pars[idx[40]], pars[idx[41]], pars[idx[42]], model[43]));
		lk = lk - LHT[i];
	}
	*llh=lk;
}
void starxfilters2s(int* model, double *pars, int *idx, double *x, double *res, double *z,
		double *mexdata, double *prob, double *constm, double *condm,  double *marx, int *T, double *LHT, double *llh)
{
	int i;
	double lk=0;
	for(i=0; i<*T; i++)
	{
		starxfilter2(model, pars, idx, x, res, mexdata, prob, constm, condm, marx, i, *T);
		z[i] = res[i]/pars[idx[29]];
		LHT[i] = log(garchdistribution(z[i], pars[idx[29]], pars[idx[40]], pars[idx[41]], pars[idx[42]], model[43]));
		lk = lk - LHT[i];
	}
	*llh=lk;
}

void starxfilters3s(int* model, double *pars, int *idx, double *x, double *res, double *z,
		double *mexdata, double *prob, double *constm, double *condm,  double *marx, int *T, double *LHT, double *llh)
{
	int i;
	double lk=0;
	for(i=0; i<*T; i++)
	{
		starxfilter3(model, pars, idx, x, res, mexdata, prob, constm, condm, marx, i, *T);
		z[i] = res[i]/pars[idx[29]];
		LHT[i] = log(garchdistribution(z[i], pars[idx[29]], pars[idx[40]], pars[idx[41]], pars[idx[42]], model[43]));
		lk = lk - LHT[i];
	}
	*llh=lk;
}

void starxfilters4s(int* model, double *pars, int *idx, double *x, double *res, double *z,
		double *mexdata, double *prob, double *constm, double *condm,  double *marx, int *T, double *LHT, double *llh)
{
	int i;
	double lk=0;
	for(i=0; i<*T; i++)
	{
		starxfilter4(model, pars, idx, x, res, mexdata, prob, constm, condm, marx, i, *T);
		z[i] = res[i]/pars[idx[29]];
		LHT[i] = log(garchdistribution(z[i], pars[idx[29]], pars[idx[40]], pars[idx[41]], pars[idx[42]], model[43]));
		lk = lk - LHT[i];
	}
	*llh=lk;
}

void starxfilters1sx(int* model, double *pars, int *idx, double *x, double *res, double *z,
		double *mexdata, double *prob, double *constm, double *condm,  double *marx, int *T)
{
	int i;
	for(i=0; i<*T; i++)
	{
		starxfilter1(model, pars, idx, x, res, mexdata, prob, constm, condm, marx, i, *T);
	}
}
void starxfilters2sx(int* model, double *pars, int *idx, double *x, double *res, double *z,
		double *mexdata, double *prob, double *constm, double *condm,  double *marx, int *T)
{
	int i;
	for(i=0; i<*T; i++)
	{
		starxfilter2(model, pars, idx, x, res, mexdata, prob, constm, condm, marx, i, *T);
	}
}
void starxfilters3sx(int* model, double *pars, int *idx, double *x, double *res, double *z,
		double *mexdata, double *prob, double *constm, double *condm,  double *marx, int *T)
{
	int i;
	for(i=0; i<*T; i++)
	{
		starxfilter3(model, pars, idx, x, res, mexdata, prob, constm, condm, marx, i, *T);
	}
}
void starxfilters4sx(int* model, double *pars, int *idx, double *x, double *res, double *z,
		double *mexdata, double *prob, double *constm, double *condm,  double *marx, int *T)
{
	int i;
	for(i=0; i<*T; i++)
	{
		starxfilter4(model, pars, idx, x, res, mexdata, prob, constm, condm, marx, i, *T);
	}
}

void mixturefilterC(int* model, double *pars, int *idx, double *prob, double *res, double *z,
		double *h, int *T, double *LHT, double *llh)
{
	int i;
	double lk=0;
	int TT = 2*(*T);
	int TTT = 3*(*T);
	switch(model[45]){
	case 2:
		for(i=0; i<*T; i++)
		{
			h[i] = prob[i]*pars[idx[29]] + prob[*T+i]*pars[idx[30]];
			z[i] = res[i]/h[i];
			LHT[i] = log(garchdistribution(z[i], h[i], pars[idx[40]], pars[idx[41]], pars[idx[42]], model[43]));
			lk = lk - LHT[i];
		}
		*llh = lk;
		break;
	case 3:
		for(i=0; i<*T; i++)
		{
			h[i] = prob[i]*pars[idx[29]] + prob[*T+i]*pars[idx[30]] + prob[TT+i]*pars[idx[31]];
			z[i] = res[i]/h[i];
			LHT[i] = log(garchdistribution(z[i], h[i], pars[idx[40]], pars[idx[41]], pars[idx[42]], model[43]));
			lk = lk - LHT[i];
		}
		*llh = lk;
		break;
	case 4:
		for(i=0; i<*T; i++)
		{
			h[i] = prob[i]*pars[idx[29]] + prob[*T+i]*pars[idx[30]] + prob[TT+i]*pars[idx[31]] + prob[TTT+i]*pars[idx[32]];
			z[i] = res[i]/h[i];
			LHT[i] = log(garchdistribution(z[i], h[i], pars[idx[40]], pars[idx[41]], pars[idx[42]], model[43]));
			lk = lk - LHT[i];
		}
		*llh = lk;
		break;
	default:
		for(i=0; i<*T; i++)
		{
			h[i] = prob[i]*pars[idx[29]] + prob[*T+i]*pars[idx[30]];
			z[i] = res[i]/h[i];
			LHT[i] = log(garchdistribution(z[i], h[i], pars[idx[40]], pars[idx[41]], pars[idx[42]], model[43]));
			lk = lk - LHT[i];
		}
		*llh = lk;
		break;
	}
}

void garchfilterC(int* model, double *pars, int *idx, double *prob, double *res, double *e,
		double *z, double *vexdata, double *hEst, double *h, double *nres, double *meanz,
		int *m, int *T, double *LHT, double *llh)
{
	int i;
	double lk=0;
	switch(model[49]){
	case 1:
		for(i=0; i<*m; i++)
		{
			h[i] = *hEst;
			e[i] = res[i] * res[i];
			z[i] = res[i]/sqrt(fabs(h[i]));
			LHT[i] = log(garchdistribution(z[i], sqrt(fabs(h[i])), pars[idx[40]], pars[idx[41]], pars[idx[42]], model[43]));
			lk = lk - LHT[i];
		}
		for (i=*m; i<*T; i++)
		{
			sgarchfilter(model, pars, idx, vexdata, e, *T, i, h);
			e[i] = res[i] * res[i];
			z[i] = res[i]/sqrt(fabs(h[i]));
			LHT[i] = log(garchdistribution(z[i], sqrt(fabs(h[i])), pars[idx[40]], pars[idx[41]], pars[idx[42]], model[43]));
			lk = lk - LHT[i];
		}
		*llh=lk;
		break;
	case 2:
		for( i=0; i<*m; i++ )
		{
			h[i] = *hEst;
			e[i] = res[i] * res[i];
			nres[i] = res[i] < 0.0 ? e[i] : 0.0;
			z[i] = res[i]/sqrt(fabs(h[i]));
			LHT[i] = log(garchdistribution(z[i], sqrt(fabs(h[i])), pars[idx[40]], pars[idx[41]], pars[idx[42]], model[43]));
			lk = lk - LHT[i];
		}
		// Main Loop
		for ( i=*m; i<*T; i++ )
		{
			gjrgarchfilter(model, pars, idx, vexdata, nres, e, *T, i, h);
			e[i] = res[i] * res[i];
			nres[i] = res[i] < 0.0 ? e[i] : 0.0;
			z[i] = res[i]/sqrt(fabs(h[i]));
			LHT[i] = log(garchdistribution(z[i], sqrt(fabs(h[i])), pars[idx[40]], pars[idx[41]], pars[idx[42]], model[43]));
			lk = lk - LHT[i];
		}
		*llh = lk;
		break;
	case 3:
		for( i=0; i<*m; i++ )
		{
			h[i] = *hEst;
			e[i] = res[i] * res[i];
			z[i] = res[i]/sqrt(h[i]);
			LHT[i] = log(garchdistribution(z[i], sqrt(fabs(h[i])), pars[idx[40]], pars[idx[41]], pars[idx[42]], model[43]));
			lk = lk - LHT[i];
		}
		for ( i=*m; i<*T; i++ )
		{
			egarchfilter(model, pars, idx, *meanz, z, vexdata, *T, i, h);
			e[i] = res[i] * res[i];
			z[i] = res[i]/sqrt(h[i]);
			LHT[i] = log(garchdistribution(z[i], sqrt(fabs(h[i])), pars[idx[40]], pars[idx[41]], pars[idx[42]], model[43]));
			lk = lk - LHT[i];
		}
		*llh = lk;
		break;
	case 4:
		for(i=0; i<*T; i++)
		{
			h[i] = prob[i]*pars[idx[29]] + prob[*T+i]*pars[idx[30]];
			z[i] = res[i]/h[i];
			LHT[i] = log(garchdistribution(z[i], h[i], pars[idx[40]], pars[idx[41]], pars[idx[42]], model[43]));
			lk = lk - LHT[i];
		}
		*llh = lk;
		break;
	default:
		for(i=0; i<*m; i++)
		{
			h[i] = *hEst;
			e[i] = res[i] * res[i];
			z[i] = res[i]/sqrt(fabs(h[i]));
			LHT[i] = log(garchdistribution(z[i], sqrt(fabs(h[i])), pars[idx[40]], pars[idx[41]], pars[idx[42]], model[43]));
			lk = lk - LHT[i];
		}
		for (i=*m; i<*T; i++)
		{
			sgarchfilter(model, pars, idx, vexdata, e, *T, i, h);
			e[i] = res[i] * res[i];
			z[i] = res[i]/sqrt(fabs(h[i]));
			LHT[i] = log(garchdistribution(z[i], sqrt(fabs(h[i])), pars[idx[40]], pars[idx[41]], pars[idx[42]], model[43]));
			lk = lk - LHT[i];
		}
		*llh=lk;
		break;
	}
}

void garchsimC(int *model, double *pars, int *idx, double *h, double *z, double *res, double *e,
		double *nres, double *meanz, double *vexdata, int *T, int *m)
{
	int i;
	switch(model[49]){
	case 1:
		for (i=*m; i<*T; i++)
		{
			sgarchfilter(model, pars, idx, vexdata, e, *T, i, h);
			res[i]=pow(h[i], 0.5)*z[i];
			e[i] = res[i]*res[i];
		}
		break;
	case 2:
		for ( i=*m; i<*T; i++ )
		{
			gjrgarchfilter(model, pars, idx, vexdata, nres, e, *T, i, h);
			res[i] = pow( h[i],0.5 ) * z[i];
			e[i] = res[i] * res[i];
			nres[i] = res[i] < 0.0 ? e[i] : 0.0;
		}
		break;
	case 3:
		for ( i=*m; i<*T; i++ )
		{
			egarchfilter(model, pars, idx, *meanz, z, vexdata, *T, i, h);
			res[i] = pow(h[i], 0.5)*z[i];
		}
		break;
	default:
		for (i=*m; i<*T; i++)
		{
			sgarchfilter(model, pars, idx, vexdata, e, *T, i, h);
			res[i]=pow(h[i], 0.5)*z[i];
			e[i] = res[i]*res[i];
		}
		break;
	}
}

void starxsim1(int* model, double *pars, int *idx, double *x, double *res,
		double *prob, double *constm, int *m, int *T)
{
	int j, i;
	double im = (double) model[50];
	for(i=*m; i<*T; i++)
	{
		x[i] = prob[i]*constm[i];
		for( j=0; j<model[1]; j++ )
		{
			if(i>=j){
				x[i]+= prob[i]*(pars[idx[1]+j] * (x[i-(j+1)] - im*constm[i-(j+1)]));
			}
		}
		for( j=0; j<model[3]; j++ )
		{
			if(i>=j){
				x[i]+= prob[i]*(pars[idx[3]+j] * res[i-(j+1)]);
			}
		}
		for( j=0; j<model[16]; j++ )
		{
			if(i>=j){
				x[i]+= pars[idx[16]+j] * res[i-(j+1)];
			}
		}
		x[i]+= res[i];
	}
}

void starxsim2(int* model, double *pars, int *idx, double *x, double *s, double *res,
		double *prob, double *constm, int *m, int *T)
{
	int j, i;
	double im = (double) model[50];
	for(i=*m; i<*T; i++)
	{
		s[i]   = constm[i];
		s[i+*T] = constm[i+*T];
		for( j=0; j<model[1]; j++ )
		{
			if(i>=j) s[i]+= pars[idx[1]+j] * (x[i-(j+1)] - im*constm[i-(j+1)]);
		}
		for( j=0; j<model[5]; j++ )
		{
			if(i>=j) s[i+*T]+= pars[idx[5]+j] * (x[i-(j+1)] - im*constm[*T+i-(j+1)]);
		}
		for( j=0; j<model[3]; j++ )
		{
			if(i>=j) s[i]+= pars[idx[3]+j] * res[i-(j+1)];
		}
		for( j=0; j<model[7]; j++ )
		{
			if(i>=j) s[i+*T]+= pars[idx[7]+j] * res[i-(j+1)];
		}
		for( j=0; j<model[16]; j++ )
		{
			if(i>=j){
				x[i]+= pars[idx[16]+j] * res[i-(j+1)];
			}
		}
		x[i]+= prob[i]*s[i] + prob[i+*T]*s[i+*T];
		x[i]+= res[i];
	}
}
//continue here
void starxsim3(int* model, double *pars, int *idx, double *x, double *s, double *res,
		double *prob, double *constm, int *m, int *T)
{
	int j, i;
	double im = (double) model[50];
	int TT = 2*(*T);
	for(i=*m; i<*T; i++)
	{
		s[i]   = constm[i];
		s[i+*T] = constm[i+*T];
		s[i+TT] = constm[i+TT];
		for( j=0; j<model[1]; j++ )
		{
			if(i>=j) s[i]+= pars[idx[1]+j] * (x[i-(j+1)] - im*constm[i-(j+1)]);
		}
		for( j=0; j<model[5]; j++ )
		{
			if(i>=j) s[i+*T]+= pars[idx[5]+j] * (x[i-(j+1)] - im*constm[*T+i-(j+1)]);
		}
		for( j=0; j<model[9]; j++ )
		{
			if(i>=j) s[i+TT]+= pars[idx[9]+j] * (x[i-(j+1)] - im*constm[TT+i-(j+1)]);
		}
		for( j=0; j<model[3]; j++ )
		{
			if(i>=j) s[i]+= pars[idx[3]+j] * res[i-(j+1)];
		}
		for( j=0; j<model[7]; j++ )
		{
			if(i>=j) s[i+*T]+= pars[idx[7]+j] * res[i-(j+1)];
		}
		for( j=0; j<model[11]; j++ )
		{
			if(i>=j) s[i+TT]+= pars[idx[11]+j] * res[i-(j+1)];
		}
		for( j=0; j<model[16]; j++ )
		{
			if(i>=j){
				x[i]+= pars[idx[16]+j] * res[i-(j+1)];
			}
		}
		x[i]+= prob[i]*s[i] + prob[i+*T]*s[i+*T] + prob[i+TT]*s[i+TT];
		x[i]+= res[i];
	}
}

void starxsim4(int* model, double *pars, int *idx, double *x, double *s, double *res,
		double *prob, double *constm, int *m, int *T)
{
	int j, i;
	double im = (double) model[50];
	int TT = 2*(*T);
	int TTT = 3*(*T);
	for(i=*m; i<*T; i++)
	{
		s[i]   = constm[i];
		s[i+*T] = constm[i+*T];
		s[i+TT] = constm[i+TT];
		for( j=0; j<model[1]; j++ )
		{
			if(i>=j) s[i]+= pars[idx[1]+j] * (x[i-(j+1)] - im*constm[i-(j+1)]);
		}
		for( j=0; j<model[5]; j++ )
		{
			if(i>=j) s[i+*T]+= pars[idx[5]+j] * (x[i-(j+1)] - im*constm[*T+i-(j+1)]);
		}
		for( j=0; j<model[9]; j++ )
		{
			if(i>=j) s[i+TT]+= pars[idx[9]+j] * (x[i-(j+1)] - im*constm[TT+i-(j+1)]);
		}
		for( j=0; j<model[13]; j++ )
		{
			if(i>=j) s[i+TTT]+= pars[idx[13]+j] * (x[i-(j+1)] - im*constm[TTT+i-(j+1)]);
		}
		for( j=0; j<model[3]; j++ )
		{
			if(i>=j) s[i]+= pars[idx[3]+j] * res[i-(j+1)];
		}
		for( j=0; j<model[7]; j++ )
		{
			if(i>=j) s[i+*T]+= pars[idx[7]+j] * res[i-(j+1)];
		}
		for( j=0; j<model[11]; j++ )
		{
			if(i>=j) s[i+TT]+= pars[idx[11]+j] * res[i-(j+1)];
		}
		for( j=0; j<model[15]; j++ )
		{
			if(i>=j) s[i+TTT]+= pars[idx[15]+j] * res[i-(j+1)];
		}
		for( j=0; j<model[16]; j++ )
		{
			if(i>=j){
				x[i]+= pars[idx[16]+j] * res[i-(j+1)];
			}
		}
		x[i]+= prob[i]*s[i] + prob[i+*T]*s[i+*T] + prob[i+TT]*s[i+TT] + prob[i+TTT]*s[i+TTT];
		x[i]+= res[i];
	}
}
