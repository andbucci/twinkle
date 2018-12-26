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
#ifndef TWINKLE_H
#define TWINKLE_H
void starxfilters1s(int* model, double *pars, int *idx, double *x, double *res, double *z,
		double *mexdata, double *prob, double *constm, double *condm, double *marx, int *T,
		double *LHT, double *llh);
void starxfilters2s(int* model, double *pars, int *idx, double *x, double *res, double *z,
		double *mexdata, double *prob, double *constm, double *condm, double *marx,  int *T,
		double *LHT, double *llh);
void starxfilters3s(int* model, double *pars, int *idx, double *x, double *res, double *z,
		double *mexdata, double *prob, double *constm, double *condm, double *marx,  int *T,
		double *LHT, double *llh);
void starxfilters4s(int* model, double *pars, int *idx, double *x, double *res, double *z,
		double *mexdata, double *prob, double *constm, double *condm, double *marx,  int *T,
		double *LHT, double *llh);
void starxfilters1sx(int* model, double *pars, int *idx, double *x, double *res, double *z,
		double *mexdata, double *prob, double *constm, double *condm, double *marx,  int *T);
void starxfilters2sx(int* model, double *pars, int *idx, double *x, double *res, double *z,
		double *mexdata, double *prob, double *constm, double *condm, double *marx,  int *T);
void starxfilters3sx(int* model, double *pars, int *idx, double *x, double *res, double *z,
		double *mexdata, double *prob, double *constm, double *condm, double *marx,  int *T);
void starxfilters4sx(int* model, double *pars, int *idx, double *x, double *res, double *z,
		double *mexdata, double *prob, double *constm, double *condm, double *marx,  int *T);
void mixturefilterC(int* model, double *pars, int *idx, double *prob, double *res, double *z,
		double *h, int *T, double *LHT, double *llh);
void garchfilterC(int* model, double *pars, int *idx, double *prob, double *res, double *e,
		double *z, double *vexdata, double *hEst, double *h, double *nres, double *meanz,
		int *m, int *T, double *LHT, double *llh);
void garchsimC(int *model, double *pars, int *idx, double *h, double *z, double *res,
		double *e, double *nres, double *meanz, double *vexdata, int *T, int *m);
void starxsim1(int* model, double *pars, int *idx, double *x, double *res,
		double *prob, double *constm, int *m, int *T);
void starxsim2(int* model, double *pars, int *idx, double *x, double *s, double *res,
		double *prob, double *constm, int *m, int *T);
void starxsim3(int* model, double *pars, int *idx, double *x, double *s, double *res,
		double *prob, double *constm, int *m, int *T);
void starxsim4(int* model, double *pars, int *idx, double *x, double *s, double *res,
		double *prob, double *constm, int *m, int *T);
void stars2normalderiv1(int* model, double *pars, double *dpars, int *idx, double *x, double *res, double *mexdata,
		double *prob, double *y, double *mpt, double *mpinit, double *constm, double *condm, int *T,
		double *dm1, double *dm2, double *dxi1, double *dxi2, double *dphi1, double *dphi2, double *dc,
		double *dalpha, double *dbeta, double *dsigma);
#endif /* TWINKLE_H */
