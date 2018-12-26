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
#ifndef STARFILTERS_H
#define STARFILTERS_H
void starxfilter1(int* model, double *pars, int *idx, double *x, double *res,
		double *mexdata, double *prob, double *constm, double *condm, double *marx,
		int i, int T);
void starxfilter2(int* model, double *pars, int *idx, double *x, double *res,
		double *mexdata, double *prob, double *constm, double *condm, double *marx,
		int i, int T);
void starxfilter3(int* model, double *pars, int *idx, double *x, double *res,
		double *mexdata, double *prob, double *constm, double *condm, double *marx,
		int i, int T);
void starxfilter4(int* model, double *pars, int *idx, double *x, double *res,
		double *mexdata, double *prob, double *constm, double *condm, double *marx,
		int i, int T);
void sgarchfilter(int *model, double *pars, int *idx, double *vexdata,
		double *e, int T, int i, double *h);
void gjrgarchfilter(int *model, double *pars, int *idx, double *vexdata,
		double *nres, double *e, int T, int i, double *h);
void egarchfilter(int *model, double *pars, int *idx, double meanz,
		double *z, double *vexdata, int T, int i, double *h);
void star3dfun(int *n, double *pars, double *zmat, double *v);
void rfilter(double *x, double *filter, double *out, int *n);
#endif /* STARFILTERS_H */
