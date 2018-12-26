\name{twinkle-package}
\alias{twinkle-package}
\alias{twinkle}
\title{The twinkle package}
\description{
The twinkle package provides an S4 framework for estimating and working with 
logistic smooth transition ARMAX models with a maximum of 4 states. Models can
be specified via the \code{\link{starspec}} method which returns a 
\code{\linkS4class{STARspec}} object which in turn can be passed to the estimation
method \code{link{starfit}}, filter method \code{\link{starfilter}} or path 
simulation method \code{\link{starpath}}. The estimated object of class 
\code{linkS4{STARfit}} can be passed to the forecast method 
\code{\link{starforecast}} or simulation method \code{\link{starsim}}.
A number of extraction methods are available, uniform across most returned
classes, following the same convention as in related packages such as 
\link{rugarch}.\cr
While the conditional mean equation is defined by the STARMAX dynamics, the
conditional variance may follow either a static model, a mixture (only for 2 
states) or one of 3 GARCH models (\sQuote{sGARCH}, \sQuote{eGARCH} and 
\sQuote{gjrGARCH}).\cr
The same number of conditional distributions are implemented as in the 
rugarch package.\cr
The package vignette provides for a detailed and formal treatment of the 
STARMAX model as well as the finer details of the estimation and forecast
methodology implemented.
}
\details{
\tabular{ll}{
Package: \tab twinkle\cr
Type: \tab Package\cr
Title: \tab Dynamic Smooth Transition AutoRegressive Models\cr
Version: \tab 0.9-1\cr
Date: \tab 2014-03-30\cr
Author: \tab alexios ghalanos <alexios@4dscape.com>\cr
Maintainer: \tab alexios ghalanos <alexios@4dscape.com>\cr
Description: \tab Smooth Transition ARMAX models and extensions\cr
Depends: \tab R (>= 2.10), methods, rugarch\cr
Imports: \tab numDeriv, ks, KernSmooth, Rsolnp, lattice, parallel, zoo, xts, parma, RcppDE, SkewHyperbolic\cr
LazyLoad: \tab yes\cr
License: \tab GPL-3\cr}
}
\section{How to cite this package}{
  Whenever using this package, please cite as\cr
 \preformatted{@Manual{Ghalanos_2014,
 author       = {Alexios Ghalanos},
 title        = {{twinkle}: Smooth Transition ARMAX models and extensions.},
 year         = {2014},
 note         = {R package version 0.9-1.},}}
\emph{}
}
\section{License}{
  The releases of this package is licensed under GPL version 3.
}
\author{Alexios Ghalanos}
\references{
Chan, F. and McAleer, M., 2002, Maximum likelihood estimation of STAR and 
STAR-GARCH models: theory and Monte Carlo evidence, 
\emph{Journal of Applied Econometrics}, \bold{17}, 509--534.\cr
Terasvirta, T., 1994, Specification, estimation, and evaluation of smooth 
transition autoregressive models, \emph{Journal of the American Statistical 
Association}, \bold{89}, 208--218.\cr
Tong, H. and Lim, K., 1980, Threshold autoregression, limit cycles 
and cyclical data, \emph{Journal of the Royal Statistical Society}, 
\bold{Series B(42)}, 245--292.\cr
van Dijk, D. and Terasvirta, T. and Franses, P. H., 2002, Smooth transition 
autoregressive models a survey of recent developments, \emph{Econometric Reviews}, 
\bold{21}, 1--47.\cr
}