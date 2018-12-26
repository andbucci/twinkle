\name{starpath-methods}
\docType{methods}
\alias{starpath}
\alias{starpath-methods}
\alias{starpath,ANY-method}
\alias{starpath,STARspec-method}
\title{function: Univariate STAR Path Simulation}
\description{
Method for simulating the path of a STAR model. This is a convenience function 
which does not require a fitted object (see note below).
}
\usage{
starpath(spec, n.sim = 1000, n.start = 0, m.sim = 1, presigma = NA, 
    prereturns = NA, preresiduals = NA, rseed = NA, 
    custom.dist = list(name = NA, distfit = NA), xregsim = NULL, 
    vregsim = NULL, ssim = NULL, probsim = NULL, ...)
}
\arguments{
\item{spec}{
A univariate STAR spec object of class \code{\linkS4class{STARspec}} with
the required parameters of the model supplied via the fixed.pars list 
argument or \code{\link{setfixed<-}} method.}
\item{n.sim}{
The simulation horizon.}
\item{n.start}{
The burn-in sample.}
\item{m.sim}{
The number of independent simulations.}
\item{presigma}{
Initialization values for sigma in the STAR-GARCH model (required in this model).}
\item{prereturns}{
Initialization value for returns data (required in all models).}
\item{preresiduals}{
Initialization values for residuals in the STAR-GARCH model (required in this model).}
\item{rseed}{
Optional seeding value(s) for the random number generator. For m.sim>1, it is 
possible to provide either a single seed to initialize all values, or one seed 
per separate simulation (i.e. m.sim seeds). However, in the latter case this may 
result in some slight overhead depending on how large m.sim is.}
\item{custom.dist}{ Optional density with fitted object from which to simulate. 
See notes below for details.}
\item{xregsim}{
List of matrices (size of list m.sim, with each matrix having n.sim rows) of 
simulated external regressor-in-mean data. If the specification object contains 
external regressors in the mean equation, this must be provided else will be
assumed zero.}
\item{vregsim}{
List of matrices (size of list m.sim, with each matrix having n.sim rows) of 
simulated external regressor-in-variance data. If the specification object contains 
external regressors in the variance equation, this must be provided else will be
assumed zero.}
\item{ssim}{
List of matrices (size of list m.sim, with each matrix having n.sim rows) of 
simulated external regressors in state equation data (when \sQuote{type} s was 
used in specification, this is required).}
\item{probsim}{
List of matrices (size of list m.sim, with each matrix having n.sim rows and columns
the number of states) of the state simulated probabilities in the case that 
\sQuote{fixed.prob} was used in the specification.}
\item{...}{.}
}
\value{
A \code{\linkS4class{STARpath}} object containing details of the STAR path 
simulation.
}
\details{
This is a convenience method to allow path simulation of various STAR models 
without the need to supply a fit object as in the \code{\link{starsim}} method. 
Instead, a STAR spec object is required with the fixed model parameters.
}
\author{Alexios Ghalanos}
\seealso{
For specification \code{\link{starspec}}, estimation \code{\link{starfit}}, 
filtering \code{\link{starfilter}}, forecasting \code{\link{starforecast}}, 
simulation from estimated object \code{\link{starsim}}, 
rolling forecast and estimation \code{\link{rollstar}}.
}
\keyword{methods}