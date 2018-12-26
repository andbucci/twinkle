\name{starsim-methods}
\docType{methods}
\alias{starsim}
\alias{starsim-methods}
\alias{starsim,ANY-method}
\alias{starsim,STARfit-method}
\title{function: Univariate STAR Simulation}
\description{
Method for simulation from univariate STAR models.
}
\usage{
starsim(fit, n.sim = 1000, n.start = 0, m.sim = 1, presigma = NA, 
prereturns = NA, preresiduals = NA, rseed = NA, 
custom.dist = list(name = NA, distfit = NA), xregsim = NULL, 
vregsim = NULL, ssim = NULL, probsim = NULL, ...)
}
\arguments{
\item{fit}{
A univariate STAR fit object of class \code{\linkS4class{STARfit}}.}
\item{n.sim}{
The simulation horizon.}
\item{n.start}{
The burn-in sample.}
\item{m.sim}{
The number of simulations.}
\item{presigma}{
Allows the starting sigma values to be provided by the user.}
\item{prereturns}{
Allows the starting return data to be provided by the user.}
\item{preresiduals}{
Allows the starting residuals to be provided by the user.}
\item{rseed}{
Optional seeding value(s) for the random number generator. For m.sim>1, it is 
possible to provide either a single seed to initialize all values, or one seed 
per separate simulation (i.e. m.sim seeds). However, in the latter case this may 
result in some slight overhead depending on how large m.sim is.}
\item{custom.dist}{
Optional density with fitted object from which to simulate. See notes below for 
details.}
\item{xregsim}{
List of matrices (size of list m.sim, with each matrix having n.sim rows) of 
simulated external regressor-in-mean data. If the fit object contains 
external regressors in the mean equation, this must be provided else will be
assumed zero.}
\item{vregsim}{
List of matrices (size of list m.sim, with each matrix having n.sim rows) of 
simulated external regressor-in-variance data. If the fit object contains 
external regressors in the mean equation, this must be provided else will be
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
A \code{\linkS4class{STARsim}} object containing details of the STAR simulation.
}
\details{
The custom.dist option allows for defining a custom density which exists in the
users workspace with methods for \dQuote{r} (sampling, e.g. rnorm) and \dQuote{d} 
(density e.g. dnorm). It must take a single fit object as its second argument.
Alternatively, custom.dist can take any name in the name slot (e.g.\dQuote{sample}) 
and a matrix in the fit slot with dimensions equal to m.sim (columns) and n.sim (rows).
Note that standardized (0,1) innovations are required.
}
\author{Alexios Ghalanos}
\seealso{
For specification \code{\link{starspec}}, estimation \code{\link{starfit}}, 
filtering \code{\link{starfilter}}, forecasting \code{\link{starforecast}}, 
simulation from spec \code{\link{starpath}}, rolling forecast and 
estimation \code{\link{rollstar}}.
}
\keyword{methods}