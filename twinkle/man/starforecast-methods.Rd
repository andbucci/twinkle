\name{starforecast-methods}
\docType{methods}
\alias{starforecast}
\alias{starforecast-methods}
\alias{starforecast,ANY-method}
\alias{starforecast,STARfit-method}
\alias{starforecast,STARspec-method}
\title{function: Univariate STAR Forecasting}
\description{
Method for forecasting from a STAR models.
}
\usage{
starforecast(fitORspec, data = NULL, n.ahead = 1, n.roll = 0, out.sample = 0, 
    external.forecasts = list(xregfor = NULL, vregfor = NULL, 
        sfor = NULL, probfor = NULL), method = c("an.parametric", 
        "an.kernel", "mc.empirical", "mc.parametric", "mc.kernel"), 
    mc.sims = NULL, ...)
}
\arguments{
\item{fitORspec}{
Either a univariate STAR fit object of class \code{\linkS4class{STARfit}} or 
alternatively a univariate STAR specification object of class 
\code{\linkS4class{STARspec}} with valid fixed parameters.}
\item{data}{
Required if a specification rather than a fit object is supplied.}
\item{n.ahead}{
The forecast horizon.}
\item{n.roll}{
The no. of rolling forecasts to create beyond the first one (see details).}
\item{out.sample}{
Optional. If a specification object is supplied, indicates how many data points
to keep for out of sample testing.}
\item{external.forecasts}{
A list with forecasts for the external regressors in the mean and/or variance 
equations if specified, the state dynamics (sfor) in the case when \sQuote{type}
used was \dQuote{s} in the model, and forecast state probabilities in the case when
\sQuote{fixed.probs} was used in the model.}
\item{method}{
The nonlinear nature of the model means that for n.ahead>1 there is no simple closed form
way of obtaining the forecasts. In that case, a number of methods are provided including
the numerical integration using the parametric density (\dQuote{an.parametric}) or a kernel
estimated density on the residuals (\dQuote{an.kernel}), monte carlo integration using the
empirical residuals (\dQuote{mc.empirical}), parametric density (\dQuote{mc.parametric}) and
kernel estimated density (\dQuote{mc.kernel}). See the vignette for further details.
}
\item{mc.sims}{The number of simulations to perform for the monte carlo integration.}
\item{...}{.}
}
\value{
A \code{\linkS4class{STARforecast}} object containing details of the STAR 
forecast. See the class for details on the returned object and methods for 
accessing it and performing some tests.
}
\details{
The forecast function has two dispatch methods allowing the user to call it with 
either a fitted object (in which case the data argument is ignored), or a 
specification object (in which case the data is required) with fixed parameters.\cr
The ability to roll the forecast 1 step at a time is implemented with the 
\code{n.roll} argument which controls how many times to roll the n.ahead 
forecast. The default argument of n.roll = 0 denotes no rolling and returns the 
standard n.ahead forecast. Critically, since n.roll depends on data being 
available from which to base the rolling forecast, the \code{\link{starfit}} 
function needs to be called with the argument \code{out.sample} being at least 
as large as the n.roll argument, or in the case of a specification being used 
instead of a fit object, the \code{out.sample} argument directly in the forecast 
function.\cr
The vignette contains the details on the n.ahead>1 implementation and the interested
reader should definitely consult that.
}
\author{Alexios Ghalanos}
\seealso{
For specification \code{\link{starspec}}, estimation \code{\link{starfit}}, 
filtering \code{\link{starfilter}}, forecasting \code{\link{starforecast}}, 
simulation from estimated object \code{\link{starsim}}, 
simulation from spec \code{\link{starpath}}, rolling forecast and 
estimation \code{\link{rollstar}}.
}
\keyword{methods}