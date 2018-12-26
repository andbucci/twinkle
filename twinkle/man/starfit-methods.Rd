\name{starfit-methods}
\docType{methods}
\alias{starfit}
\alias{starfit-methods}
\alias{starfit,ANY-method}
\alias{starfit,STARspec-method}
\title{function: Univariate STARMAX Fitting}
\description{
Method for estimating a variety of univariate STARMAX models.
}
\usage{
starfit(spec, data, out.sample = 0, solver = "optim", solver.control = list(), 
    fit.control = list(stationarity = 0, fixed.se = 0, rec.init = "all"), 
    cluster = NULL, n=25, ...) 
}
\arguments{
\item{data}{A univariate xts data object.}
\item{spec}{A univariate STARspec object of class \code{\linkS4class{STARspec}}.}
\item{out.sample}{
A positive integer indicating the number of periods before the last to keep for 
out of sample forecasting (see details).}
\item{solver}{
One of either \dQuote{optim}, \dQuote{cmaes}, \dQuote{nlminb}, 
or \dQuote{solnp}. All solvers except \sQuote{cmaes} can be prepended 
with \dQuote{ms} to indicate a multistart approach (see notes).}
\item{solver.control}{
Control arguments list passed to relevant optimizer.}
\item{fit.control}{
Control arguments passed to the fitting routine. Stationarity explicitly imposes 
the variance stationarity constraint during optimization of a STAR-GARCH model. 
The fixed.se argument controls whether standard errors should be calculated for 
those parameters which were fixed (through the fixed.pars argument of the 
\code{\link{starspec}} function). The rec.init option determines the type of 
initialization for the variance recursion in the STAR-GARCH model. Valid options 
are \sQuote{all} which uses all the values for the unconditional variance 
calculation, an integer greater than or equal to 1 denoting the number of data 
points to use for the calculation, or a positive numeric value less than one 
which determines the weighting for use in an exponential smoothing backcast.}
\item{cluster}{A pre-initialized cluster object from the parallel package to
use with parameter search (when that is turned on in the \sQuote{solver.control} 
and the multi-start strategy (i.e. running the solvers in parallel).}
\item{n}{If \dQuote{strategy} was chosen for \sQuote{solver}, this indicates
the number of outer iterations to perform.}
\item{...}{.}
}
\value{
A \code{\linkS4class{STARfit}} object containing details of the STAR fit.
}
\details{
The optimization routine first calculates a set of feasible starting 
points which are used to initiate the STAR recursion. The main part of the 
likelihood calculation is performed in C-code for speed. If the 
\dQuote{parsearch} argument is passed and set to TRUE in 
the \sQuote{solver.control} list, the parameter space is searched with respect
to the likelihood function (see the startpars function in the Rsolnp package),
else a the uniform random number generator is used. The \dQuote{rseed} argument
can be passed to the \sQuote{solver.control} list to enable replication of 
starting parameters. When using a multi-start strategy (prepend \dQuote{ms}
to the solver e.g. msoptim) the \dQuote{n.restarts} argument passed to the 
\sQuote{solver.control} list indicates the number of multi-start solutions to
attempt from different starting points. Generally the recommended solver for 
such problems is the BFGS solver in optim (set argument \dQuote{method} 
to \sQuote{BFGS} in the \sQuote{solver.control} list), preferably run using
a multi-start strategy. The upper and lower bounds in the unconstrained solvers 
are imposed using a logistic transformation approach which is found to be quite
robust.\cr
The out.sample option is provided in order to carry out forecast performance 
testing against actual data. If the out.sample option is strictly positive, 
then the routine will fit only N - out.sample (where N is the total data length) 
data points, leaving out.sample points for forecasting and testing using 
forecast performance measures. In the \code{\link{starforecast}} routine 
the n.ahead may also be greater than the out.sample number resulting in a 
combination of out of sample data points matched against actual data and 
some without.\cr
}
\note{
STAR models are highly nonlinear in the state dynamic parameters (and linear when
the state probabilities are fixed). Depending on the type of model chosen (STAR, 
STAR-GARCH or STAR-mixture) and choice of transition variable, lags etc, the model
may require a carefully considered optimization, checking the solution by 
re-estimating the model from different starting points (or just using the multi-start
strategy), and tweaking the convergence criteria of the solver used.
}
\author{Alexios Ghalanos}
\seealso{
For specification \code{\link{starspec}}, filtering \code{\link{starfilter}}, 
forecasting \code{\link{starforecast}}, simulation from fitted object 
\code{\link{starsim}}, simulation from spec \code{\link{starpath}},
rolling forecast and estimation \code{\link{rollstar}}.
}
\examples{
\dontrun{
# van Dijk and Frances Dutch Guilder dataset and model
data(forex)
require(quantmod)
# equivalent ways of using type="y".
# First method is the best since fun can then be used in the simulation (path depedence)
fx = na.locf(forex, fromLast = TRUE)
fxw = fx[which(weekdays(index(forex))=="Wednesday"),4]
dx = ROC(fxw, na.pad=FALSE)*100
	fun = function(x){
	x = as.numeric(x)
	y = runMean(abs(x), n=4)
	y[1:3] = c(abs(x[1]), mean(abs(x[1:2])), mean(abs(x[1:3])))
	return(y)
}	
spec = starspec(
mean.model = list(states = 2, include.intercept = c(1,1), arOrder = c(0,2), 
type = c("y","x")[1], lags = 1, fun = fun))
ctrl=list(maxit=10000, alpha=1, beta=0.4, gamma=1.4, reltol=1e-12, 
trace=1,method="BFGS",n.restarts=3,rseed=5)
fit = starfit(spec, data = dx[1:521], solver = "msoptim", solver.control = ctrl)
}
}
\keyword{methods}