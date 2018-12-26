\name{starspec-methods}
\docType{methods}
\alias{starspec}
\alias{starspec-methods}
\alias{starspec,ANY-method}
\title{function: STARX Specification}
\description{
Method for creating a STARX specification object prior to estimation.
}
\usage{
starspec(mean.model = list(states = 2, include.intercept = c(1,1), 
arOrder = c(1,1), maOrder = c(0,0), matype = "linear", 
statevar = c("y","s"), s = NULL, ylags=1, statear = FALSE, 
yfun = NULL, xreg = NULL, transform = "logis"), 
variance.model = list(dynamic = FALSE, model = "sGARCH", 
garchOrder = c(1, 1), submodel = NULL, vreg = NULL, 
variance.targeting = FALSE), distribution.model = "norm", 
start.pars = list(), fixed.pars = list(), fixed.prob = NULL, ...)
}
\arguments{
\item{mean.model}{
List containing the conditional mean specification:\cr
\code{states} The number of states (1-4) to estimate.\cr
\code{include.intercept} Whether to include an intercept in each state (a vector
of length equal to the number of states).\cr
\code{arOrder} The autoregressive order (ar) for each state (A vector
of length equal to the number of states).\cr
\code{maOrder} The moving average order (ma) for each state (A vector
of length equal to the number of states) or a single number if used 
outside the states.\cr
\code{matype} One of \dQuote{linear} or \dQuote{state} denoting whether
the MA dynamics are state based or not.\cr
\code{statevar} Whether the switching variable(s) is an external matrix of regressors
(in which case \code{s} must be supplied) or the lagged value(s) of the returns (y).\cr
\code{s} The xts matrix of external regressors for the state switching dynamics which
must be aligned with the returns data (y) which will be passed in the estimation or filter
routines.\cr
\code{yfun} A user supplied function which applies a transformation to the value of
\sQuote{y} before being passed to the switching dynamics equation. The function must
take a numeric vector and return a numeric vector of the same size (in future this may 
allow a matrix with a fixed number of columns to be returned).\cr
\code{xreg} A matrix object containing the (common) external 
regressors to include in each state (which is passed in the estimation function).\cr
\code{statear} Whether to use autoregressive dynamics in the switching equation.\cr
\code{ylags} A positive integer vector of lags to use when y is chosen as the state 
variable.\cr
}
\item{variance.model}{
List containing the conditional variance model:\cr
\code{dynamic} Whether to use a GARCH model for the variance (i.e. time varying).\cr
\code{model} The following models are allowed: \dQuote{sGARCH}, \dQuote{eGARCH}
and \dQuote{gjrGARCH}.\cr
\code{garchOrder} The ARCH (q) and GARCH (p) orders.\cr
\code{submodel} NULL (nor currently used).\cr
\code{vreg} A matrix object containing the external regressors to 
include in the variance equation with as many rows as will be included in the
data (which is passed in the estimation function).
\code{variance.targeting} (Logical or Numeric) If logical, indicates whether to use 
variance targeting for the conditional variance intercept \dQuote{omega}, else
if numeric, the value provided is used instead of the unconditional variance for 
the calculation of the intercept (in combination with the persistence value).
}
\item{distribution.model}{
The conditional density to use for the innovations. Valid choices are 
\dQuote{norm} for the normal distibution, \dQuote{snorm} for the skew-normal 
distribution, \dQuote{std} for the student-t, \dQuote{sstd} for the skew-student, 
\dQuote{ged} for the generalized error distribution, \dQuote{sged} for the 
skew-generalized error distribution, \dQuote{nig} for the normal inverse 
gaussian distribution, \dQuote{ghyp} for the Generalized  Hyperbolic, and \dQuote{jsu} 
for Johnson's SU distribution. Note that some of the distributions are taken 
from the fBasics package and implenented locally here for convenience. 
The \dQuote{jsu} distribution is the reparametrized version from the 
\dQuote{gamlss} package.}
\item{start.pars}{
List of staring parameters for the optimization routine. These are not usually 
required unless the optimization has problems converging.}
\item{fixed.pars}{
List of parameters which are to be kept fixed during the optimization. It is 
possible that you designate all parameters as fixed so as to quickly recover 
just the results of some previous work or published work. The optional argument
\dQuote{fixed.se} in the \code{\link{starfit}} function indicates whether to
calculate standard errors for those parameters fixed during the post 
optimization stage.}
\item{fixed.prob}{A matrix of fixed probabilities to use (instead of estimating
the switching equation) with columns equal to the number of states. For the
case of only 1 state this can be interpreted as a weighted ARX model where the
probabilities are time weights.}
\item{...}{ . }
}
\value{
A \code{\linkS4class{STARspec}} object containing details of the STARX
specification.
}
\details{
Details of the model are best understood by reading the vignette which also includes
the naming convention of the parameter and finer details of the options available.
}
\author{Alexios Ghalanos}
\examples{
# a standard specification
spec = starspec()
spec
}
\keyword{methods}