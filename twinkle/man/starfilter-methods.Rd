\name{starfilter-methods}
\docType{methods}
\alias{starfilter}
\alias{starfilter-methods}
\alias{starfilter,ANY-method}
\alias{starfilter,STARspec-method}
\title{function: Univariate STAR Filtering}
\description{
Method for filtering univariate STAR models.
}
\usage{
starfilter(spec, data, out.sample = 0, n.old=NULL, rec.init = 'all', ...)
}
\arguments{
\item{data}{
A univariate xts data object.}
\item{spec}{
A univariate STAR spec object of class \code{\linkS4class{STARspec}} with the 
fixed.pars argument having the model parameters on which the filtering is
to take place.}
\item{out.sample}{
A positive integer indicating the number of periods before the last to keep for 
out of sample forecasting (as in \code{\link{starfit}} function).}
\item{n.old}{
For comparison with STARfit models using the out.sample argument, this is the 
length of the original dataset (see details).}
\item{rec.init}{ The recursion initialization method for the STAR-GARCH model 
(see \code{\link{starfit}} for explanation).}
\item{...}{.}
}
\value{
A \code{\linkS4class{STARfilter}} object containing details of the STAR filter.
}
\details{
The n.old argument is optional and indicates the length of the original data (in 
cases when this represents a series augmented by newer data). The reason for 
using this is so that the old and new datasets agree since the original 
recursion uses the sum of the residuals (in the STAR-GARCH model) to start 
the recursion and therefore is influenced by new data. For a small augmentation 
the values converge after x periods, but it is sometimes preferable to have this 
option so that there is no forward looking information contaminating the study.
}
\author{Alexios Ghalanos}
\seealso{
For specification \code{\link{starspec}}, estimation \code{\link{starfit}}, 
forecasting \code{\link{starforecast}}, simulation from fitted object 
\code{\link{starsim}}, simulation from spec \code{\link{starpath}},
rolling forecast and estimation \code{\link{rollstar}}.
}
\keyword{methods}