\name{STARsim-class}
\docType{class}
\alias{STARsim-class}
\alias{sigma,STARsim-method}
\alias{fitted,STARsim-method}
\alias{quantile,STARsim-method}
\alias{states,STARsim-method}
\alias{show,STARsim-method}
\title{class: Univariate STAR Simulation Class}
\description{Class for the univariate STAR simulation.}
\section{Slots}{
  \describe{
    \item{\code{simulation}:}{Object of class \code{"vector"}
    Holds data on the simulation.}
    \item{\code{model}:}{Object of class \code{"vector"}
    The model specification common to all objects.}
    \item{\code{seed}:}{Object of class \code{"integer"}
    The random seed used.}
  }
}
\section{Methods}{
  \describe{
    \item{sigma}{\code{signature(object = "STARsim")}: 
	Extracts the conditional sigma simulated values as a matrix of size n.sim x m.sim (assuming
	a dynamic model for the variance was used).}
	\item{fitted}{\code{signature(object = "STARsim")}: 
	Extracts the conditional mean simulated values as a matrix of size n.sim x m.sim.}
	\item{quantile}{\code{signature(object = "STARsim", probs="numeric")}: 
	Calculates and returns, given a scalar for the probability (additional argument 
	\dQuote{probs}), the conditional quantile of the simulated object as an
	n.sim by m.sim matrix (with the same type of headings as the sigma and
	fitted methods).}
	\item{states}{\code{signature(object = "STARsim")}: 
	Extracts the state probabilities as a matrix. Optional argument \sQuote{type}
	takes values \dQuote{prob} for the probabilities (default) and \dQuote{condm} for
	the simulated values of each state, while \sQuote{sim} denotes the simulation run
	to extract values from (based on argument m.sim in the \code{\link{starsim}} 
	method).}
    \item{show}{\code{signature(object = "STARsim")}: 
    Simulation summary. }
	 }
}
\note{
The \code{sigma} and \code{fitted} methods are used to extract the matrix of 
simulated conditional sigma and mean values.
}
\author{Alexios Ghalanos}
\seealso{
Class \code{\linkS4class{STARpath}}.
}
\keyword{classes}
