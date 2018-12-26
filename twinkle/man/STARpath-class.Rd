\name{STARpath-class}
\docType{class}
\alias{STARpath-class}
\alias{sigma,STARpath-method}
\alias{fitted,STARpath-method}
\alias{quantile,STARpath-method}
\alias{states,STARpath-method}
\alias{show,STARpath-method}
\title{class: Univariate STAR Path Simulation Class}
\description{Class for the univariate STAR Path simulation.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Methods}{
  \describe{
    \item{sigma}{\code{signature(object = "STARpath")}: 
	Extracts the conditional sigma simulated values as a matrix of size n.sim x m.sim
	(for the STAR-GARCH model).}
	\item{fitted}{\code{signature(object = "STARpath")}: 
	Extracts the conditional mean simulated values as a matrix of size n.sim x m.sim.}
	\item{quantile}{\code{signature(x = "STARpath")}: 
	Calculates and returns, given a scalar for the probability (additional argument 
	\dQuote{probs}), the conditional quantile of the simulated object as an
	n.sim by m.sim matrix (with the same type of headings as the sigma and
	fitted methods).}
	\item{states}{\code{signature(object = "STARpath")}: 
	Extracts the state probabilities as a matrix. Optional argument \sQuote{type}
	takes values \dQuote{prob} for the probabilities (default) and \dQuote{condm} for
	the simulated values of each state, while \sQuote{sim} denotes the simulation run
	to extract values from (based on argument m.sim in the \code{\link{starpath}} 
	method).}
    \item{show}{\code{signature(object = "STARpath")}: path simulation summary. }
	 }
}
\author{Alexios Ghalanos}
\seealso{
Class \code{\linkS4class{STARsim}}.
}
\keyword{classes}