\name{STARroll-class}
\docType{class}
\alias{STARroll-class}
\alias{as.data.frame,STARroll-method}
\alias{show,STARroll-method}
\alias{resume,STARroll-method}
\alias{quantile,STARroll-method}
\alias{pit,STARroll-method}
\alias{convergence,STARroll-method}
\alias{coef,STARroll-method}
\alias{infocriteria,STARroll-method}
\title{class: Univariate STAR Rolling Forecast Class}
\description{Class for the univariate STAR rolling forecast.}
\section{Slots}{
  \describe{
    \item{\code{forecast}:}{Object of class \code{"vector"} }
    \item{\code{model}:}{Object of class \code{"vector"} }
  }
}
\section{Methods}{
  \describe{
    \item{as.data.frame}{\code{signature(x = "STARroll")}: 
    Extracts various values from object (see note).}
	\item{resume}{\code{signature(object = "STARroll")}:
	Resumes a rolling backtest which has non-converged windows using 
	alternative solver and control parameters.}
    \item{coef}{\code{signature(object = "STARroll")}: 
    Extracts the list of coefficients for each estimated window in the 
    rolling backtest.}
    \item{quantile}{\code{signature(x = "STARroll")}: 
	Calculates and returns, given a vector of probabilities (additional argument 
	\dQuote{probs}), the conditional quantiles of the rolling object as an
	xts matrix.}
	\item{pit}{\code{signature(object = "STARroll")}: 
	Calculates and returns the conditional probability integral transform given the 
	realized data and forecast density.}
    \item{convergence}{\code{signature(object = "STARroll")}: 
    Returns the convergence code for the estimation windows, with 0 indicating
    that all have converged and 1 that there were non-converged windows. In 
    the latter case the \sQuote{nonconverged} attribute is also printed of those 
    windows which failed to converge.}
    \item{show}{\code{signature(object = "STARroll")}:
	Summary.}
	 }
}
\note{
The \code{as.data.frame} extractor method allows the extraction of either the
conditional forecast density or the VaR. It takes additional argument
\code{which} with valid values either \dQuote{density} or \dQuote{VaR}.\cr
The \code{coef} method will return a list of the coefficients and their robust 
standard errors (assuming the keep.coef argument was set to TRUE in the 
\code{\link{rollstar}} function), and the ending date of each estimation window.\cr
}
\author{Alexios Ghalanos}
\keyword{classes}