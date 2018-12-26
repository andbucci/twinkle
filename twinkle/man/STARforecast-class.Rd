\name{STARforecast-class}
\docType{class}
\alias{STARforecast-class}
\alias{sigma,STARforecast-method}
\alias{fitted,STARforecast-method}
\alias{states,STARforecast-method}
\alias{show,STARforecast-method}
\alias{quantile,STARforecast-method}
\alias{plot,STARforecast,missing-method}
\title{class: Univariate STAR Forecast Class}
\description{Class for the univariate STAR forecast.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Methods}{
  \describe{
   \item{sigma}{\code{signature(x = "STARforecast")}: The n.ahead by n.roll+1
   matrix of conditional sigma forecasts (for GARCH and mixture models). 
   The column names are the T[0] dates.}
   \item{fitted}{\code{signature(x = "STARforecast")}: The n.ahead by n.roll+1
   matrix of conditional mean forecasts. The column names are the T[0] dates.}
   	\item{quantile}{\code{signature(x = "STARforecast")}: 
	Calculates and returns, given a scalar for the probability (additional argument 
	\dQuote{probs}), the conditional quantile of the forecast object as an
	n.ahead by n.roll+1 matrix (with the same type of headings as the sigma and
	fitted methods).}
	\item{states}{\code{signature(object = "STARforecast")}: 
	Extracts the state probabilities as an xts matrix. Optional argument \sQuote{type}
	takes values \dQuote{prob} for the probabilities (default), \dQuote{condm} for
	the forecast values of each state, and anything else will return the untransformed state
	dynamics.}
    \item{show}{\code{signature(object = "STARforecast")}: 
    Forecast summary returning the 0-roll frame only.}
    \item{plot}{\code{signature(x = "STARforecast", y = "missing")}: 
    Forecast plots taking optional argument \sQuote{roll} with valid options
    \dQuote{all} for all the rolling periods or an integer for choosing one of the rolling
    windows (when n.ahead>1 was used in combination with rolling forecasts).}
	 }
}
\author{Alexios Ghalanos}
\keyword{classes}