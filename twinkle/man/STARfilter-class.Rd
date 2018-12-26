\name{STARfilter-class}
\docType{class}
\alias{STARfilter-class}
\alias{coef,STARfilter-method}
\alias{fitted,STARfilter-method}
\alias{states,STARfilter-method}
\alias{modelmatrix,STARfilter-method}
\alias{likelihood,STARfilter-method}
\alias{infocriteria,STARfilter-method}
\alias{plot,STARfilter,missing-method}
\alias{residuals,STARfilter-method}
\alias{sigma,STARfilter-method}
\alias{quantile,STARfilter-method}
\alias{pit,STARfilter-method}
\alias{show,STARfilter-method}
\title{class: Univariate STAR Filter Class}
\description{Class for the univariate STAR filter.}
\section{Methods}{
 \describe{
\item{fitted}{\code{signature(object = "STARfilter")}: 
	Extracts the fitted values.}
\item{residuals}{\code{signature(object = "STARfilter")}: 
	Extracts the residuals. Optional logical argument \code{standardize} 
	(default is FALSE) allows to extract the standardized residuals.}
\item{sigma}{\code{signature(object = "STARfilter")}:
	Extracts the conditional sigma values.}
\item{states}{\code{signature(object = "STARfilter")}: 
	Extracts the state probabilities as an xts matrix. Optional argument \sQuote{type}
	takes values \dQuote{prob} for the probabilities (default), \dQuote{condm} for
	the filtered values of each state, and anything else will return the untransformed state
	dynamics.}
\item{modelmatrix}{\code{signature(object = "STARfilter")}: 
	S4 extractor for the model.matrix.}
\item{coef}{\code{signature(object = "STARfilter")}:
	Extracts the coefficients.}
\item{infocriteria}{\code{signature(object = "STARfilter")}: 
	Calculates and returns various information criteria. }
\item{likelihood}{\code{signature(object = "STARfilter")}:
	Extracts the likelihood.}
\item{quantile}{\code{signature(x = "STARfilter")}: 
	Calculates and returns, given a vector of probabilities (additional argument 
	\dQuote{probs}), the conditional quantiles of the filtered object (x).}
\item{pit}{\code{signature(object = "STARfilter")}: 
	Calculates and returns the conditional probability integral transform given the
	data and estimated density.}
\item{plot}{\code{signature(x = "STARfilter", y = "missing")}:
	Filter plots}
\item{show}{\code{signature(object = "STARfilter")}: 
	Filter summary.}
}}
\author{Alexios Ghalanos}
\keyword{classes}