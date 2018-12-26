\name{STARfit-class}
\docType{class}
\alias{STARfit-class}
\alias{residuals,STARfit-method}
\alias{modelmatrix,STARfit,xts-method}
\alias{modelmatrix,STARfit,missing-method}
\alias{modelmatrix,ANY,ANY-method}
\alias{modelmatrix}
\alias{states,ANY-method}
\alias{states,STARfit-method}
\alias{states}
\alias{score,ANY-method}
\alias{score,STARfit-method}
\alias{score}
\alias{fitted,STARfit-method}
\alias{vcov,STARfit-method}
\alias{coef,STARfit-method}
\alias{infocriteria,STARfit-method}
\alias{getspec,STARfit-method}
\alias{sigma,STARfit-method}
\alias{quantile,STARfit-method}
\alias{pit,STARfit-method}
\alias{likelihood,STARfit-method}
\alias{convergence,STARfit-method}
\alias{plot,STARfit,missing-method}
\alias{show,STARfit-method}
\alias{nonlinearTest,STARfit,xts-method}
\alias{nonlinearTest,STARfit,missing-method}
\alias{nonlinearTest}
\title{class: Univariate STAR Fit Class}
\description{Class for the univariate STAR fit.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{fit}:}{Object of class \code{"vector"}
    Holds data on the fitted model.}
    \item{\code{model}:}{Object of class \code{"vector"}
    The model specification common to all objects.}
  }
}
\section{Methods}{
\describe{
\item{coef}{\code{signature(object = "STARfit")}: 
Extracts the coefficients. }
\item{vcov}{\code{signature(object = "STARfit")}: 
Extracts the covariance matrix of the parameters. Additional logical option of 
\sQuote{robust} indicates whether to extract the robust based covariance matrix. }
\item{infocriteria}{\code{signature(object = "STARfit")}: 
Calculates and returns various information criteria. }
\item{likelihood}{\code{signature(object = "STARfit")}: 
Extracts the likelihood. }
\item{states}{\code{signature(object = "STARfit")}: 
Extracts the state probabilities as an xts matrix. Optional argument \sQuote{type}
takes values \dQuote{prob} for the probabilities (default), \dQuote{condm} for
the fitted values of each state, and anything else will return the untransformed state
dynamics.}
\item{sigma}{\code{signature(object = "STARfit")}: 
Extracts the conditional sigma values (depending on the type of model).}
\item{modelmatrix}{\code{signature(object = "STARfit",data="missing")}: 
S4 extractor for the model.matrix.}
\item{score}{\code{signature(object = "STARfit")}: 
S4 extractor for the score (jacobian) matrix.}
\item{fitted}{\code{signature(object = "STARfit")}: 
Extracts the fitted values.}
\item{residuals}{\code{signature(object = "STARfit")}: 
Extracts the residuals. Optional logical argument \code{standardize} 
(default is FALSE) allows to extract the standardized residuals.}
\item{getspec}{\code{signature(object = "STARfit")}: 
Extracts and returns the STAR specification from a fit object.}
\item{convergence}{\code{signature(object = "STARfit")}: 
Returns the solver convergence code for the fitted object (zero denotes
convergence). }
\item{quantile}{\code{signature(x = "STARfit")}: 
Calculates and returns, given a vector of probabilities (additional argument 
\dQuote{probs}), the conditional quantiles of the fitted object (x).}
\item{pit}{\code{signature(object = "STARfit")}: 
Calculates and returns the conditional probability integral transform given the
data and estimated density.}
\item{nonlinearTest}{\code{signature(object = "STARfit", data = "missing")}: 
Luukkonen, Saikkonen and Terasvirta (1988) test of STAR nonlinearity against
linearity. Takes additional argument \sQuote{robust} indicating whether to 
adjust the test to account for heteroscedasticity in the errors.}
\item{plot}{\code{signature(x = "STARfit", y = "missing")}: 
Fit plots. }
\item{show}{\code{signature(object = "STARfit")}: 
Fit summary. }
}
}
\author{Alexios Ghalanos}
\keyword{classes}