\name{STARspec-class}
\docType{class}
\alias{STARspec-class}
\alias{setfixed<-,STARspec,vector-method}
\alias{modelmatrix,STARspec-method}
\alias{setstart<-,STARspec,vector-method}
\alias{setbounds<-,STARspec,vector-method}
\alias{nonlinearTest,STARspec,xts-method}
\alias{modelmatrix,STARspec,xts-method}
\alias{show,STARspec-method}
\title{class: Univariate STAR Specification Class}
\description{Class for the univariate STAR specification.}
\section{Slots}{
  \describe{
    \item{\code{model}:}{Object of class \code{"vector"}
     The model specification common to all objects.}
    }
}
\section{Methods}{
  \describe{
	\item{show}{\code{signature(object = "STARspec")}: 
	Specification summary.}
	\item{modelmatrix}{\code{signature(object = "STARspec",data="xts")}: 
	S4 model.matrix extractor.}
	\item{setfixed<-}{\code{signature(object = "STARspec", value = "vector")}: 
	Sets the fixed parameters (which must be supplied as a named list). }
	\item{setstart<-}{\code{signature(object = "STARspec", value = "vector")}: 
	Sets the starting parameters (which must be supplied as a named list). }
	\item{setbounds<-}{\code{signature(object = "STARspec", value = "vector")}: 
	Sets the parameters lower and upper bounds, which must be supplied as a named 
	list with each parameter being a numeric vector of length 2 i.e. 
	"s1.phi0"=c(-1,1)). If the vector is of length 1, then this is assumed to be 
	the lower bound, and the upper bound will be set to its default value prior
	to estimation.}
	\item{nonlinearTest}{\code{signature(object = "STARspec", data = "xts")}: 
	Luukkonen, Saikkonen and Terasvirta (1988) test of STAR nonlinearity against
	linearity. Takes additional argument \sQuote{robust} indicating whether to 
	adjust the test to account for heteroscedasticity in the errors.}
	}
}
\author{Alexios Ghalanos}
\examples{
# Default 2-state STAR Model
spec = starspec()
spec
}
\keyword{classes}