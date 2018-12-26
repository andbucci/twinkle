\name{starplus}
\alias{trans2fun3d}
\alias{trans2fun2d}
\alias{trans3fun2d}
\title{
  Some extra methods for use with STAR models.
}
\description{
A miscellaneous set of functions for use with the twinkle package.
}
\usage{
trans2fun2d(object, colidx = 1, fixed.values = NULL, doplot = TRUE, ...)
trans2fun3d(object, colidx = c(1,2), fixed.values = NULL, doplot = TRUE, ...)
trans3fun2d(object, colidx = 1, fixed.values = NULL, doplot = TRUE, ...)
}
\arguments{
\item{object}{
A \code{\linkS4class{STARfit}} estimated object.}
\item{colidx}{
The colidx of the state variable (see notes).}
\item{fixed.values}{
The variables (if any) which are not used are set to their mean values. Instead, it is possible
for the user to provide a vector of such values instead (excluding colidx).}
\item{doplot}{
Whether to plot the resulting analysis.
}
\item{...}{.}
}
\details{
These plot methods provide what can best be described as the marginal effect of an explanatory variable in 
the state equation on the probability, leaving all other explanatory variables in that equation fixed 
at their mean value (or any other fixed value passed by the user). The \code{trans2fun2d} plot returns
the 2-state, 2-dimensional plot of 1-variable, the \code{trans2dun3d} the 2-state, 3-dimensional plot
of any 2 variables, and the trans3fun2d the 3-state, 2-dimensional plot of 1-variable.\cr
Investigating how a change in any one variable, ceteris paribus, leads to a change in the state probabilities
should provide useful post-estimation insights.
}
\value{
A list with the estimated x-y (and possibly z) values.
}
\author{Alexios Ghalanos}