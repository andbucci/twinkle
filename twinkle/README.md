# twinkle #

The twinkle package provides an S4 framework for estimating and working with logistic smooth transition ARMAX models with a maximum of 4 states. Models can
be specified via the *starspec* method which returns a **STARspec** object which in turn can be passed to the estimation method *starfit*, filter method *starfilter* or path simulation method *starpath*. The estimated object of class **STARfit** can be passed to the forecast method *starforecast* or simulation method *starsim*. A number of extraction methods are available, uniform across most returned classes, following the same convention as in related packages such as [rugarch](https://bitbucket.org/alexiosg/rugarch).

While the conditional mean equation is defined by the STARMAX dynamics, the conditional variance may follow either a static model, a mixture (upto 4 states) or one of 3 GARCH models (vanilla, exponential or GJR). 

The same number of conditional distributions are implemented as in the rugarch package.

The package vignette provides for a detailed and formal treatment of the STARMAX model as well as the finer details of the estimation and forecast methodology implemented.