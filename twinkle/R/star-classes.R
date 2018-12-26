#################################################################################
##
##   R package twinkle by Alexios Ghalanos Copyright (C) 2014.
##   This file is part of the R package twinkle.
##
##   The R package twinkle is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package twinkle is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################
#----------------------------------------------------------------------------------
# univariate spec class
#----------------------------------------------------------------------------------
setClass("STARspec",
		representation(
				model = "vector"))
#----------------------------------------------------------------------------------
# univariate fit class
#----------------------------------------------------------------------------------
setClass("STARfit",
		representation(
				fit = "vector",
				model = "vector"))
#----------------------------------------------------------------------------------
# univariate filter class
#----------------------------------------------------------------------------------
setClass("STARfilter",
		representation(
				filter = "vector",
				model = "vector"))
#----------------------------------------------------------------------------------
# univariate forecast class
#----------------------------------------------------------------------------------
setClass("STARforecast",
		representation(
				forecast = "vector",
				model = "vector"))
#----------------------------------------------------------------------------------
# univariate simulation class
#----------------------------------------------------------------------------------
setClass("STARsim",
		representation(
				simulation = "vector",
				model = "vector",
				seed = "integer"))
#----------------------------------------------------------------------------------
# univariate path simulation class
#----------------------------------------------------------------------------------
setClass("STARpath",
		representation(path = "vector",
				model = "vector",
				seed = "integer"))
#----------------------------------------------------------------------------------
# univariate roll class
#----------------------------------------------------------------------------------
setClass("STARroll",
		representation(
				model = "vector",
				forecast = "vector"))