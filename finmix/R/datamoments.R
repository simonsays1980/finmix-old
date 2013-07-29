## Copyright (C) 2013 Lars Simon Zehnder
#
# This file is part of finmix.
#
# finmix is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# finmix is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.

## 'datamoments' is a virtual classes from which the corresponding      ##
## datamoments for 'continuous' and 'discrete' inherit	 		## 
setClass("datamoments",
	representation(
		mean = "matrix",
		var = "matrix",
		data = "data",
		sdatamoments = "sdatamoments",
		"VIRTUAL")
)

## Getters ##
setGeneric("getData", function(.Object) standardGeneric("getData"))
setGeneric("getSDatamoments", function(.Object) standardGeneric("getSDatamoments"))

## Setters ##
## No setters as users should not manipulate a 'datamoments' object ##

## mutual constructor for all type of datamoments ##
"datamoments" <- function(dataset) {
		if(getType(dataset) == "continuous") {
			.Object <- new("cdatamoments", data = dataset)
		}
		else { 
			.Object <- new("ddatamoments", data = dataset)
		}
		return(.Object)
}

