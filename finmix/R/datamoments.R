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
setGeneric("getSDataMoments", function(.Object) standardGeneric("getSDataMoments"))

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

