## 'datamoments' is a virtual classes from which the corresponding           ##
## datamoments for 'continuous' and 'discrete' and 'classifications' inherit ## 
setClass("datamoments",
	representation(
		mean = "matrix",
		variance = "matrix",
		"VIRTUAL")

setGeneric("show", function(.Object) standardGeneric("show"))

setGeneric("mean", function(.Object) standardGeneric("mean"))
setGeneric("variance", function(.Object) standardGeneric("variance"))
setGeneric("mean<-", function(.Object, value) standardGeneric("mean<-"))
setGeneric("variance<-", function(.Object, value) standardGeneric("variance<-"))



