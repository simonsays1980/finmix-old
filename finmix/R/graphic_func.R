".check.grDevice" <- function() {
	## title argument ##
	any(names(formals(getOption("device")))
		== "title")
}
