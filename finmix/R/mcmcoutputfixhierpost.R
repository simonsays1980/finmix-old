setClass("mcmcoutputfixhierpost",
	contains = c("mcmcoutputfixhier", "mcmcoutputfixpost"),
	validity = function(object) {
			## else: OK
			TRUE
	}
)


