".setUp.data" <- function(withInd = FALSE) {
	## Get path ##
	data.path <- paste(path.package("finmix"),
			"/tests/data/poisson.data.csv", sep = "")
	data <- read.csv(data.path, header = FALSE, sep = ",")
	if(withInd) {
		ind.path <- paste(path.package("finmix"),
			"tests/data/poisson.ind.csv", sep = "")
		ind <- read.csv(ind.path, header = FALSE, sep = ",")
		data <- data(y. = data, S. = ind, type. = "discrete", 
				r. = 1, N. = nrow(data), sim. = TRUE, 
				bycolumn. = TRUE)
		return(data)
	}
	else {
		data <- data(y. = data, type. = "discrete", r. = 1,
				N. = nrow(data), sim. = TRUE,
				bycolumn. = TRUE)
		return(data)	
	}
}

".setUp.model" <- function() {
	model <- model(.dist = "poisson", .K = 2)
	return(model)
} 

".setUp.mcmc" <- function() {
	mcmc <- mcmc(burnin. = 0, M. = 100, startpar. = FALSE, 
			storeS. = 0, storepost. = FALSE, 
			.ranperm = FALSE)
	return(mcmc)
}

"test.FIX.2" <- function() {
	## Setting:
	## 	indicfix: TRUE
	data <- .setUp.data(withInd = TRUE)
	model <- .setUp.model()
	setIndicFix(model) <- TRUE
	prior <- priordefine(data, model)
	mcmc <- .setUp.mcmc()
	mcmcout <- mixturemcmc(data,model,prior,mcmc)
	## Test cases ##
	checkTRUE(mcmcout@ranperm == FALSE, "check1")
	checkTRUE(mcmcout@M == mcmc@M, "check2")
	checkEqualsNumeric(ncol(mcmcout@par$lambda), model@K)
	checkEqualsNumeric(nrow(mcmcout@par$lambda), mcmc@M)
}
