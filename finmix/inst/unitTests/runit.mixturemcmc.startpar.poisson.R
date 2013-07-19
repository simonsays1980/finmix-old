### --- Test Setup --- ###

if(TRUE) {
	## Not really needed, but can be handy 
	## when writing tests 
	library("RUnit")
	library("finmix")
}

".setUp.data.startpar" <- function(withInd = FALSE) {
        ## Get path ##
        data.path <- paste(getwd(), 
			"/../data/poisson.data.csv", sep= "")
        data <- read.csv(data.path, header = FALSE, sep = ",")
        if(withInd) {
                ind.path <- paste(getwd(), 
				"/../data/poisson.ind.csv", sep = "")
                ind <- read.csv(ind.path, header = FALSE, sep = ",")
                data <- data(y. = as.matrix(data), S. = as.matrix(ind), type. = "discrete",
                                r. = 1, N. = nrow(data), sim. = TRUE,
                                bycolumn. = TRUE)
                return(data)
        }
        else {
                data <- data(y. = as.matrix(data), type. = "discrete", r. = 1,
                                N. = nrow(data), sim. = TRUE,
                                bycolumn. = TRUE)
                return(data)
        }
}

".setUp.model.startpar" <- function() {
        model <- model(dist. = "poisson", K. = 2)
        return(model)
}

".setUp.mcmc.startpar" <- function() {
        mcmc <- mcmc(burnin. = 0, M. = 100, startpar. = TRUE,
                        storeS. = 2, storepost. = FALSE,
                        ranperm. = FALSE)
        return(mcmc)
}

### --- Test functions --- ###
"test.FIX.STARTPAR.1"  <- function() {
        ## Setting:
        ##      indicfix: TRUE
       	set.seed(0)       
	data <- .setUp.data.startpar(withInd = FALSE)
        model <- .setUp.model.startpar()
	setK(model) <- 1
        setIndicFix(model) <- TRUE
	prior <- prior(hier = FALSE)
        prior <- priordefine(data, model, varargin = prior)
        mcmc <- .setUp.mcmc.startpar()
	(data ~ model ~ mcmc) %=% mcmcstart(data,model,mcmc)
        mcmcout <- mixturemcmc(data,model,prior,mcmc)
        ## Test cases ##
	checkTrue(identical(class(mcmcout)[1], "mcmcoutputfix"), "check1")
        checkTrue(mcmcout@ranperm == FALSE, "check2")
        checkEquals(mcmcout@M, mcmc@M)
	checkTrue(is.list(mcmcout@par), "check3")
	checkTrue("lambda" %in% names(mcmcout@par), "check4")
        checkEquals(ncol(mcmcout@par$lambda), model@K)
        checkEquals(nrow(mcmcout@par$lambda), mcmc@M)
	checkTrue(!any(is.na(mcmcout@par$lambda)), "check5")
	for(k in 1:model@K) {
		checkTrue(all(abs(mcmcout@par$lambda[,k] - 
			mean(mcmcout@par$lambda[,k])) > 0), 
			paste("check", 5 + k, sep =""))
	}
	checkTrue(is.list(mcmcout@log), "check7")
	checkTrue("mixlik" %in% names(mcmcout@log),"check8")
	checkTrue("mixprior" %in% names(mcmcout@log), "check9")
	checkTrue(!("cdpost" %in% names(mcmcout@log)), "check10")
	checkEquals(nrow(mcmcout@log$mixlik), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixlik), 1)
	checkEquals(nrow(mcmcout@log$mixprior), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixprior), 1)
	checkTrue(!any(is.na(mcmcout@log$mixlik)), "check11")
	checkTrue(!any(is.na(mcmcout@log$mixprior)), "check12")
	checkTrue(all(abs(mcmcout@log$mixlik - mean(mcmcout@log$mixlik)) 
			> 0), "check13")
	checkTrue(all(abs(mcmcout@log$mixprior - mean(mcmcout@log$mixprior)) 
			> 0), "check14")
	checkTrue(!.hasSlot(mcmcout, "weight"), "check15")
	checkTrue(!.hasSlot(mcmcout, "hyper"), "check16")
	checkTrue(!.hasSlot(mcmcout, "post"), "check17")
	checkTrue(!.hasSlot(mcmcout, "entropy"), "check18")
	checkTrue(!.hasSlot(mcmcout, "ST"), "check19")
	checkTrue(!.hasSlot(mcmcout, "S"), "check20")
	checkTrue(!.hasSlot(mcmcout, "NK"), "check21")
	checkTrue(!.hasSlot(mcmcout, "clust"), "check22")
	checkTrue(identical(mcmcout@model, model), "check23")
	checkTrue(identical(mcmcout@prior, prior), "check24")	
}

"test.FIX.STARTPAR.2" <- function() {
        ## Setting:
        ##      indicfix: TRUE
       	set.seed(0)       
	data <- .setUp.data.startpar(withInd = TRUE)
        model <- .setUp.model.startpar()
        setIndicFix(model) <- TRUE
	prior <- prior(hier = FALSE)
        prior <- priordefine(data, model, varargin = prior)
        mcmc <- .setUp.mcmc.startpar()
       	(data ~ model ~ mcmc) %=% mcmcstart(data,model,mcmc)
	mcmcout <- mixturemcmc(data,model,prior,mcmc)
        ## Test cases ##
	checkTrue(identical(class(mcmcout)[1], "mcmcoutputfix"), "check1")
        checkTrue(mcmcout@ranperm == FALSE, "check2")
        checkEquals(mcmcout@M, mcmc@M)
	checkTrue(is.list(mcmcout@par), "check3")
	checkTrue("lambda" %in% names(mcmcout@par), "check4")
        checkEquals(ncol(mcmcout@par$lambda), model@K)
        checkEquals(nrow(mcmcout@par$lambda), mcmc@M)
	checkTrue(!any(is.na(mcmcout@par$lambda)), "check5")
	for(k in 1:model@K) {
		checkTrue(all(abs(mcmcout@par$lambda[,k] - 
			mean(mcmcout@par$lambda[,k])) > 0), 
			paste("check", 5 + k, sep =""))
	}
	checkTrue(is.list(mcmcout@log), "check8")
	checkTrue("mixlik" %in% names(mcmcout@log),"check9")
	checkTrue("mixprior" %in% names(mcmcout@log), "check10")
	checkTrue(!("cdpost" %in% names(mcmcout@log)), "check11")
	checkEquals(nrow(mcmcout@log$mixlik), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixlik), 1)
	checkEquals(nrow(mcmcout@log$mixprior), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixprior), 1)
	checkTrue(!any(is.na(mcmcout@log$mixlik)), "check12")
	checkTrue(!any(is.na(mcmcout@log$mixprior)), "check13")
	checkTrue(all(abs(mcmcout@log$mixlik - mean(mcmcout@log$mixlik)) 
			> 0), "check14")
	checkTrue(all(abs(mcmcout@log$mixprior - mean(mcmcout@log$mixprior)) 
			> 0), "check15")
	checkTrue(!.hasSlot(mcmcout, "weight"), "check16")
	checkTrue(!.hasSlot(mcmcout, "hyper"), "check17")
	checkTrue(!.hasSlot(mcmcout, "post"), "check18")
	checkTrue(!.hasSlot(mcmcout, "entropy"), "check19")
	checkTrue(!.hasSlot(mcmcout, "ST"), "check20")
	checkTrue(!.hasSlot(mcmcout, "S"), "check21")
	checkTrue(!.hasSlot(mcmcout, "NK"), "check22")
	checkTrue(!.hasSlot(mcmcout, "clust"), "check23")
	checkTrue(identical(mcmcout@model, model), "check24")
	checkTrue(identical(mcmcout@prior, prior), "check25")	
}

"test.FIX.STARTPAR.3" <- function() {
        ## Setting:
        ##      indicfix: TRUE
       	set.seed(0)       
	data <- .setUp.data.startpar(withInd = TRUE)
        model <- .setUp.model.startpar()
	setK(model) <- 3
        setIndicFix(model) <- TRUE
	prior <- prior(hier = FALSE)
        prior <- priordefine(data, model, varargin = prior)
        mcmc <- .setUp.mcmc.startpar()
 	(data ~ model ~ mcmc) %=% mcmcstart(data,model,mcmc)
       	mcmcout <- mixturemcmc(data,model,prior,mcmc)
        ## Test cases ##
	checkTrue(identical(class(mcmcout)[1], "mcmcoutputfix"), "check1")
        checkTrue(mcmcout@ranperm == FALSE, "check2")
        checkEquals(mcmcout@M, mcmc@M)
	checkTrue(is.list(mcmcout@par), "check3")
	checkTrue("lambda" %in% names(mcmcout@par), "check4")
        checkEquals(ncol(mcmcout@par$lambda), model@K)
        checkEquals(nrow(mcmcout@par$lambda), mcmc@M)
	checkTrue(!any(is.na(mcmcout@par$lambda)), "check5")
	m <- 0
	for(k in 1:model@K) {
		m <- m + k
		checkTrue(all(abs(mcmcout@par$lambda[,k] - 
			mean(mcmcout@par$lambda[,k])) > 0), 
			paste("check", 5 + m, sep =""))
	}
	checkTrue(is.list(mcmcout@log), "check9")
	checkTrue("mixlik" %in% names(mcmcout@log),"check10")
	checkTrue("mixprior" %in% names(mcmcout@log), "check11")
	checkTrue(!("cdpost" %in% names(mcmcout@log)), "check12")
	checkEquals(nrow(mcmcout@log$mixlik), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixlik), 1)
	checkEquals(nrow(mcmcout@log$mixprior), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixprior), 1)
	checkTrue(!any(is.na(mcmcout@log$mixlik)), "check13")
	checkTrue(!any(is.na(mcmcout@log$mixprior)), "check14")
	checkTrue(all(abs(mcmcout@log$mixlik - mean(mcmcout@log$mixlik)) 
			> 0), "check15")
	checkTrue(all(abs(mcmcout@log$mixprior - mean(mcmcout@log$mixprior)) 
			> 0), "check16")
	checkTrue(!.hasSlot(mcmcout, "weight"), "check17")
	checkTrue(!.hasSlot(mcmcout, "hyper"), "check18")
	checkTrue(!.hasSlot(mcmcout, "post"), "check19")
	checkTrue(!.hasSlot(mcmcout, "entropy"), "check20")
	checkTrue(!.hasSlot(mcmcout, "ST"), "check21")
	checkTrue(!.hasSlot(mcmcout, "S"), "check22")
	checkTrue(!.hasSlot(mcmcout, "NK"), "check23")
	checkTrue(!.hasSlot(mcmcout, "clust"), "check24")
	checkTrue(identical(mcmcout@model, model), "check25")
	checkTrue(identical(mcmcout@prior, prior), "check26")	
}

".test.FIX.HIER.STARTPAR.1"  <- function() {
        ## Setting:
        ##      indicfix: TRUE
       	set.seed(0)       
	data <- .setUp.data.startpar(withInd = FALSE)
        model <- .setUp.model.startpar()
	setK(model) <- 1
        setIndicFix(model) <- TRUE
        prior <- priordefine(data, model)
        mcmc <- .setUp.mcmc.startpar()
       	(data ~ model ~ mcmc) %=% mcmcstart(data,model,mcmc)
	mcmcout <- mixturemcmc(data,model,prior,mcmc)
        ## Test cases ##
	checkTrue(identical(class(mcmcout)[1], "mcmcoutputfixhier"), "check1")
        checkTrue(mcmcout@ranperm == FALSE, "check2")
        checkEquals(mcmcout@M, mcmc@M)
	checkTrue(is.list(mcmcout@par), "check3")
	checkTrue("lambda" %in% names(mcmcout@par), "check4")
        checkEquals(ncol(mcmcout@par$lambda), model@K)
        checkEquals(nrow(mcmcout@par$lambda), mcmc@M)
	checkTrue(!any(is.na(mcmcout@par$lambda)), "check5")
	for(k in 1:model@K) {
		checkTrue(all(abs(mcmcout@par$lambda[,k] - 
			mean(mcmcout@par$lambda[,k])) > 0), 
			paste("check", 5 + k, sep =""))
	}
	checkTrue(is.list(mcmcout@log), "check7")
	checkTrue("mixlik" %in% names(mcmcout@log),"check8")
	checkTrue("mixprior" %in% names(mcmcout@log), "check9")
	checkTrue(!("cdpost" %in% names(mcmcout@log)), "check10")
	checkEquals(nrow(mcmcout@log$mixlik), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixlik), 1)
	checkEquals(nrow(mcmcout@log$mixprior), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixprior), 1)
	checkTrue(!any(is.na(mcmcout@log$mixlik)), "check11")
	checkTrue(!any(is.na(mcmcout@log$mixprior)), "check12")
	checkTrue(all(abs(mcmcout@log$mixlik - mean(mcmcout@log$mixlik)) 
			> 0), "check13")
	checkTrue(all(abs(mcmcout@log$mixprior - mean(mcmcout@log$mixprior)) 
			> 0), "check14")
	checkTrue(!.hasSlot(mcmcout, "weight"), "check15")
	checkTrue(is.list(mcmcout@hyper), "check16")
	checkTrue("b" %in% names(mcmcout@hyper), "check17")
	checkEquals(nrow(mcmcout@hyper$b), mcmc@M)
	checkEquals(ncol(mcmcout@hyper$b), 1)
	checkTrue(!any(is.na(mcmcout@hyper$b)), "check18")
	checkTrue(all(abs(mcmcout@hyper$b - mean(mcmcout@hyper$b)) 
		> 0), "check19")
	checkTrue(!.hasSlot(mcmcout, "post"), "check20")
	checkTrue(!.hasSlot(mcmcout, "entropy"), "check21")
	checkTrue(!.hasSlot(mcmcout, "ST"), "check22")
	checkTrue(!.hasSlot(mcmcout, "S"), "check23")
	checkTrue(!.hasSlot(mcmcout, "NK"), "check24")
	checkTrue(!.hasSlot(mcmcout, "clust"), "check25")
	checkTrue(identical(mcmcout@model, model), "check26")
	checkTrue(identical(mcmcout@prior, prior), "check27")	
}

".test.FIX.HIER.STARTPAR.2" <- function() {
        ## Setting:
        ##      indicfix: TRUE
 	set.seed(0)       
      	data <- .setUp.data.startpar(withInd = TRUE)
        model <- .setUp.model.startpar()
        setIndicFix(model) <- TRUE
        prior <- priordefine(data, model)
        mcmc <- .setUp.mcmc.startpar()
       	(data ~ model ~ mcmc) %=% mcmcstart(data,model,mcmc)
	 mcmcout <- mixturemcmc(data,model,prior,mcmc)
        ## Test cases ##
	checkTrue(identical(class(mcmcout)[1], "mcmcoutputfixhier"), "check1")
        checkTrue(mcmcout@ranperm == FALSE, "check2")
        checkEquals(mcmcout@M, mcmc@M)
	checkTrue(is.list(mcmcout@par), "check3")
	checkTrue("lambda" %in% names(mcmcout@par), "check4")
        checkEquals(ncol(mcmcout@par$lambda), model@K)
        checkEquals(nrow(mcmcout@par$lambda), mcmc@M)
	checkTrue(!any(is.na(mcmcout@par$lambda)), "check5")
	for(k in 1:model@K) {
		checkTrue(all(abs(mcmcout@par$lambda[,k] - 
			mean(mcmcout@par$lambda[,k])) > 0), 
			paste("check", 5 + k, sep =""))
	}
	checkTrue(is.list(mcmcout@log), "check8")
	checkTrue("mixlik" %in% names(mcmcout@log),"check9")
	checkTrue("mixprior" %in% names(mcmcout@log), "check10")
	checkTrue(!("cdpost" %in% names(mcmcout@log)), "check11")
	checkEquals(nrow(mcmcout@log$mixlik), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixlik), 1)
	checkEquals(nrow(mcmcout@log$mixprior), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixprior), 1)
	checkTrue(!any(is.na(mcmcout@log$mixlik)), "check12")
	checkTrue(!any(is.na(mcmcout@log$mixprior)), "check13")
	checkTrue(all(abs(mcmcout@log$mixlik - mean(mcmcout@log$mixlik)) 
			> 0), "check14")
	checkTrue(all(abs(mcmcout@log$mixprior - mean(mcmcout@log$mixprior)) 
			> 0), "check15")
	checkTrue(!.hasSlot(mcmcout, "weight"), "check16")
	checkTrue(is.list(mcmcout@hyper), "check17")
	checkTrue("b" %in% names(mcmcout@hyper), "check18")
	checkEquals(nrow(mcmcout@hyper$b), mcmc@M)
	checkEquals(ncol(mcmcout@hyper$b), 1)
	checkTrue(!any(is.na(mcmcout@hyper$b)), "check19")
	checkTrue(all(abs(mcmcout@hyper$b - mean(mcmcout@hyper$b)) 
			> 0), "check20")
	checkTrue(!.hasSlot(mcmcout, "post"), "check21")
	checkTrue(!.hasSlot(mcmcout, "entropy"), "check22")
	checkTrue(!.hasSlot(mcmcout, "ST"), "check23")
	checkTrue(!.hasSlot(mcmcout, "S"), "check24")
	checkTrue(!.hasSlot(mcmcout, "NK"), "check25")
	checkTrue(!.hasSlot(mcmcout, "clust"), "check26")
	checkTrue(identical(mcmcout@model, model), "check27")
	checkTrue(identical(mcmcout@prior, prior), "check28")	
}

".test.FIX.HIER.STARTPAR.3" <- function() {
        ## Setting:
        ##      indicfix: TRUE
 	set.seed(0)       
        data <- .setUp.data.startpar(withInd = TRUE)
        model <- .setUp.model.startpar()
	setK(model) <- 3
        setIndicFix(model) <- TRUE
        prior <- priordefine(data, model)
        mcmc <- .setUp.mcmc.startpar()
       	(data ~ model ~ mcmc) %=% mcmcstart(data,model,mcmc)
	mcmcout <- mixturemcmc(data,model,prior,mcmc)
        ## Test cases ##
	checkTrue(identical(class(mcmcout)[1], "mcmcoutputfixhier"), "check1")
        checkTrue(mcmcout@ranperm == FALSE, "check2")
        checkEquals(mcmcout@M, mcmc@M)
	checkTrue(is.list(mcmcout@par), "check3")
	checkTrue("lambda" %in% names(mcmcout@par), "check4")
        checkEquals(ncol(mcmcout@par$lambda), model@K)
        checkEquals(nrow(mcmcout@par$lambda), mcmc@M)
	checkTrue(!any(is.na(mcmcout@par$lambda)), "check5")
	for(k in 1:model@K) {
		checkTrue(all(abs(mcmcout@par$lambda[,k] - 
			mean(mcmcout@par$lambda[,k])) > 0), 
			paste("check", 5 + k, sep =""))
	}
	checkTrue(is.list(mcmcout@log), "check9")
	checkTrue("mixlik" %in% names(mcmcout@log),"check10")
	checkTrue("mixprior" %in% names(mcmcout@log), "check11")
	checkTrue(!("cdpost" %in% names(mcmcout@log)), "check12")
	checkEquals(nrow(mcmcout@log$mixlik), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixlik), 1)
	checkEquals(nrow(mcmcout@log$mixprior), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixprior), 1)
	checkTrue(!any(is.na(mcmcout@log$mixlik)), "check13")
	checkTrue(!any(is.na(mcmcout@log$mixprior)), "check14")
	checkTrue(all(abs(mcmcout@log$mixlik - mean(mcmcout@log$mixlik)) 
			> 0), "check15")
	checkTrue(all(abs(mcmcout@log$mixprior - mean(mcmcout@log$mixprior)) 
			> 0), "check16")
	checkTrue(!.hasSlot(mcmcout, "weight"), "check17")
	checkTrue(is.list(mcmcout@hyper), "check18")
	checkTrue("b" %in% names(mcmcout@hyper), "check19")
	checkEquals(nrow(mcmcout@hyper$b), mcmc@M)
	checkEquals(ncol(mcmcout@hyper$b), 1)
	checkTrue(!any(is.na(mcmcout@hyper$b)), "check20")
	checkTrue(all(abs(mcmcout@hyper$b - mean(mcmcout@hyper$b)) 
			> 0), "check21")
	checkTrue(!.hasSlot(mcmcout, "post"), "check22")
	checkTrue(!.hasSlot(mcmcout, "entropy"), "check23")
	checkTrue(!.hasSlot(mcmcout, "ST"), "check24")
	checkTrue(!.hasSlot(mcmcout, "S"), "check25")
	checkTrue(!.hasSlot(mcmcout, "NK"), "check26")
	checkTrue(!.hasSlot(mcmcout, "clust"), "check27")
	checkTrue(identical(mcmcout@model, model), "check28")
	checkTrue(identical(mcmcout@prior, prior), "check29")	
}

"test.FIX.POST.STARTPAR.1"  <- function() {
        ## Setting:
        ##      indicfix: TRUE
 	set.seed(0)       
        data <- .setUp.data.startpar(withInd = FALSE)
        model <- .setUp.model.startpar()
	setK(model) <- 1
        setIndicFix(model) <- TRUE
	prior <- prior(hier = FALSE)
        prior <- priordefine(data, model, varargin = prior)
        mcmc <- .setUp.mcmc.startpar()
	(data ~ model ~ mcmc) %=% mcmcstart(data,model,mcmc)
	setStorepost(mcmc) <- TRUE
        mcmcout <- mixturemcmc(data,model,prior,mcmc)
        ## Test cases ##
	checkTrue(identical(class(mcmcout)[1], "mcmcoutputfixpost"), "check1")
        checkTrue(mcmcout@ranperm == FALSE, "check2")
        checkEquals(mcmcout@M, mcmc@M)
	checkTrue(is.list(mcmcout@par), "check3")
	checkTrue("lambda" %in% names(mcmcout@par), "check4")
        checkEquals(ncol(mcmcout@par$lambda), model@K)
        checkEquals(nrow(mcmcout@par$lambda), mcmc@M)
	checkTrue(!any(is.na(mcmcout@par$lambda)), "check5")
	for(k in 1:model@K) {
		checkTrue(all(abs(mcmcout@par$lambda[,k] - 
			mean(mcmcout@par$lambda[,k])) > 0), 
			paste("check", 5 + k, sep =""))
	}
	checkTrue(is.list(mcmcout@log), "check7")
	checkTrue("mixlik" %in% names(mcmcout@log),"check8")
	checkTrue("mixprior" %in% names(mcmcout@log), "check9")
	checkTrue(!("cdpost" %in% names(mcmcout@log)), "check10")
	checkEquals(nrow(mcmcout@log$mixlik), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixlik), 1)
	checkEquals(nrow(mcmcout@log$mixprior), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixprior), 1)
	checkTrue(!any(is.na(mcmcout@log$mixlik)), "check11")
	checkTrue(!any(is.na(mcmcout@log$mixprior)), "check12")
	checkTrue(all(abs(mcmcout@log$mixlik - mean(mcmcout@log$mixlik)) 
			> 0), "check13")
	checkTrue(all(abs(mcmcout@log$mixprior - mean(mcmcout@log$mixprior)) 
			> 0), "check14")
	checkTrue(!.hasSlot(mcmcout, "weight"), "check15")
	checkTrue(!.hasSlot(mcmcout, "hyper"), "check16")
	checkTrue(is.list(mcmcout@post), "check17")
	checkTrue("par" %in% names(mcmcout@post), "check18")
	checkTrue(is.list(mcmcout@post$par), "check19")
	checkTrue("a" %in% names(mcmcout@post$par), "check20")
	checkTrue("b" %in% names(mcmcout@post$par), "check21")
	checkEquals(nrow(mcmcout@post$par$a), mcmc@M)
	checkEquals(nrow(mcmcout@post$par$b), mcmc@M)
	checkEquals(ncol(mcmcout@post$par$a), model@K)
	checkEquals(ncol(mcmcout@post$par$b), model@K)
	for(k in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@post$par$a[, k])), 
			paste("check", 21 + k, sep = ""))
		checkTrue(!any(is.na(mcmcout@post$par$b[, k])), 
			paste("check", 22 + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$par$a[, k] - 
			mean(mcmcout@post$par$a[, k])) == 0),
			paste("check", 23 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$par$b[, k] - 
			mean(mcmcout@post$par$b[, k])) == 0),
			paste("check", 24 + m * (k - 1) + k, sep = ""))

	}
	checkTrue(!.hasSlot(mcmcout, "entropy"), "check26")
	checkTrue(!.hasSlot(mcmcout, "ST"), "check27")
	checkTrue(!.hasSlot(mcmcout, "S"), "check28")
	checkTrue(!.hasSlot(mcmcout, "NK"), "check29")
	checkTrue(!.hasSlot(mcmcout, "clust"), "check30")
	checkTrue(identical(mcmcout@model, model), "check31")
	checkTrue(identical(mcmcout@prior, prior), "check32")	
}

"test.FIX.POST.STARTPAR.2" <- function() {
        ## Setting:
        ##      indicfix: TRUE
       	set.seed(0)       
	data <- .setUp.data.startpar(withInd = TRUE)
        model <- .setUp.model.startpar()
        setIndicFix(model) <- TRUE
	prior <- prior(hier = FALSE)
        prior <- priordefine(data, model, varargin = prior)
        mcmc <- .setUp.mcmc.startpar()
	setStorepost(mcmc) <- TRUE
       	(data ~ model ~ mcmc) %=% mcmcstart(data,model,mcmc)
	mcmcout <- mixturemcmc(data,model,prior,mcmc)
        ## Test cases ##
	checkTrue(identical(class(mcmcout)[1], "mcmcoutputfixpost"), "check1")
        checkTrue(mcmcout@ranperm == FALSE, "check2")
        checkEquals(mcmcout@M, mcmc@M)
	checkTrue(is.list(mcmcout@par), "check3")
	checkTrue("lambda" %in% names(mcmcout@par), "check4")
        checkEquals(ncol(mcmcout@par$lambda), model@K)
        checkEquals(nrow(mcmcout@par$lambda), mcmc@M)
	checkTrue(!any(is.na(mcmcout@par$lambda)), "check5")
	for(k in 1:model@K) {
		checkTrue(all(abs(mcmcout@par$lambda[,k] - 
			mean(mcmcout@par$lambda[,k])) > 0), 
			paste("check", 5 + k, sep =""))
	}
	checkTrue(is.list(mcmcout@log), "check8")
	checkTrue("mixlik" %in% names(mcmcout@log),"check9")
	checkTrue("mixprior" %in% names(mcmcout@log), "check10")
	checkTrue(!("cdpost" %in% names(mcmcout@log)), "check11")
	checkEquals(nrow(mcmcout@log$mixlik), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixlik), 1)
	checkEquals(nrow(mcmcout@log$mixprior), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixprior), 1)
	checkTrue(!any(is.na(mcmcout@log$mixlik)), "check12")
	checkTrue(!any(is.na(mcmcout@log$mixprior)), "check13")
	checkTrue(all(abs(mcmcout@log$mixlik - mean(mcmcout@log$mixlik)) 
			> 0), "check14")
	checkTrue(all(abs(mcmcout@log$mixprior - mean(mcmcout@log$mixprior)) 
			> 0), "check15")
	checkTrue(!.hasSlot(mcmcout, "weight"), "check16")
	checkTrue(!.hasSlot(mcmcout, "hyper"), "check17")
	checkTrue(is.list(mcmcout@post), "check18")	
	checkTrue("par" %in% names(mcmcout@post), "check19")
	checkTrue(is.list(mcmcout@post$par), "check20")
	checkTrue("a" %in% names(mcmcout@post$par), "check21")
	checkTrue("b" %in% names(mcmcout@post$par), "check22")
	checkEquals(nrow(mcmcout@post$par$a), mcmc@M)
	checkEquals(nrow(mcmcout@post$par$b), mcmc@M)
	checkEquals(ncol(mcmcout@post$par$a), model@K)
	checkEquals(ncol(mcmcout@post$par$b), model@K)
	m <- 4
	for(k in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@post$par$a[, k])), 
			paste("check", 22 + m * (k - 1) + k, sep = ""))
		checkTrue(!any(is.na(mcmcout@post$par$b[, k])), 
			paste("check", 23 + m * (k - 1)+ k, sep = ""))
		checkTrue(all(abs(mcmcout@post$par$a[, k] - 
			mean(mcmcout@post$par$a[, k])) == 0),
			paste("check", 24 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$par$b[, k] - 
			mean(mcmcout@post$par$b[, k])) == 0),
			paste("check", 25 + m * (k - 1) + k, sep = ""))


	}
	checkTrue(!.hasSlot(mcmcout, "entropy"), "check31")
	checkTrue(!.hasSlot(mcmcout, "ST"), "check32")
	checkTrue(!.hasSlot(mcmcout, "S"), "check33")
	checkTrue(!.hasSlot(mcmcout, "NK"), "check34")
	checkTrue(!.hasSlot(mcmcout, "clust"), "check35")
	checkTrue(identical(mcmcout@model, model), "check36")
	checkTrue(identical(mcmcout@prior, prior), "check37")	
}

"test.FIX.POST.STARTPAR.3" <- function() {
        ## Setting:
        ##      indicfix: TRUE
       	set.seed(0)       
	data <- .setUp.data.startpar(withInd = TRUE)
        model <- .setUp.model.startpar()
	setK(model) <- 3
        setIndicFix(model) <- TRUE
	prior <- prior(hier = FALSE)
        prior <- priordefine(data, model, varargin = prior)
        mcmc <- .setUp.mcmc.startpar()
	setStorepost(mcmc) <- TRUE
 	(data ~ model ~ mcmc) %=% mcmcstart(data,model,mcmc)
      	 mcmcout <- mixturemcmc(data,model,prior,mcmc)
        ## Test cases ##
	checkTrue(identical(class(mcmcout)[1], "mcmcoutputfixpost"), "check1")
        checkTrue(mcmcout@ranperm == FALSE, "check2")
        checkEquals(mcmcout@M, mcmc@M)
	checkTrue(is.list(mcmcout@par), "check3")
	checkTrue("lambda" %in% names(mcmcout@par), "check4")
        checkEquals(ncol(mcmcout@par$lambda), model@K)
        checkEquals(nrow(mcmcout@par$lambda), mcmc@M)
	checkTrue(!any(is.na(mcmcout@par$lambda)), "check5")
	for(k in 1:model@K) {
		checkTrue(all(abs(mcmcout@par$lambda[,k] - 
			mean(mcmcout@par$lambda[,k])) > 0), 
			paste("check", 5 + k, sep =""))
	}
	checkTrue(is.list(mcmcout@log), "check9")
	checkTrue("mixlik" %in% names(mcmcout@log),"check10")
	checkTrue("mixprior" %in% names(mcmcout@log), "check11")
	checkTrue(!("cdpost" %in% names(mcmcout@log)), "check12")
	checkEquals(nrow(mcmcout@log$mixlik), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixlik), 1)
	checkEquals(nrow(mcmcout@log$mixprior), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixprior), 1)
	checkTrue(!any(is.na(mcmcout@log$mixlik)), "check13")
	checkTrue(!any(is.na(mcmcout@log$mixprior)), "check14")
	checkTrue(all(abs(mcmcout@log$mixlik - mean(mcmcout@log$mixlik)) 
			> 0), "check15")
	checkTrue(all(abs(mcmcout@log$mixprior - mean(mcmcout@log$mixprior)) 
			> 0), "check16")
	checkTrue(!.hasSlot(mcmcout, "weight"), "check17")
	checkTrue(!.hasSlot(mcmcout, "hyper"), "check18")
	checkTrue(is.list(mcmcout@post), "check19")
	checkTrue("par" %in% names(mcmcout@post), "check20")
	checkTrue(is.list(mcmcout@post$par), "check21")
	checkTrue("a" %in% names(mcmcout@post$par), "check22")
	checkTrue("b" %in% names(mcmcout@post$par), "check23")
	checkEquals(nrow(mcmcout@post$par$a), mcmc@M)
	checkEquals(nrow(mcmcout@post$par$b), mcmc@M)
	checkEquals(ncol(mcmcout@post$par$a), model@K)
	checkEquals(ncol(mcmcout@post$par$b), model@K)
	m <- 4
	for(k in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@post$par$a[, k])), 
			paste("check", 23 + m * (k - 1) + k, sep = ""))
		checkTrue(!any(is.na(mcmcout@post$par$b[, k])), 
			paste("check", 24 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$par$a[, k] - 
			mean(mcmcout@post$par$a[, k])) == 0),
			paste("check", 25 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$par$b[, k] - 
			mean(mcmcout@post$par$b[, k])) == 0),
			paste("check", 26 + m * (k - 1) + k, sep = ""))


	}
	checkTrue(!.hasSlot(mcmcout, "entropy"), "check36")
	checkTrue(!.hasSlot(mcmcout, "ST"), "check37")
	checkTrue(!.hasSlot(mcmcout, "S"), "check38")
	checkTrue(!.hasSlot(mcmcout, "NK"), "check39")
	checkTrue(!.hasSlot(mcmcout, "clust"), "check40")
	checkTrue(identical(mcmcout@model, model), "check41")
	checkTrue(identical(mcmcout@prior, prior), "check42")	
}

"test.FIX.HIER.POST.STARTPAR.1"  <- function() {
        ## Setting:
        ##      indicfix: TRUE
       	set.seed(0)       
	data <- .setUp.data.startpar(withInd = FALSE)
        model <- .setUp.model.startpar()
	setK(model) <- 1
        setIndicFix(model) <- TRUE
        prior <- priordefine(data, model)
        mcmc <- .setUp.mcmc.startpar()
	setStorepost(mcmc) <- TRUE
 	(data ~ model ~ mcmc) %=% mcmcstart(data,model,mcmc)
      	mcmcout <- mixturemcmc(data,model,prior,mcmc)
        ## Test cases ##
	checkTrue(identical(class(mcmcout)[1], "mcmcoutputfixhierpost"), "check1")
        checkTrue(mcmcout@ranperm == FALSE, "check2")
        checkEquals(mcmcout@M, mcmc@M)
	checkTrue(is.list(mcmcout@par), "check3")
	checkTrue("lambda" %in% names(mcmcout@par), "check4")
        checkEquals(ncol(mcmcout@par$lambda), model@K)
        checkEquals(nrow(mcmcout@par$lambda), mcmc@M)
	checkTrue(!any(is.na(mcmcout@par$lambda)), "check5")
	for(k in 1:model@K) {
		checkTrue(all(abs(mcmcout@par$lambda[,k] - 
			mean(mcmcout@par$lambda[,k])) > 0), 
			paste("check", 5 + k, sep =""))
	}
	checkTrue(is.list(mcmcout@log), "check7")
	checkTrue("mixlik" %in% names(mcmcout@log),"check8")
	checkTrue("mixprior" %in% names(mcmcout@log), "check9")
	checkTrue(!("cdpost" %in% names(mcmcout@log)), "check10")
	checkEquals(nrow(mcmcout@log$mixlik), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixlik), 1)
	checkEquals(nrow(mcmcout@log$mixprior), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixprior), 1)
	checkTrue(!any(is.na(mcmcout@log$mixlik)), "check11")
	checkTrue(!any(is.na(mcmcout@log$mixprior)), "check12")
	checkTrue(all(abs(mcmcout@log$mixlik - mean(mcmcout@log$mixlik)) 
			> 0), "check13")
	checkTrue(all(abs(mcmcout@log$mixprior - mean(mcmcout@log$mixprior)) 
			> 0), "check14")
	checkTrue(!.hasSlot(mcmcout, "weight"), "check15")
	checkTrue(is.list(mcmcout@hyper), "check17")
	checkTrue("b" %in% names(mcmcout@hyper), "check18")
	checkEquals(nrow(mcmcout@hyper$b), mcmc@M)
	checkEquals(ncol(mcmcout@hyper$b), 1)
	checkTrue(!any(is.na(mcmcout@hyper$b)), "check19")
	checkTrue(all(abs(mcmcout@hyper$b - mean(mcmcout@hyper$b)) 
			> 0), "check20")
	checkTrue(is.list(mcmcout@post), "check21")
	checkTrue("par" %in% names(mcmcout@post), "check22")
	checkTrue(is.list(mcmcout@post$par), "check23")
	checkTrue("a" %in% names(mcmcout@post$par), "check24")
	checkTrue("b" %in% names(mcmcout@post$par), "check25")
	checkEquals(nrow(mcmcout@post$par$a), mcmc@M)
	checkEquals(nrow(mcmcout@post$par$b), mcmc@M)
	checkEquals(ncol(mcmcout@post$par$a), model@K)
	checkEquals(ncol(mcmcout@post$par$b), model@K)
	m <- 4
	for(k in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@post$par$a[, k])), 
			paste("check", 25 + m * ( k - 1) + k, sep = ""))
		checkTrue(!any(is.na(mcmcout@post$par$b[, k])), 
			paste("check", 26 + m * (k - 1) + k, sep = ""))	
		checkTrue(all(abs(mcmcout@post$par$a[, k] - 
			mean(mcmcout@post$par$a[, k])) == 0),
			paste("check", 27 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$par$b[, k] -
			 mean(mcmcout@post$par$b[, k])) > 0), 
			paste("check", 28 + m * (k - 1) + k, sep = ""))
	}
	checkTrue(!.hasSlot(mcmcout, "entropy"), "check29")
	checkTrue(!.hasSlot(mcmcout, "ST"), "check30")
	checkTrue(!.hasSlot(mcmcout, "S"), "check31")
	checkTrue(!.hasSlot(mcmcout, "NK"), "check31")
	checkTrue(!.hasSlot(mcmcout, "clust"), "check32")
	checkTrue(identical(mcmcout@model, model), "check33")
	checkTrue(identical(mcmcout@prior, prior), "check34")	
}

"test.FIX.HIER.POST.STARTPAR.2" <- function() {
        ## Setting:
        ##      indicfix: TRUE
       	set.seed(0)       
	data <- .setUp.data.startpar(withInd = TRUE)
        model <- .setUp.model.startpar()
        setIndicFix(model) <- TRUE
        prior <- priordefine(data, model)
        mcmc <- .setUp.mcmc.startpar()
	setStorepost(mcmc) <- TRUE
       	(data ~ model ~ mcmc) %=% mcmcstart(data,model,mcmc)
	mcmcout <- mixturemcmc(data,model,prior,mcmc)
        ## Test cases ##
	checkTrue(identical(class(mcmcout)[1], "mcmcoutputfixhierpost"), "check1")
        checkTrue(mcmcout@ranperm == FALSE, "check2")
        checkEquals(mcmcout@M, mcmc@M)
	checkTrue(is.list(mcmcout@par), "check3")
	checkTrue("lambda" %in% names(mcmcout@par), "check4")
        checkEquals(ncol(mcmcout@par$lambda), model@K)
        checkEquals(nrow(mcmcout@par$lambda), mcmc@M)
	checkTrue(!any(is.na(mcmcout@par$lambda)), "check5")
	for(k in 1:model@K) {
		checkTrue(all(abs(mcmcout@par$lambda[,k] - 
			mean(mcmcout@par$lambda[,k])) > 0), 
			paste("check", 5 + k, sep =""))
	}
	checkTrue(is.list(mcmcout@log), "check8")
	checkTrue("mixlik" %in% names(mcmcout@log),"check9")
	checkTrue("mixprior" %in% names(mcmcout@log), "check10")
	checkTrue(!("cdpost" %in% names(mcmcout@log)), "check11")
	checkEquals(nrow(mcmcout@log$mixlik), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixlik), 1)
	checkEquals(nrow(mcmcout@log$mixprior), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixprior), 1)
	checkTrue(!any(is.na(mcmcout@log$mixlik)), "check12")
	checkTrue(!any(is.na(mcmcout@log$mixprior)), "check13")
	checkTrue(all(abs(mcmcout@log$mixlik - mean(mcmcout@log$mixlik)) 
			> 0), "check14")
	checkTrue(all(abs(mcmcout@log$mixprior - mean(mcmcout@log$mixprior)) 
			> 0), "check15")
	checkTrue(!.hasSlot(mcmcout, "weight"), "check16")
	checkTrue(is.list(mcmcout@hyper), "check17")
	checkTrue("b" %in% names(mcmcout@hyper), "check18")
	checkEquals(nrow(mcmcout@hyper$b), mcmc@M)
	checkEquals(ncol(mcmcout@hyper$b), 1)
	checkTrue(!any(is.na(mcmcout@hyper$b)), "check19")
	checkTrue(all(abs(mcmcout@hyper$b - mean(mcmcout@hyper$b)) 
			> 0), "check20")
	checkTrue(is.list(mcmcout@post), "check21")	
	checkTrue("par" %in% names(mcmcout@post), "check22")
	checkTrue(is.list(mcmcout@post$par), "check23")
	checkTrue("a" %in% names(mcmcout@post$par), "check24")
	checkTrue("b" %in% names(mcmcout@post$par), "check25")
	checkEquals(nrow(mcmcout@post$par$a), mcmc@M)
	checkEquals(nrow(mcmcout@post$par$b), mcmc@M)
	checkEquals(ncol(mcmcout@post$par$a), model@K)
	checkEquals(ncol(mcmcout@post$par$b), model@K)
	m <- 4
	for(k in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@post$par$a[, k])), 
			paste("check", 25 + m * (k - 1) + k, sep = ""))
		checkTrue(!any(is.na(mcmcout@post$par$b[, k])), 
			paste("check", 26 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$par$a[, k] - 
			mean(mcmcout@post$par$a[, k])) == 0),
			paste("check", 27 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$par$b[, k] -
			mean(mcmcout@post$par$b[, k])) > 0), 
			paste("check", 28 + m * (k - 1) + k, sep = ""))
	}
	checkTrue(!.hasSlot(mcmcout, "entropy"), "check34")
	checkTrue(!.hasSlot(mcmcout, "ST"), "check35")
	checkTrue(!.hasSlot(mcmcout, "S"), "check36")
	checkTrue(!.hasSlot(mcmcout, "NK"), "check37")
	checkTrue(!.hasSlot(mcmcout, "clust"), "check38")
	checkTrue(identical(mcmcout@model, model), "check39")
	checkTrue(identical(mcmcout@prior, prior), "check40")	
}

"test.FIX.HIER.POST.STARTPAR.3" <- function() {
        ## Setting:
        ##      indicfix: TRUE
       	set.seed(0)       
	data <- .setUp.data.startpar(withInd = TRUE)
        model <- .setUp.model.startpar()
	setK(model) <- 3
        setIndicFix(model) <- TRUE
        prior <- priordefine(data, model)
        mcmc <- .setUp.mcmc.startpar()
	setStorepost(mcmc) <- TRUE
       	(data ~ model ~ mcmc) %=% mcmcstart(data,model,mcmc)
	mcmcout <- mixturemcmc(data,model,prior,mcmc)
        ## Test cases ##
	checkTrue(identical(class(mcmcout)[1], "mcmcoutputfixhierpost"), "check1")
        checkTrue(mcmcout@ranperm == FALSE, "check2")
        checkEquals(mcmcout@M, mcmc@M)
	checkTrue(is.list(mcmcout@par), "check3")
	checkTrue("lambda" %in% names(mcmcout@par), "check4")
        checkEquals(ncol(mcmcout@par$lambda), model@K)
        checkEquals(nrow(mcmcout@par$lambda), mcmc@M)
	checkTrue(!any(is.na(mcmcout@par$lambda)), "check5")
	for(k in 1:model@K) {
		checkTrue(all(abs(mcmcout@par$lambda[,k] - 
			mean(mcmcout@par$lambda[,k])) > 0), 
			paste("check", 5 + k, sep =""))
	}
	checkTrue(is.list(mcmcout@log), "check9")
	checkTrue("mixlik" %in% names(mcmcout@log),"check10")
	checkTrue("mixprior" %in% names(mcmcout@log), "check11")
	checkTrue(!("cdpost" %in% names(mcmcout@log)), "check12")
	checkEquals(nrow(mcmcout@log$mixlik), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixlik), 1)
	checkEquals(nrow(mcmcout@log$mixprior), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixprior), 1)
	checkTrue(!any(is.na(mcmcout@log$mixlik)), "check13")
	checkTrue(!any(is.na(mcmcout@log$mixprior)), "check14")
	checkTrue(all(abs(mcmcout@log$mixlik - mean(mcmcout@log$mixlik)) 
			> 0), "check15")
	checkTrue(all(abs(mcmcout@log$mixprior - mean(mcmcout@log$mixprior)) 
			> 0), "check16")
	checkTrue(!.hasSlot(mcmcout, "weight"), "check17")
	checkTrue(is.list(mcmcout@hyper), "check18")
	checkTrue("b" %in% names(mcmcout@hyper), "check19")
	checkEquals(nrow(mcmcout@hyper$b), mcmc@M)
	checkEquals(ncol(mcmcout@hyper$b), 1)
	checkTrue(!any(is.na(mcmcout@hyper$b)), "check20")
	checkTrue(all(abs(mcmcout@hyper$b - mean(mcmcout@hyper$b)) 
			> 0), "check21")
	checkTrue(is.list(mcmcout@post), "check22")
	checkTrue("par" %in% names(mcmcout@post), "check23")
	checkTrue(is.list(mcmcout@post$par), "check24")
	checkTrue("a" %in% names(mcmcout@post$par), "check25")
	checkTrue("b" %in% names(mcmcout@post$par), "check26")
	checkEquals(nrow(mcmcout@post$par$a), mcmc@M)
	checkEquals(nrow(mcmcout@post$par$b), mcmc@M)
	checkEquals(ncol(mcmcout@post$par$a), model@K)
	checkEquals(ncol(mcmcout@post$par$b), model@K)
	m <- 4
	for(k in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@post$par$a[, k])), 
			paste("check", 26 + m * (k - 1) + k, sep = ""))
		checkTrue(!any(is.na(mcmcout@post$par$b[, k])), 
			paste("check", 27 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$par$a[, k] - 
			mean(mcmcout@post$par$a[, k])) == 0),
			paste("check", 28 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$par$b[, k] -
			 mean(mcmcout@post$par$b[, k])) > 0), 
			paste("check", 29 + m * (k - 1) + k, sep = ""))
	}

	checkTrue(!.hasSlot(mcmcout, "entropy"), "check39")
	checkTrue(!.hasSlot(mcmcout, "ST"), "check40")
	checkTrue(!.hasSlot(mcmcout, "S"), "check41")
	checkTrue(!.hasSlot(mcmcout, "NK"), "check42")
	checkTrue(!.hasSlot(mcmcout, "clust"), "check43")
	checkTrue(identical(mcmcout@model, model), "check44")
	checkTrue(identical(mcmcout@prior, prior), "check45")	
}

"test.IND.STARTPAR.1"  <- function() {
        ## Setting:
        ##      indicfix: FALSE
       	set.seed(0)       
	data <- .setUp.data.startpar(withInd = FALSE)
        model <- .setUp.model.startpar()
	setK(model) <- 1
        setIndicFix(model) <- FALSE
	prior <- prior(hier = FALSE)
        prior <- priordefine(data, model, varargin = prior)
        mcmc <- .setUp.mcmc.startpar()
       	(data ~ model ~ mcmc) %=% mcmcstart(data,model,mcmc)
	mcmcout <- mixturemcmc(data,model,prior,mcmc)
        ## Test cases ##
	checkTrue(identical(class(mcmcout)[1], "mcmcoutputfix"), "check1")
        checkTrue(mcmcout@ranperm == FALSE, "check2")
        checkEquals(mcmcout@M, mcmc@M)
	checkTrue(is.list(mcmcout@par), "check3")
	checkTrue("lambda" %in% names(mcmcout@par), "check4")
        checkEquals(ncol(mcmcout@par$lambda), model@K)
        checkEquals(nrow(mcmcout@par$lambda), mcmc@M)
	checkTrue(!any(is.na(mcmcout@par$lambda)), "check5")
	for(k in 1:model@K) {
		checkTrue(all(abs(mcmcout@par$lambda[,k] - 
			mean(mcmcout@par$lambda[,k])) > 0), 
			paste("check", 5 + k, sep =""))
	}
	checkTrue(is.list(mcmcout@log), "check7")
	checkTrue("mixlik" %in% names(mcmcout@log),"check8")
	checkTrue("mixprior" %in% names(mcmcout@log), "check9")
	checkTrue(!("cdpost" %in% names(mcmcout@log)), "check10")
	checkEquals(nrow(mcmcout@log$mixlik), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixlik), 1)
	checkEquals(nrow(mcmcout@log$mixprior), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixprior), 1)
	checkTrue(!any(is.na(mcmcout@log$mixlik)), "check11")
	checkTrue(!any(is.na(mcmcout@log$mixprior)), "check12")
	checkTrue(all(abs(mcmcout@log$mixlik - mean(mcmcout@log$mixlik)) 
			> 0), "check13")
	checkTrue(all(abs(mcmcout@log$mixprior - mean(mcmcout@log$mixprior)) 
			> 0), "check14")
	checkTrue(!.hasSlot(mcmcout, "weight"), "check15")
	checkTrue(!.hasSlot(mcmcout, "hyper"), "check16")
	checkTrue(!.hasSlot(mcmcout, "post"), "check17")
	checkTrue(!.hasSlot(mcmcout, "entropy"), "check18")
	checkTrue(!.hasSlot(mcmcout, "ST"), "check19")
	checkTrue(!.hasSlot(mcmcout, "S"), "check20")
	checkTrue(!.hasSlot(mcmcout, "NK"), "check21")
	checkTrue(!.hasSlot(mcmcout, "clust"), "check22")
	checkTrue(identical(mcmcout@model, model), "check23")
	checkTrue(identical(mcmcout@prior, prior), "check24")	
}

"test.IND.STARTPAR.2" <- function() {
        ## Setting:
        ##      indicfix: FALSE
	set.seed(0)
        data <- .setUp.data.startpar(withInd = TRUE)
        model <- .setUp.model.startpar()
        setIndicFix(model) <- FALSE
	prior <- prior(hier = FALSE)
        prior <- priordefine(data, model, varargin = prior)
        mcmc <- .setUp.mcmc.startpar()
 	(data ~ model ~ mcmc) %=% mcmcstart(data,model,mcmc)
      	mcmcout <- mixturemcmc(data,model,prior,mcmc)
        ## Test cases ##
	checkTrue(identical(class(mcmcout)[1], "mcmcoutputbase"), "check1")
        checkTrue(mcmcout@ranperm == FALSE, "check2")
        checkEquals(mcmcout@M, mcmc@M)
	checkTrue(is.list(mcmcout@par), "check3")
	checkTrue("lambda" %in% names(mcmcout@par), "check4")
        checkEquals(ncol(mcmcout@par$lambda), model@K)
        checkEquals(nrow(mcmcout@par$lambda), mcmc@M)
	checkTrue(!any(is.na(mcmcout@par$lambda)), "check5")
	for(k in 1:model@K) {
		checkTrue(all(abs(mcmcout@par$lambda[,k] - 
			mean(mcmcout@par$lambda[,k])) > 0), 
			paste("check", 5 + k, sep =""))
	}
	checkTrue(is.list(mcmcout@log), "check8")
	checkTrue("mixlik" %in% names(mcmcout@log),"check9")
	checkTrue("mixprior" %in% names(mcmcout@log), "check10")
	checkTrue("cdpost" %in% names(mcmcout@log), "check11")
	checkEquals(nrow(mcmcout@log$mixlik), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixlik), 1)
	checkEquals(nrow(mcmcout@log$mixprior), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixprior), 1)
	checkEquals(nrow(mcmcout@log$cdpost), mcmc@M)
	checkEquals(ncol(mcmcout@log$cdpost), 1)
	checkTrue(!any(is.na(mcmcout@log$mixlik)), "check12")
	checkTrue(!any(is.na(mcmcout@log$mixprior)), "check13")
	checkTrue(!any(is.na(mcmcout@log$cdpost)), "check14")
	checkTrue(all(abs(mcmcout@log$mixlik - mean(mcmcout@log$mixlik)) 
			> 0), "check15")
	checkTrue(all(abs(mcmcout@log$mixprior - mean(mcmcout@log$mixprior)) 
			> 0), "check16")
	checkTrue(all(abs(mcmcout@log$cdpost - mean(mcmcout@log$cdpost))
			> 0), "check17")
	checkEquals(nrow(mcmcout@weight), mcmc@M)
	checkEquals(ncol(mcmcout@weight), model@K)
	m <- 0
	for(k in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@weight[, k])), 
			paste("check", 17 + m + k, sep = ""))
		checkTrue(all(abs(mcmcout@weight[, k] 
			- mean(mcmcout@weight[, k])) > 0),
			paste("check", 18 + m + k, sep = ""))
		m <- m * k
	}
	checkTrue(!.hasSlot(mcmcout, "hyper"), "check22")
	checkTrue(!.hasSlot(mcmcout, "post"), "check23")
	checkEquals(nrow(mcmcout@entropy), mcmc@M)
	checkEquals(ncol(mcmcout@entropy), 1)
	checkTrue(!any(is.na(mcmcout@entropy)), "check24")
	checkTrue(all(abs(mcmcout@entropy 
		- mean(mcmcout@entropy)) > 0), "check25")
	checkEquals(nrow(mcmcout@ST), mcmc@M)
	checkEquals(ncol(mcmcout@ST), 1)	
	checkTrue(!any(is.na(mcmcout@ST)), "check26")
	checkTrue(all(abs(mcmcout@ST - mean(mcmcout@ST)) 
			> 0), "check27")
	checkEquals(nrow(mcmcout@S), data@N)
	checkEquals(ncol(mcmcout@S), mcmc@storeS)
	checkTrue(!any(is.na(mcmcout@S[, 1])), "check28")
	checkTrue(all(abs(mcmcout@S[, 1] 
		- mean(mcmcout@S[, 1])) == 0), "check29")
	m <- 2
	for(s in 2:mcmc@storeS) {
		checkTrue(!any(is.na(mcmcout@S[, s])), 
			paste("check", 29 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@S[, s] 
			- mean(mcmcout@S[, s])) > 0),
			paste("check", 30 + m * (k - 1) + k, sep = ""))
	}
	checkEquals(nrow(mcmcout@NK), mcmc@M)
	checkEquals(ncol(mcmcout@NK), model@K)
	m <- 2
	for(n in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@NK[, n])),
			paste("check", 31 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@NK[, n]
			- mean(mcmcout@NK[, n])) > 0),
			paste("check", 32 + m * (k - 1) + k, sep = ""))
	}
	checkEquals(nrow(mcmcout@clust), data@N)
	checkEquals(ncol(mcmcout@clust), 1)
	checkTrue(!any(is.na(mcmcout@clust)), "check35")
	checkTrue(all(abs(mcmcout@clust - mean(mcmcout@clust)) 
			> 0), "check36")
	checkTrue(identical(mcmcout@model, model), "check37")
	checkTrue(identical(mcmcout@prior, prior), "check38")	
}

"test.IND.STARTPAR.3" <- function() {
        ## Setting:
        ##      indicfix: FALSE
	set.seed(0)
        data <- .setUp.data.startpar(withInd = TRUE)
        model <- .setUp.model.startpar()
	setK(model) <- 3
        setIndicFix(model) <- FALSE
	prior <- prior(hier = FALSE)
        prior <- priordefine(data, model, varargin = prior)
        mcmc <- .setUp.mcmc.startpar()
       	(data ~ model ~ mcmc) %=% mcmcstart(data,model,mcmc)
	mcmcout <- mixturemcmc(data,model,prior,mcmc)
        ## Test cases ##
	checkTrue(identical(class(mcmcout)[1], "mcmcoutputbase"), "check1")
        checkTrue(mcmcout@ranperm == FALSE, "check2")
        checkEquals(mcmcout@M, mcmc@M)
	checkTrue(is.list(mcmcout@par), "check3")
	checkTrue("lambda" %in% names(mcmcout@par), "check4")
        checkEquals(ncol(mcmcout@par$lambda), model@K)
        checkEquals(nrow(mcmcout@par$lambda), mcmc@M)
	checkTrue(!any(is.na(mcmcout@par$lambda)), "check5")
	for(k in 1:model@K) {
		checkTrue(all(abs(mcmcout@par$lambda[,k] - 
			mean(mcmcout@par$lambda[,k])) > 0), 
			paste("check", 5 + k, sep =""))
	}
	checkTrue(is.list(mcmcout@log), "check9")
	checkTrue("mixlik" %in% names(mcmcout@log),"check10")
	checkTrue("mixprior" %in% names(mcmcout@log), "check11")
	checkTrue("cdpost" %in% names(mcmcout@log), "check12")
	checkEquals(nrow(mcmcout@log$mixlik), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixlik), 1)
	checkEquals(nrow(mcmcout@log$mixprior), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixprior), 1)
	checkEquals(nrow(mcmcout@log$cdpost), mcmc@M)
	checkEquals(ncol(mcmcout@log$cdpost), 1)
	checkTrue(!any(is.na(mcmcout@log$mixlik)), "check13")
	checkTrue(!any(is.na(mcmcout@log$mixprior)), "check14")
	checkTrue(!any(is.na(mcmcout@log$cdpost)), "check15")
	checkTrue(all(abs(mcmcout@log$mixlik - mean(mcmcout@log$mixlik)) 
			> 0), "check16")
	checkTrue(all(abs(mcmcout@log$mixprior - mean(mcmcout@log$mixprior)) 
			> 0), "check17")
	checkTrue(all(abs(mcmcout@log$cdpost - mean(mcmcout@log$cdpost))
			> 0), "check18")
	checkEquals(nrow(mcmcout@weight), mcmc@M)
	checkEquals(ncol(mcmcout@weight), model@K)
	m <- 2
	for(k in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@weight[, k])), 
			paste("check", 18 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@weight[, k] 
			- mean(mcmcout@weight[, k])) > 0),
			paste("check", 19 + m * (k - 1)  + k, sep = ""))
	}
	checkTrue(!.hasSlot(mcmcout, "hyper"), "check25")
	checkTrue(!.hasSlot(mcmcout, "post"), "check26")
	checkEquals(nrow(mcmcout@entropy), mcmc@M)
	checkEquals(ncol(mcmcout@entropy), 1)
	checkTrue(!any(is.na(mcmcout@entropy)), "check27")
	checkTrue(all(abs(mcmcout@entropy 
		- mean(mcmcout@entropy)) > 0), "check28")
	checkEquals(nrow(mcmcout@ST), mcmc@M)
	checkEquals(ncol(mcmcout@ST), 1)	
	checkTrue(!any(is.na(mcmcout@ST)), "check29")
	checkTrue(!all(abs(mcmcout@ST - mean(mcmcout@ST)) 
			== 0), "check30")
	checkEquals(nrow(mcmcout@S), data@N)
	checkEquals(ncol(mcmcout@S), mcmc@storeS)
	checkTrue(!any(is.na(mcmcout@S[, 1])), "check28")
	checkTrue(all(abs(mcmcout@S[, 1] 
		- mean(mcmcout@S[, 1])) == 0), "check29")
	m <- 2
	for(s in 2:mcmc@storeS) {
		checkTrue(!any(is.na(mcmcout@S[, s])), 
			paste("check", 30 + m + k, sep = ""))
		checkTrue(all(abs(mcmcout@S[, s] 
			- mean(mcmcout@S[, s])) > 0),
			paste("check", 31 + m + k, sep = ""))
		m <- m * k
	}
	checkEquals(nrow(mcmcout@NK), mcmc@M)
	checkEquals(ncol(mcmcout@NK), model@K)
	m <- 2
	for(n in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@NK[, n])),
			paste("check", 36 + m + k))
		checkTrue(all(abs(mcmcout@NK[, n]
			- mean(mcmcout@NK[, n])) > 0),
			paste("check", 37 + m + k, sep = ""))
	}
	checkEquals(nrow(mcmcout@clust), data@N)
	checkEquals(ncol(mcmcout@clust), 1)
	checkTrue(!any(is.na(mcmcout@clust)), "check43")
	checkTrue(all(abs(mcmcout@clust - mean(mcmcout@clust)) 
			> 0), "check44")
	checkTrue(identical(mcmcout@model, model), "check45")
	checkTrue(identical(mcmcout@prior, prior), "check46")	
}

"test.IND.HIER.STARTPAR.1"  <- function() {
        ## Setting:
        ##      indicfix: FALSE
      	set.seed(0)       
	data <- .setUp.data.startpar(withInd = TRUE)
        model <- .setUp.model.startpar()
	setK(model) <- 1
        setIndicFix(model) <- FALSE
        prior <- priordefine(data, model)
        mcmc <- .setUp.mcmc.startpar()
       	(data ~ model ~ mcmc) %=% mcmcstart(data,model,mcmc)
	mcmcout <- mixturemcmc(data,model,prior,mcmc)
        ## Test cases ##
	checkTrue(identical(class(mcmcout)[1], "mcmcoutputfixhier"), "check1")
        checkTrue(mcmcout@ranperm == FALSE, "check2")
        checkEquals(mcmcout@M, mcmc@M)
	checkTrue(is.list(mcmcout@par), "check3")
	checkTrue("lambda" %in% names(mcmcout@par), "check4")
        checkEquals(ncol(mcmcout@par$lambda), model@K)
        checkEquals(nrow(mcmcout@par$lambda), mcmc@M)
	checkTrue(!any(is.na(mcmcout@par$lambda)), "check5")
	for(k in 1:model@K) {
		checkTrue(all(abs(mcmcout@par$lambda[,k] - 
			mean(mcmcout@par$lambda[,k])) > 0), 
			paste("check", 5 + k, sep =""))
	}
	checkTrue(is.list(mcmcout@log), "check7")
	checkTrue("mixlik" %in% names(mcmcout@log),"check8")
	checkTrue("mixprior" %in% names(mcmcout@log), "check9")
	checkTrue(!("cdpost" %in% names(mcmcout@log)), "check10")
	checkEquals(nrow(mcmcout@log$mixlik), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixlik), 1)
	checkEquals(nrow(mcmcout@log$mixprior), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixprior), 1)
	checkTrue(!any(is.na(mcmcout@log$mixlik)), "check11")
	checkTrue(!any(is.na(mcmcout@log$mixprior)), "check12")
	checkTrue(all(abs(mcmcout@log$mixlik - mean(mcmcout@log$mixlik)) 
			> 0), "check13")
	checkTrue(all(abs(mcmcout@log$mixprior - mean(mcmcout@log$mixprior)) 
			> 0), "check14")
	checkTrue(!.hasSlot(mcmcout, "weight"), "check15")
	checkTrue(is.list(mcmcout@hyper), "check16")
	checkTrue("b" %in% names(mcmcout@hyper), "check17")
	checkEquals(nrow(mcmcout@hyper$b), mcmc@M)
	checkEquals(ncol(mcmcout@hyper$b), 1)
	checkTrue(!any(is.na(mcmcout@hyper$b)), "check18")
	checkTrue(all(abs(mcmcout@hyper$b - mean(mcmcout@hyper$b)) 
		> 0), "check19")
	checkTrue(!.hasSlot(mcmcout, "post"), "check20")
	checkTrue(!.hasSlot(mcmcout, "entropy"), "check21")
	checkTrue(!.hasSlot(mcmcout, "ST"), "check22")
	checkTrue(!.hasSlot(mcmcout, "S"), "check23")
	checkTrue(!.hasSlot(mcmcout, "NK"), "check24")
	checkTrue(!.hasSlot(mcmcout, "clust"), "check25")
	checkTrue(identical(mcmcout@model, model), "check26")
	checkTrue(identical(mcmcout@prior, prior), "check27")	
}

"test.IND.HIER.STARTPAR.2" <- function() {
        ## Setting:
        ##      indicfix: FALSE
	set.seed(0)       
	data <- .setUp.data.startpar(withInd = TRUE)
        model <- .setUp.model.startpar()
        setIndicFix(model) <- FALSE
        prior <- priordefine(data, model)
        mcmc <- .setUp.mcmc.startpar()
       	(data ~ model ~ mcmc) %=% mcmcstart(data,model,mcmc)
	mcmcout <- mixturemcmc(data,model,prior,mcmc)
        ## Test cases ##
	checkTrue(identical(class(mcmcout)[1], "mcmcoutputhier"), "check1")
        checkTrue(mcmcout@ranperm == FALSE, "check2")
        checkEquals(mcmcout@M, mcmc@M)
	checkTrue(is.list(mcmcout@par), "check3")
	checkTrue("lambda" %in% names(mcmcout@par), "check4")
        checkEquals(ncol(mcmcout@par$lambda), model@K)
        checkEquals(nrow(mcmcout@par$lambda), mcmc@M)
	checkTrue(!any(is.na(mcmcout@par$lambda)), "check5")
	for(k in 1:model@K) {
		checkTrue(all(abs(mcmcout@par$lambda[,k] - 
			mean(mcmcout@par$lambda[,k])) > 0), 
			paste("check", 5 + k, sep =""))
	}
	checkTrue(is.list(mcmcout@log), "check8")
	checkTrue("mixlik" %in% names(mcmcout@log),"check9")
	checkTrue("mixprior" %in% names(mcmcout@log), "check10")
	checkTrue("cdpost" %in% names(mcmcout@log), "check11")
	checkEquals(nrow(mcmcout@log$mixlik), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixlik), 1)
	checkEquals(nrow(mcmcout@log$mixprior), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixprior), 1)
	checkEquals(nrow(mcmcout@log$cdpost), mcmc@M)
	checkEquals(ncol(mcmcout@log$cdpost), 1)
	checkTrue(!any(is.na(mcmcout@log$mixlik)), "check12")
	checkTrue(!any(is.na(mcmcout@log$mixprior)), "check13")
	checkTrue(!any(is.na(mcmcout@log$cdpost)), "check14")
	checkTrue(all(abs(mcmcout@log$mixlik - mean(mcmcout@log$mixlik)) 
			> 0), "check15")
	checkTrue(all(abs(mcmcout@log$mixprior - mean(mcmcout@log$mixprior)) 
			> 0), "check16")
	checkTrue(all(abs(mcmcout@log$cdpost - mean(mcmcout@log$cdpost))
			> 0), "check17")
	checkEquals(nrow(mcmcout@weight), mcmc@M)
	checkEquals(ncol(mcmcout@weight), model@K)
	m <- 2
	for(k in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@weight[, k])), 
			paste("check", 17 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@weight[, k] 
			- mean(mcmcout@weight[, k])) > 0),
			paste("check", 18 + m * (k - 1) + k, sep = ""))
	}
	checkTrue(is.list(mcmcout@hyper), "check22")
	checkTrue("b" %in% names(mcmcout@hyper), "check23")
	checkEquals(nrow(mcmcout@hyper$b), mcmc@M)
	checkEquals(ncol(mcmcout@hyper$b), 1)
	checkTrue(!any(is.na(mcmcout@hyper$b)), "check24")
	checkTrue(all(abs(mcmcout@hyper$b - mean(mcmcout@hyper$b)) 
			> 0), "check25")
	checkTrue(!.hasSlot(mcmcout, "post"), "check26")
	checkEquals(nrow(mcmcout@entropy), mcmc@M)
	checkEquals(ncol(mcmcout@entropy), 1)
	checkTrue(!any(is.na(mcmcout@entropy)), "check27")
	checkTrue(all(abs(mcmcout@entropy 
		- mean(mcmcout@entropy)) > 0), "check28")
	checkEquals(nrow(mcmcout@ST), mcmc@M)
	checkEquals(ncol(mcmcout@ST), 1)	
	checkTrue(!any(is.na(mcmcout@ST)), "check29")
	checkTrue(all(abs(mcmcout@ST - mean(mcmcout@ST)) 
			> 0), "check30")
	checkEquals(nrow(mcmcout@S), data@N)
	checkEquals(ncol(mcmcout@S), mcmc@storeS)
	checkTrue(!any(is.na(mcmcout@S[, 1])), "check28")
	checkTrue(all(abs(mcmcout@S[, 1] 
		- mean(mcmcout@S[, 1])) == 0), "check29")
	m <- 2
	for(s in 2:mcmc@storeS) {
		checkTrue(!any(is.na(mcmcout@S[, s])), 
			paste("check", 30 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@S[, s] 
			- mean(mcmcout@S[, s])) > 0),
			paste("check", 31 + m * (k - 1) + k, sep = ""))
	}
	checkEquals(nrow(mcmcout@NK), mcmc@M)
	checkEquals(ncol(mcmcout@NK), model@K)
	for(n in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@NK[, n])),
			paste("check", 34 + m * (k - 1) + k))
		checkTrue(all(abs(mcmcout@NK[, n]
			- mean(mcmcout@NK[, n])) > 0),
			paste("check", 35 + m * (k - 1) + k, sep = ""))
	}
	checkEquals(nrow(mcmcout@clust), data@N)
	checkEquals(ncol(mcmcout@clust), 1)
	checkTrue(!any(is.na(mcmcout@clust)), "check38")
	checkTrue(all(abs(mcmcout@clust - mean(mcmcout@clust)) 
			> 0), "check39")
	checkTrue(identical(mcmcout@model, model), "check40")
	checkTrue(identical(mcmcout@prior, prior), "check41")	
}

"test.IND.HIER.STARTPAR.3" <- function() {
        ## Setting:
        ##      indicfix: FALSE
       	set.seed(0)
	data <- .setUp.data.startpar(withInd = TRUE)
        model <- .setUp.model.startpar()
        setIndicFix(model) <- FALSE
        prior <- priordefine(data, model)
        mcmc <- .setUp.mcmc.startpar()
       	(data ~ model ~ mcmc) %=% mcmcstart(data,model,mcmc)
	mcmcout <- mixturemcmc(data,model,prior,mcmc)
        ## Test cases ##
	checkTrue(identical(class(mcmcout)[1], "mcmcoutputhier"), "check1")
        checkTrue(mcmcout@ranperm == FALSE, "check2")
        checkEquals(mcmcout@M, mcmc@M)
	checkTrue(is.list(mcmcout@par), "check3")
	checkTrue("lambda" %in% names(mcmcout@par), "check4")
        checkEquals(ncol(mcmcout@par$lambda), model@K)
        checkEquals(nrow(mcmcout@par$lambda), mcmc@M)
	checkTrue(!any(is.na(mcmcout@par$lambda)), "check5")
	for(k in 1:model@K) {
		checkTrue(all(abs(mcmcout@par$lambda[,k] - 
			mean(mcmcout@par$lambda[,k])) > 0), 
			paste("check", 5 + k, sep =""))
	}
	checkTrue(is.list(mcmcout@log), "check9")
	checkTrue("mixlik" %in% names(mcmcout@log),"check10")
	checkTrue("mixprior" %in% names(mcmcout@log), "check11")
	checkTrue("cdpost" %in% names(mcmcout@log), "check12")
	checkEquals(nrow(mcmcout@log$mixlik), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixlik), 1)
	checkEquals(nrow(mcmcout@log$mixprior), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixprior), 1)
	checkEquals(nrow(mcmcout@log$cdpost), mcmc@M)
	checkEquals(ncol(mcmcout@log$cdpost), 1)
	checkTrue(!any(is.na(mcmcout@log$mixlik)), "check13")
	checkTrue(!any(is.na(mcmcout@log$mixprior)), "check14")
	checkTrue(!any(is.na(mcmcout@log$cdpost)), "check15")
	checkTrue(all(abs(mcmcout@log$mixlik - mean(mcmcout@log$mixlik)) 
			> 0), "check16")
	checkTrue(all(abs(mcmcout@log$mixprior - mean(mcmcout@log$mixprior)) 
			> 0), "check17")
	checkTrue(all(abs(mcmcout@log$cdpost - mean(mcmcout@log$cdpost))
			> 0), "check18")
	checkEquals(nrow(mcmcout@weight), mcmc@M)
	checkEquals(ncol(mcmcout@weight), model@K)
	m <- 2
	for(k in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@weight[, k])), 
			paste("check", 18 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@weight[, k] 
			- mean(mcmcout@weight[, k])) > 0),
			paste("check", 19 + m * (k - 1) + k, sep = ""))
	}
	checkTrue(is.list(mcmcout@hyper), "check25")
	checkTrue("b" %in% names(mcmcout@hyper), "check26")
	checkEquals(nrow(mcmcout@hyper$b), mcmc@M)
	checkEquals(ncol(mcmcout@hyper$b), 1)
	checkTrue(!any(is.na(mcmcout@hyper$b)), "check27")
	checkTrue(all(abs(mcmcout@hyper$b - mean(mcmcout@hyper$b)) 
			> 0), "check28")
	checkTrue(!.hasSlot(mcmcout, "post"), "check29")
	checkEquals(nrow(mcmcout@entropy), mcmc@M)
	checkEquals(ncol(mcmcout@entropy), 1)
	checkTrue(!any(is.na(mcmcout@entropy)), "check30")
	checkTrue(all(abs(mcmcout@entropy 
		- mean(mcmcout@entropy)) > 0), "check31")
	checkEquals(nrow(mcmcout@ST), mcmc@M)
	checkEquals(ncol(mcmcout@ST), 1)	
	checkTrue(!any(is.na(mcmcout@ST)), "check32")
	checkTrue(all(abs(mcmcout@ST - mean(mcmcout@ST)) 
			> 0), "check33")
	checkEquals(nrow(mcmcout@S), data@N)
	checkEquals(ncol(mcmcout@S), mcmc@storeS)
	checkTrue(!any(is.na(mcmcout@S[, 1])), "check28")
	checkTrue(all(abs(mcmcout@S[, 1] 
		- mean(mcmcout@S[, 1])) == 0), "check29")
	m <- 2
	for(s in 2:mcmc@storeS) {
		checkTrue(!any(is.na(mcmcout@S[, s])), 
			paste("check", 33 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@S[, s] 
			- mean(mcmcout@S[, s])) > 0),
			paste("check", 34 + m * (k - 1) + k, sep = ""))
	}
	checkEquals(nrow(mcmcout@NK), mcmc@M)
	checkEquals(ncol(mcmcout@NK), model@K)
	m <- 2
	for(n in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@NK[, n])),
			paste("check", 40 + m * (k - 1) + k))
		checkTrue(all(abs(mcmcout@NK[, n]
			- mean(mcmcout@NK[, n])) > 0),
			paste("check", 41 + m * (k - 1) + k, sep = ""))
	}
	checkEquals(nrow(mcmcout@clust), data@N)
	checkEquals(ncol(mcmcout@clust), 1)
	checkTrue(!any(is.na(mcmcout@clust)), "check47")
	checkTrue(all(abs(mcmcout@clust - mean(mcmcout@clust)) 
			> 0), "check48")
	checkTrue(identical(mcmcout@model, model), "check49")
	checkTrue(identical(mcmcout@prior, prior), "check50")	
}

"test.IND.POST.1"  <- function() {
        ## Setting:
        ##      indicfix: FALSE
       	set.seed(0)       
	data <- .setUp.data.startpar(withInd = FALSE)
        model <- .setUp.model.startpar()
	setK(model) <- 1
        setIndicFix(model) <- FALSE
	prior <- prior(hier = FALSE)
        prior <- priordefine(data, model, varargin = prior)
        mcmc <- .setUp.mcmc.startpar()
	setStorepost(mcmc) <- TRUE
 	(data ~ model ~ mcmc) %=% mcmcstart(data,model,mcmc)
      	mcmcout <- mixturemcmc(data,model,prior,mcmc)
        ## Test cases ##
	checkTrue(identical(class(mcmcout)[1], "mcmcoutputfixpost"), "check1")
        checkTrue(mcmcout@ranperm == FALSE, "check2")
        checkEquals(mcmcout@M, mcmc@M)
	checkTrue(is.list(mcmcout@par), "check3")
	checkTrue("lambda" %in% names(mcmcout@par), "check4")
        checkEquals(ncol(mcmcout@par$lambda), model@K)
        checkEquals(nrow(mcmcout@par$lambda), mcmc@M)
	checkTrue(!any(is.na(mcmcout@par$lambda)), "check5")
	for(k in 1:model@K) {
		checkTrue(all(abs(mcmcout@par$lambda[,k] - 
			mean(mcmcout@par$lambda[,k])) > 0), 
			paste("check", 5 + k, sep =""))
	}
	checkTrue(is.list(mcmcout@log), "check7")
	checkTrue("mixlik" %in% names(mcmcout@log),"check8")
	checkTrue("mixprior" %in% names(mcmcout@log), "check9")
	checkTrue(!("cdpost" %in% names(mcmcout@log)), "check10")
	checkEquals(nrow(mcmcout@log$mixlik), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixlik), 1)
	checkEquals(nrow(mcmcout@log$mixprior), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixprior), 1)
	checkTrue(!any(is.na(mcmcout@log$mixlik)), "check11")
	checkTrue(!any(is.na(mcmcout@log$mixprior)), "check12")
	checkTrue(all(abs(mcmcout@log$mixlik - mean(mcmcout@log$mixlik)) 
			> 0), "check13")
	checkTrue(all(abs(mcmcout@log$mixprior - mean(mcmcout@log$mixprior)) 
			> 0), "check14")
	checkTrue(!.hasSlot(mcmcout, "weight"), "check15")
	checkTrue(!.hasSlot(mcmcout, "hyper"), "check16")
	checkTrue(is.list(mcmcout@post), "check17")
	checkTrue("par" %in% names(mcmcout@post), "check18")
	checkTrue(is.list(mcmcout@post$par), "check19")
	checkTrue("a" %in% names(mcmcout@post$par), "check20")
	checkTrue("b" %in% names(mcmcout@post$par), "check21")
	checkEquals(nrow(mcmcout@post$par$a), mcmc@M)
	checkEquals(nrow(mcmcout@post$par$b), mcmc@M)
	checkEquals(ncol(mcmcout@post$par$a), model@K)
	checkEquals(ncol(mcmcout@post$par$b), model@K)
	for(k in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@post$par$a[, k])), 
			paste("check", 21 + k, sep = ""))
		checkTrue(!any(is.na(mcmcout@post$par$b[, k])), 
			paste("check", 22 + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$par$a[, k] - 
			mean(mcmcout@post$par$a[, k])) == 0), 
			paste("check", 23 + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$par$b[, k] -
			 mean(mcmcout@post$par$b[, k])) == 0), 
			paste("check", 24 + k, sep = ""))
	}
	checkTrue(!.hasSlot(mcmcout, "entropy"), "check26")
	checkTrue(!.hasSlot(mcmcout, "ST"), "check27")
	checkTrue(!.hasSlot(mcmcout, "S"), "check28")
	checkTrue(!.hasSlot(mcmcout, "NK"), "check29")
	checkTrue(!.hasSlot(mcmcout, "clust"), "check30")
	checkTrue(identical(mcmcout@model, model), "check31")
	checkTrue(identical(mcmcout@prior, prior), "check32")	
}

"test.IND.POST.STARTPAR.2" <- function() {
        ## Setting:
        ##      indicfix: FALSE
       	set.seed(0)       
	data <- .setUp.data.startpar(withInd = TRUE)
        model <- .setUp.model.startpar()
        setIndicFix(model) <- FALSE
	prior <- prior(hier = FALSE)
        prior <- priordefine(data, model, varargin = prior)
        mcmc <- .setUp.mcmc.startpar()
	setStorepost(mcmc) <- TRUE
       	(data ~ model ~ mcmc) %=% mcmcstart(data,model,mcmc)
	mcmcout <- mixturemcmc(data,model,prior,mcmc)
        ## Test cases ##
	checkTrue(identical(class(mcmcout)[1], "mcmcoutputpost"), "check1")
        checkTrue(mcmcout@ranperm == FALSE, "check2")
        checkEquals(mcmcout@M, mcmc@M)
	checkTrue(is.list(mcmcout@par), "check3")
	checkTrue("lambda" %in% names(mcmcout@par), "check4")
        checkEquals(ncol(mcmcout@par$lambda), model@K)
        checkEquals(nrow(mcmcout@par$lambda), mcmc@M)
	checkTrue(!any(is.na(mcmcout@par$lambda)), "check5")
	for(k in 1:model@K) {
		checkTrue(all(abs(mcmcout@par$lambda[,k] - 
			mean(mcmcout@par$lambda[,k])) > 0), 
			paste("check", 5 + k, sep =""))
	}
	checkTrue(is.list(mcmcout@log), "check8")
	checkTrue("mixlik" %in% names(mcmcout@log),"check9")
	checkTrue("mixprior" %in% names(mcmcout@log), "check10")
	checkTrue("cdpost" %in% names(mcmcout@log), "check11")
	checkEquals(nrow(mcmcout@log$mixlik), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixlik), 1)
	checkEquals(nrow(mcmcout@log$mixprior), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixprior), 1)
	checkEquals(nrow(mcmcout@log$cdpost), mcmc@M)
	checkEquals(ncol(mcmcout@log$cdpost), 1)
	checkTrue(!any(is.na(mcmcout@log$mixlik)), "check12")
	checkTrue(!any(is.na(mcmcout@log$mixprior)), "check13")
	checkTrue(!any(is.na(mcmcout@log$cdpost)), "check14")
	checkTrue(all(abs(mcmcout@log$mixlik - mean(mcmcout@log$mixlik)) 
			> 0), "check15")
	checkTrue(all(abs(mcmcout@log$mixprior - mean(mcmcout@log$mixprior)) 
			> 0), "check16")
	checkTrue(all(abs(mcmcout@log$cdpost - mean(mcmcout@log$cdpost))
			> 0), "check17")
	checkEquals(nrow(mcmcout@weight), mcmc@M)
	checkEquals(ncol(mcmcout@weight), model@K)
	m <- 2
	for(k in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@weight[, k])), 
			paste("check", 17 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@weight[, k] 
			- mean(mcmcout@weight[, k])) > 0),
			paste("check", 18 + m * (k - 1) + k, sep = ""))
	}
	checkTrue(!.hasSlot(mcmcout, "hyper"), "check22")
	checkTrue(is.list(mcmcout@post), "check18")	
	checkTrue("par" %in% names(mcmcout@post), "check19")
	checkTrue("weight" %in% names(mcmcout@post), "check20")
	checkEquals(nrow(mcmcout@post$weight), mcmc@M)
	checkEquals(ncol(mcmcout@post$weight), model@K)
	checkTrue(is.list(mcmcout@post$par), "check21")
	checkTrue("a" %in% names(mcmcout@post$par), "check22")
	checkTrue("b" %in% names(mcmcout@post$par), "check23")
	checkEquals(nrow(mcmcout@post$par$a), mcmc@M)
	checkEquals(nrow(mcmcout@post$par$b), mcmc@M)
	checkEquals(ncol(mcmcout@post$par$a), model@K)
	checkEquals(ncol(mcmcout@post$par$b), model@K)
	m <- 4
	for(k in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@post$par$a[, k])), 
			paste("check", 23 + m * (k - 1) + k, sep = ""))
		checkTrue(!any(is.na(mcmcout@post$par$b[, k])), 
			paste("check", 24 + m * (k - 1) + k, sep = ""))
		checkTrue(!any(is.na(mcmcout@post$weight[, k])),
			paste("check", 25 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$par$a[, k] - 
			mean(mcmcout@post$par$a[, k])) > 0), 
			paste("check", 26 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$par$b[, k] -
			 mean(mcmcout@post$par$b[, k])) > 0), 
			paste("check", 27 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$weight[, k] 
			- mean(mcmcout@post$weight[, k])) > 0),
			paste("check", 28 + m * (k - 1) + k, sep = ""))
	}
	checkEquals(nrow(mcmcout@entropy), mcmc@M)
	checkEquals(ncol(mcmcout@entropy), 1)
	checkTrue(!any(is.na(mcmcout@entropy)), "check35")
	checkTrue(all(abs(mcmcout@entropy 
		- mean(mcmcout@entropy)) > 0), "check36")
	checkEquals(nrow(mcmcout@ST), mcmc@M)
	checkEquals(ncol(mcmcout@ST), 1)	
	checkTrue(!any(is.na(mcmcout@ST)), "check37")
	checkTrue(all(abs(mcmcout@ST - mean(mcmcout@ST)) 
			> 0), "check38")
	checkEquals(nrow(mcmcout@S), data@N)
	checkEquals(ncol(mcmcout@S), mcmc@storeS)
	checkTrue(!any(is.na(mcmcout@S[, 1])), "check39")
	checkTrue(all(abs(mcmcout@S[, 1] 
		- mean(mcmcout@S[,1])) == 0), "check40")
	m <- 2
	for(s in 2:mcmc@storeS) {
		checkTrue(!any(is.na(mcmcout@S[, s])), 
			paste("check", 40 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@S[, s] 
			- mean(mcmcout@S[, s])) > 0),
			paste("check", 41 + m * (k - 1) + k, sep = ""))
	}
	checkEquals(nrow(mcmcout@NK), mcmc@M)
	checkEquals(ncol(mcmcout@NK), model@K)
	m <- 2
	for(n in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@NK[, n])),
			paste("check", 42 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@NK[, n]
			- mean(mcmcout@NK[, n])) > 0),
			paste("check", 43 + m * (k - 1) + k, sep = ""))
	}
	checkEquals(nrow(mcmcout@clust), data@N)
	checkEquals(ncol(mcmcout@clust), 1)
	checkTrue(!any(is.na(mcmcout@clust)), "check47")
	checkTrue(all(abs(mcmcout@clust - mean(mcmcout@clust)) 
			> 0), "check48")
	checkTrue(identical(mcmcout@model, model), "check49")
	checkTrue(identical(mcmcout@prior, prior), "check50")	
}

"test.IND.POST.STARTPAR.3" <- function() {
        ## Setting:
        ##      indicfix: FALSE
       	set.seed(0)       
	data <- .setUp.data.startpar(withInd = TRUE)
        model <- .setUp.model.startpar()
        setIndicFix(model) <- FALSE
	prior <- prior(hier = FALSE)
        prior <- priordefine(data, model, varargin = prior)
        mcmc <- .setUp.mcmc.startpar()
	setStorepost(mcmc) <- TRUE
       	(data ~ model ~ mcmc) %=% mcmcstart(data,model,mcmc)
	mcmcout <- mixturemcmc(data,model,prior,mcmc)
        ## Test cases ##
	checkTrue(identical(class(mcmcout)[1], "mcmcoutputpost"), "check1")
        checkTrue(mcmcout@ranperm == FALSE, "check2")
        checkEquals(mcmcout@M, mcmc@M)
	checkTrue(is.list(mcmcout@par), "check3")
	checkTrue("lambda" %in% names(mcmcout@par), "check4")
        checkEquals(ncol(mcmcout@par$lambda), model@K)
        checkEquals(nrow(mcmcout@par$lambda), mcmc@M)
	checkTrue(!any(is.na(mcmcout@par$lambda)), "check5")
	for(k in 1:model@K) {
		checkTrue(all(abs(mcmcout@par$lambda[,k] - 
			mean(mcmcout@par$lambda[,k])) > 0), 
			paste("check", 5 + k, sep =""))
	}
	checkTrue(is.list(mcmcout@log), "check9")
	checkTrue("mixlik" %in% names(mcmcout@log),"check10")
	checkTrue("mixprior" %in% names(mcmcout@log), "check11")
	checkTrue("cdpost" %in% names(mcmcout@log), "check12")
	checkEquals(nrow(mcmcout@log$mixlik), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixlik), 1)
	checkEquals(nrow(mcmcout@log$mixprior), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixprior), 1)
	checkEquals(nrow(mcmcout@log$cdpost), mcmc@M)
	checkEquals(ncol(mcmcout@log$cdpost), 1)
	checkTrue(!any(is.na(mcmcout@log$mixlik)), "check13")
	checkTrue(!any(is.na(mcmcout@log$mixprior)), "check14")
	checkTrue(!any(is.na(mcmcout@log$cdpost)), "check15")
	checkTrue(all(abs(mcmcout@log$mixlik - mean(mcmcout@log$mixlik)) 
			> 0), "check16")
	checkTrue(all(abs(mcmcout@log$mixprior - mean(mcmcout@log$mixprior)) 
			> 0), "check17")
	checkTrue(all(abs(mcmcout@log$cdpost - mean(mcmcout@log$cdpost))
			> 0), "check18")
	checkEquals(nrow(mcmcout@weight), mcmc@M)
	checkEquals(ncol(mcmcout@weight), model@K)
	m <- 2
	for(k in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@weight[, k])), 
			paste("check", 18 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@weight[, k] 
			- mean(mcmcout@weight[, k])) > 0),
			paste("check", 19 + m * (k - 1) + k, sep = ""))
	}
	checkTrue(!.hasSlot(mcmcout, "hyper"), "check25")
	checkTrue(is.list(mcmcout@post), "check26")	
	checkTrue("par" %in% names(mcmcout@post), "check27")
	checkTrue("weight" %in% names(mcmcout@post), "check28")
	checkEquals(nrow(mcmcout@post$weight), mcmc@M)
	checkEquals(ncol(mcmcout@post$weight), model@K)
	checkTrue(is.list(mcmcout@post$par), "check29")
	checkTrue("a" %in% names(mcmcout@post$par), "check30")
	checkTrue("b" %in% names(mcmcout@post$par), "check31")
	checkEquals(nrow(mcmcout@post$par$a), mcmc@M)
	checkEquals(nrow(mcmcout@post$par$b), mcmc@M)
	checkEquals(ncol(mcmcout@post$par$a), model@K)
	checkEquals(ncol(mcmcout@post$par$b), model@K)
	m <- 4
	for(k in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@post$par$a[, k])), 
			paste("check", 31 + m * (k - 1) + k, sep = ""))
		checkTrue(!any(is.na(mcmcout@post$par$b[, k])), 
			paste("check", 32 + m * (k - 1) + k, sep = ""))
		checkTrue(!any(is.na(mcmcout@post$weight[, k])),
			paste("check", 33 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$par$a[, k] - 
			mean(mcmcout@post$par$a[, k])) > 0), 
			paste("check", 34 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$par$b[, k] -
			 mean(mcmcout@post$par$b[, k])) > 0), 
			paste("check", 35 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$weight[, k] 
			- mean(mcmcout@post$weight[, k])) > 0),
			paste("check", 36 + m * (k - 1) + k, sep = ""))
	}
	checkEquals(nrow(mcmcout@entropy), mcmc@M)
	checkEquals(ncol(mcmcout@entropy), 1)
	checkTrue(!any(is.na(mcmcout@entropy)), "check50")
	checkTrue(all(abs(mcmcout@entropy 
		- mean(mcmcout@entropy)) > 0), "check51")
	checkEquals(nrow(mcmcout@ST), mcmc@M)
	checkEquals(ncol(mcmcout@ST), 1)	
	checkTrue(!any(is.na(mcmcout@ST)), "check52")
	checkTrue(all(abs(mcmcout@ST - mean(mcmcout@ST)) 
			> 0), "check53")
	checkEquals(nrow(mcmcout@S), data@N)
	checkEquals(ncol(mcmcout@S), mcmc@storeS)
	checkTrue(!any(is.na(mcmcout@S[, 1])), "check54")
	checkTrue(all(abs(mcmcout@S[, 1] 
		- mean(mcmcout@S[,1])) == 0), "check55")
	m <- 2
	for(s in 2:mcmc@storeS) {
		checkTrue(!any(is.na(mcmcout@S[, s])), 
			paste("check", 55 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@S[, s] 
			- mean(mcmcout@S[, s])) > 0),
			paste("check", 56 + m * (k - 1) + k, sep = ""))
	}
	checkEquals(nrow(mcmcout@NK), mcmc@M)
	checkEquals(ncol(mcmcout@NK), model@K)
	m <- 2
	for(n in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@NK[, n])),
			paste("check", 59 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@NK[, n]
			- mean(mcmcout@NK[, n])) > 0),
			paste("check", 60 + m * (k - 1) + k, sep = ""))
	}
	checkEquals(nrow(mcmcout@clust), data@N)
	checkEquals(ncol(mcmcout@clust), 1)
	checkTrue(!any(is.na(mcmcout@clust)), "check67")
	checkTrue(all(abs(mcmcout@clust - mean(mcmcout@clust)) 
			> 0), "check68")
	checkTrue(identical(mcmcout@model, model), "check69")
	checkTrue(identical(mcmcout@prior, prior), "check70")	
}

"test.IND.HIER.POST.STARTPAR.1"  <- function() {
        ## Setting:
        ##      indicfix: FALSE
 	set.seed(0)       
        data <- .setUp.data.startpar(withInd = FALSE)
        model <- .setUp.model.startpar()
	setK(model) <- 1
        setIndicFix(model) <- FALSE
        prior <- priordefine(data, model)
        mcmc <- .setUp.mcmc.startpar()
	setStorepost(mcmc) <- TRUE
       	(data ~ model ~ mcmc) %=% mcmcstart(data,model,mcmc)
	mcmcout <- mixturemcmc(data,model,prior,mcmc)
        ## Test cases ##
	checkTrue(identical(class(mcmcout)[1], "mcmcoutputfixhierpost"), "check1")
        checkTrue(mcmcout@ranperm == FALSE, "check2")
        checkEquals(mcmcout@M, mcmc@M)
	checkTrue(is.list(mcmcout@par), "check3")
	checkTrue("lambda" %in% names(mcmcout@par), "check4")
        checkEquals(ncol(mcmcout@par$lambda), model@K)
        checkEquals(nrow(mcmcout@par$lambda), mcmc@M)
	checkTrue(!any(is.na(mcmcout@par$lambda)), "check5")
	for(k in 1:model@K) {
		checkTrue(all(abs(mcmcout@par$lambda[,k] - 
			mean(mcmcout@par$lambda[,k])) > 0), 
			paste("check", 5 + k, sep =""))
	}
	checkTrue(is.list(mcmcout@log), "check7")
	checkTrue("mixlik" %in% names(mcmcout@log),"check8")
	checkTrue("mixprior" %in% names(mcmcout@log), "check9")
	checkTrue(!("cdpost" %in% names(mcmcout@log)), "check10")
	checkEquals(nrow(mcmcout@log$mixlik), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixlik), 1)
	checkEquals(nrow(mcmcout@log$mixprior), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixprior), 1)
	checkTrue(!any(is.na(mcmcout@log$mixlik)), "check11")
	checkTrue(!any(is.na(mcmcout@log$mixprior)), "check12")
	checkTrue(all(abs(mcmcout@log$mixlik - mean(mcmcout@log$mixlik)) 
			> 0), "check13")
	checkTrue(all(abs(mcmcout@log$mixprior - mean(mcmcout@log$mixprior)) 
			> 0), "check14")
	checkTrue(!.hasSlot(mcmcout, "weight"), "check15")
	checkTrue(is.list(mcmcout@hyper), "check16")
	checkTrue("b" %in% names(mcmcout@hyper), "check17")
	checkEquals(nrow(mcmcout@hyper$b), mcmc@M)
	checkEquals(ncol(mcmcout@hyper$b), 1)
	checkTrue(!any(is.na(mcmcout@hyper$b)), "check18")
	checkTrue(all(abs(mcmcout@hyper$b - mean(mcmcout@hyper$b)) 
			> 0), "check18")
	checkTrue(is.list(mcmcout@post), "check19")
	checkTrue("par" %in% names(mcmcout@post), "check20")
	checkTrue(!("weight" %in% names(mcmcout@post)), "check21")
	checkTrue(is.list(mcmcout@post$par), "check22")
	checkTrue("a" %in% names(mcmcout@post$par), "check23")
	checkTrue("b" %in% names(mcmcout@post$par), "check24")
	checkEquals(nrow(mcmcout@post$par$a), mcmc@M)
	checkEquals(nrow(mcmcout@post$par$b), mcmc@M)
	checkEquals(ncol(mcmcout@post$par$a), model@K)
	checkEquals(ncol(mcmcout@post$par$b), model@K)
	for(k in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@post$par$a[, k])), 
			paste("check", 24 + k, sep = ""))
		checkTrue(!any(is.na(mcmcout@post$par$b[, k])), 
			paste("check", 25 + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$par$a[, k] - 
			mean(mcmcout@post$par$a[, k])) == 0), 
			paste("check", 26 + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$par$b[, k] -
			 mean(mcmcout@post$par$b[, k])) > 0), 
			paste("check", 27 + k, sep = ""))
	}
	checkTrue(!.hasSlot(mcmcout, "entropy"), "check29")
	checkTrue(!.hasSlot(mcmcout, "ST"), "check30")
	checkTrue(!.hasSlot(mcmcout, "S"), "check31")
	checkTrue(!.hasSlot(mcmcout, "NK"), "check32")
	checkTrue(!.hasSlot(mcmcout, "clust"), "check33")
	checkTrue(identical(mcmcout@model, model), "check34")
	checkTrue(identical(mcmcout@prior, prior), "check35")	
}

"test.IND.HIER.POST.STARTPAR.2" <- function() {
        ## Setting:
        ##      indicfix: FALSE
        set.seed(0)
	data <- .setUp.data.startpar(withInd = TRUE)
        model <- .setUp.model.startpar()
        setIndicFix(model) <- FALSE
        prior <- priordefine(data, model)
        mcmc <- .setUp.mcmc.startpar()
	setStorepost(mcmc) <- TRUE
       	(data ~ model ~ mcmc) %=% mcmcstart(data,model,mcmc)
	mcmcout <- mixturemcmc(data,model,prior,mcmc)
        ## Test cases ##
	checkTrue(identical(class(mcmcout)[1], "mcmcoutputhierpost"), "check1")
        checkTrue(mcmcout@ranperm == FALSE, "check2")
        checkEquals(mcmcout@M, mcmc@M)
	checkTrue(is.list(mcmcout@par), "check3")
	checkTrue("lambda" %in% names(mcmcout@par), "check4")
        checkEquals(ncol(mcmcout@par$lambda), model@K)
        checkEquals(nrow(mcmcout@par$lambda), mcmc@M)
	checkTrue(!any(is.na(mcmcout@par$lambda)), "check5")
	for(k in 1:model@K) {
		checkTrue(all(abs(mcmcout@par$lambda[,k] - 
			mean(mcmcout@par$lambda[,k])) > 0), 
			paste("check", 5 + k, sep =""))
	}
	checkTrue(is.list(mcmcout@log), "check8")
	checkTrue("mixlik" %in% names(mcmcout@log),"check9")
	checkTrue("mixprior" %in% names(mcmcout@log), "check10")
	checkTrue("cdpost" %in% names(mcmcout@log), "check11")
	checkEquals(nrow(mcmcout@log$mixlik), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixlik), 1)
	checkEquals(nrow(mcmcout@log$mixprior), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixprior), 1)
	checkEquals(nrow(mcmcout@log$cdpost), mcmc@M)
	checkEquals(ncol(mcmcout@log$cdpost), 1)
	checkTrue(!any(is.na(mcmcout@log$mixlik)), "check12")
	checkTrue(!any(is.na(mcmcout@log$mixprior)), "check13")
	checkTrue(!any(is.na(mcmcout@log$cdpost)), "check14")
	checkTrue(all(abs(mcmcout@log$mixlik - mean(mcmcout@log$mixlik)) 
			> 0), "check15")
	checkTrue(all(abs(mcmcout@log$mixprior - mean(mcmcout@log$mixprior)) 
			> 0), "check16")
	checkTrue(all(abs(mcmcout@log$cdpost - mean(mcmcout@log$cdpost))
			> 0), "check17")
	checkEquals(nrow(mcmcout@weight), mcmc@M)
	checkEquals(ncol(mcmcout@weight), model@K)
	m <- 2
	for(k in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@weight[, k])), 
			paste("check", 17 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@weight[, k] 
			- mean(mcmcout@weight[, k])) > 0),
			paste("check", 18 + m * (k - 1) + k, sep = ""))
	}
	checkTrue(is.list(mcmcout@hyper), "check22")
	checkTrue("b" %in% names(mcmcout@hyper), "check23")
	checkEquals(nrow(mcmcout@hyper$b), mcmc@M)
	checkEquals(ncol(mcmcout@hyper$b), 1)
	checkTrue(!any(is.na(mcmcout@hyper$b)), "check24")
	checkTrue(all(abs(mcmcout@hyper$b - mean(mcmcout@hyper$b)) 
			> 0), "check25")
	checkTrue(is.list(mcmcout@post), "check26")	
	checkTrue("par" %in% names(mcmcout@post), "check27")
	checkTrue("weight" %in% names(mcmcout@post), "check28")
	checkEquals(nrow(mcmcout@post$weight), mcmc@M)
	checkEquals(ncol(mcmcout@post$weight), model@K)
	checkTrue(is.list(mcmcout@post$par), "check29")
	checkTrue("a" %in% names(mcmcout@post$par), "check30")
	checkTrue("b" %in% names(mcmcout@post$par), "check31")
	checkEquals(nrow(mcmcout@post$par$a), mcmc@M)
	checkEquals(nrow(mcmcout@post$par$b), mcmc@M)
	checkEquals(ncol(mcmcout@post$par$a), model@K)
	checkEquals(ncol(mcmcout@post$par$b), model@K)
	m <- 4
	for(k in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@post$par$a[, k])), 
			paste("check", 31 + m * (k - 1) + k, sep = ""))
		checkTrue(!any(is.na(mcmcout@post$par$b[, k])), 
			paste("check", 32 + m * (k - 1) + k, sep = ""))
		checkTrue(!any(is.na(mcmcout@post$weight[, k])),
			paste("check", 33 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$par$a[, k] - 
			mean(mcmcout@post$par$a[, k])) > 0), 
			paste("check", 34 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$par$b[, k] -
			 mean(mcmcout@post$par$b[, k])) > 0), 
			paste("check", 35 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$weight[, k] 
			- mean(mcmcout@post$weight[, k])) > 0),
			paste("check", 36 + m * (k - 1) + k, sep = ""))
	}
	checkEquals(nrow(mcmcout@entropy), mcmc@M)
	checkEquals(ncol(mcmcout@entropy), 1)
	checkTrue(!any(is.na(mcmcout@entropy)), "check44")
	checkTrue(all(abs(mcmcout@entropy 
		- mean(mcmcout@entropy)) > 0), "check45")
	checkEquals(nrow(mcmcout@ST), mcmc@M)
	checkEquals(ncol(mcmcout@ST), 1)	
	checkTrue(!any(is.na(mcmcout@ST)), "check46")
	checkTrue(all(abs(mcmcout@ST - mean(mcmcout@ST)) 
			> 0), "check47")
	checkEquals(nrow(mcmcout@S), data@N)
	checkEquals(ncol(mcmcout@S), mcmc@storeS)
	checkTrue(!any(is.na(mcmcout@S[, 1])), "check28")
	checkTrue(all(abs(mcmcout@S[, 1] 
		- mean(mcmcout@S[, 1])) == 0), "check29")
	m <- 2
	for(s in 2:mcmc@storeS) {
		checkTrue(!any(is.na(mcmcout@S[, s])), 
			paste("check", 47 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@S[, s] 
			- mean(mcmcout@S[, s])) > 0),
			paste("check", 48 + m * (k - 1) + k, sep = ""))
	}
	checkEquals(nrow(mcmcout@NK), mcmc@M)
	checkEquals(ncol(mcmcout@NK), model@K)
	m <- 2
	for(n in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@NK[, n])),
			paste("check", 51 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@NK[, n]
			- mean(mcmcout@NK[, n])) > 0),
			paste("check", 52 + m * (k - 1) + k, sep = ""))
	}
	checkEquals(nrow(mcmcout@clust), data@N)
	checkEquals(ncol(mcmcout@clust), 1)
	checkTrue(!any(is.na(mcmcout@clust)), "check56")
	checkTrue(all(abs(mcmcout@clust - mean(mcmcout@clust)) 
			> 0), "check57")
	checkTrue(identical(mcmcout@model, model), "check58")
	checkTrue(identical(mcmcout@prior, prior), "check59")	
}

"test.IND.HIER.POST.STARTPAR.3" <- function() {
        ## Setting:
        ##      indicfix: FALSE
       	set.seed(0)       
	data <- .setUp.data.startpar(withInd = TRUE)
        model <- .setUp.model.startpar()
        setIndicFix(model) <- FALSE
        prior <- priordefine(data, model)
        mcmc <- .setUp.mcmc.startpar()
	setStorepost(mcmc) <- TRUE
 	(data ~ model ~ mcmc) %=% mcmcstart(data,model,mcmc)
      	mcmcout <- mixturemcmc(data,model,prior,mcmc)
        ## Test cases ##
	checkTrue(identical(class(mcmcout)[1], "mcmcoutputhierpost"), "check1")
        checkTrue(mcmcout@ranperm == FALSE, "check2")
        checkEquals(mcmcout@M, mcmc@M)
	checkTrue(is.list(mcmcout@par), "check3")
	checkTrue("lambda" %in% names(mcmcout@par), "check4")
        checkEquals(ncol(mcmcout@par$lambda), model@K)
        checkEquals(nrow(mcmcout@par$lambda), mcmc@M)
	checkTrue(!any(is.na(mcmcout@par$lambda)), "check5")
	for(k in 1:model@K) {
		checkTrue(all(abs(mcmcout@par$lambda[,k] - 
			mean(mcmcout@par$lambda[,k])) > 0), 
			paste("check", 5 + k, sep =""))
	}
	checkTrue(is.list(mcmcout@log), "check9")
	checkTrue("mixlik" %in% names(mcmcout@log),"check10")
	checkTrue("mixprior" %in% names(mcmcout@log), "check11")
	checkTrue("cdpost" %in% names(mcmcout@log), "check12")
	checkEquals(nrow(mcmcout@log$mixlik), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixlik), 1)
	checkEquals(nrow(mcmcout@log$mixprior), mcmc@M)
	checkEquals(ncol(mcmcout@log$mixprior), 1)
	checkEquals(nrow(mcmcout@log$cdpost), mcmc@M)
	checkEquals(ncol(mcmcout@log$cdpost), 1)
	checkTrue(!any(is.na(mcmcout@log$mixlik)), "check13")
	checkTrue(!any(is.na(mcmcout@log$mixprior)), "check14")
	checkTrue(!any(is.na(mcmcout@log$cdpost)), "check15")
	checkTrue(all(abs(mcmcout@log$mixlik - mean(mcmcout@log$mixlik)) 
			> 0), "check16")
	checkTrue(all(abs(mcmcout@log$mixprior - mean(mcmcout@log$mixprior)) 
			> 0), "check17")
	checkTrue(all(abs(mcmcout@log$cdpost - mean(mcmcout@log$cdpost))
			> 0), "check18")
	checkEquals(nrow(mcmcout@weight), mcmc@M)
	checkEquals(ncol(mcmcout@weight), model@K)
	m <- 2
	for(k in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@weight[, k])), 
			paste("check", 18 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@weight[, k] 
			- mean(mcmcout@weight[, k])) > 0),
			paste("check", 19 + m * (k - 1) + k, sep = ""))
	}
	checkTrue(is.list(mcmcout@hyper), "check23")
	checkTrue("b" %in% names(mcmcout@hyper), "check24")
	checkEquals(nrow(mcmcout@hyper$b), mcmc@M)
	checkEquals(ncol(mcmcout@hyper$b), 1)
	checkTrue(!any(is.na(mcmcout@hyper$b)), "check25")
	checkTrue(all(abs(mcmcout@hyper$b - mean(mcmcout@hyper$b)) 
			> 0), "check26")
	checkTrue(is.list(mcmcout@post), "check27")	
	checkTrue("par" %in% names(mcmcout@post), "check28")
	checkTrue("weight" %in% names(mcmcout@post), "check29")
	checkEquals(nrow(mcmcout@post$weight), mcmc@M)
	checkEquals(ncol(mcmcout@post$weight), model@K)
	checkTrue(is.list(mcmcout@post$par), "check30")
	checkTrue("a" %in% names(mcmcout@post$par), "check31")
	checkTrue("b" %in% names(mcmcout@post$par), "check32")
	checkEquals(nrow(mcmcout@post$par$a), mcmc@M)
	checkEquals(nrow(mcmcout@post$par$b), mcmc@M)
	checkEquals(ncol(mcmcout@post$par$a), model@K)
	checkEquals(ncol(mcmcout@post$par$b), model@K)
	m <- 4
	for(k in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@post$par$a[, k])), 
			paste("check", 32 + m * (k - 1) + k, sep = ""))
		checkTrue(!any(is.na(mcmcout@post$par$b[, k])), 
			paste("check", 33 + m * (k - 1) + k, sep = ""))
		checkTrue(!any(is.na(mcmcout@post$weight[, k])),
			paste("check", 34 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$par$a[, k] - 
			mean(mcmcout@post$par$a[, k])) > 0), 
			paste("check", 35 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$par$b[, k] -
			 mean(mcmcout@post$par$b[, k])) > 0), 
			paste("check", 36 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@post$weight[, k] 
			- mean(mcmcout@post$weight[, k])) > 0),
			paste("check", 37 + m * (k - 1) + k, sep = ""))
	}
	checkEquals(nrow(mcmcout@entropy), mcmc@M)
	checkEquals(ncol(mcmcout@entropy), 1)
	checkTrue(!any(is.na(mcmcout@entropy)), "check51")
	checkTrue(all(abs(mcmcout@entropy 
		- mean(mcmcout@entropy)) > 0), "check52")
	checkEquals(nrow(mcmcout@ST), mcmc@M)
	checkEquals(ncol(mcmcout@ST), 1)	
	checkTrue(!any(is.na(mcmcout@ST)), "check53")
	checkTrue(all(abs(mcmcout@ST - mean(mcmcout@ST)) 
			> 0), "check54")
	checkEquals(nrow(mcmcout@S), data@N)
	checkEquals(ncol(mcmcout@S), mcmc@storeS)
	checkTrue(!any(is.na(mcmcout@S[, 1])), "check28")
	checkTrue(all(abs(mcmcout@S[, 1] 
		- mean(mcmcout@S[, 1])) == 0), "check29")
	m <- 2
	for(s in 2:mcmc@storeS) {
		checkTrue(!any(is.na(mcmcout@S[, s])), 
			paste("check", 54 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@S[, s] 
			- mean(mcmcout@S[, s])) > 0),
			paste("check", 55 + m * (k - 1) + k, sep = ""))
	}
	checkEquals(nrow(mcmcout@NK), mcmc@M)
	checkEquals(ncol(mcmcout@NK), model@K)
	m <- 2
	for(n in 1:model@K) {
		checkTrue(!any(is.na(mcmcout@NK[, n])),
			paste("check", 60 + m * (k - 1) + k, sep = ""))
		checkTrue(all(abs(mcmcout@NK[, n]
			- mean(mcmcout@NK[, n])) > 0),
			paste("check", 61 + m * (k - 1) + k, sep = ""))
	}
	checkEquals(nrow(mcmcout@clust), data@N)
	checkEquals(ncol(mcmcout@clust), 1)
	checkTrue(!any(is.na(mcmcout@clust)), "check67")
	checkTrue(all(abs(mcmcout@clust - mean(mcmcout@clust)) 
			> 0), "check68")
	checkTrue(identical(mcmcout@model, model), "check69")
	checkTrue(identical(mcmcout@prior, prior), "check70")	
}
