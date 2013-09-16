### --- Test Setup --- ###

if(TRUE) {
	## Not really needed, but can be handy 
	## when writing tests 
	library("RUnit")
	library("finmix")
}

".setUp.data" <- function(dist = "poisson") {
        ## Get path ##
        pkg <- "finmix"
        data.name <- paste(dist, "data.csv", sep = "")
        if (Sys.getenv("RCMDCHECK") == FALSE) {
            data.path <- file.path(getwd(), "..", 
                                   "data", data.name)        
        } else {
            data.path <- system.file(package = pkg, 
                                     paste('data/', data.name))
        }
        data <- read.csv(data.path, header = FALSE, sep = ",")
        if (dist %in% c("poisson", "binomial")) {
            type <- "discrete"
        } else {
            type <- "continuous"
        }
        data <- data(y = as.matrix(data), type = type,
                     r = 1, N = nrow(data), sim = TRUE,
                     bycolumn = TRUE)
        return(data)
}

## Tests ##
"test.priordefine.poisson" <- function() 
{
    ## Setup Test ##
    data    <- .setUp.data("poisson") 
    model   <- model(dist = "poisson", K = 2)
    prior   <- priordefine(data, model)
    ##--- Check K = 2 ---##
    ## Default prior ##
    ## Type is always 'condconjugate' ## 
    checkTrue(class(prior) == "prior", "check1")
    checkTrue(prior@hier, "check2")
    checkTrue(prior@type, "condconjugate", "check3")
    checkTrue(is.list(prior@par), "check4")
    checkTrue(is.matrix(prior@weight), "check5")
    checkTrue("a" %in% names(prior@par), "check6")
    checkTrue("b" %in% names(prior@par), "check7")
    checkTrue("g" %in% names(prior@par), "check8")
    checkTrue("G" %in% names(prior@par), "check9")
    checkTrue(is.array(prior@par$a), "check10")
    checkTrue(is.array(prior@par$b), "check11")
    checkTrue(is.numeric(prior@par$g), "check12")
    checkTrue(is.numeric(prior@par$G), "check13")
    checkEquals(dim(prior@par$a)[2], 2)
    checkEquals(dim(prior@par$b)[2], 2)
    ## Non-hierarchical prior ##
    prior <- priordefine(data, model, 
                         varargin = prior(hier = FALSE))
    checkTrue(class(prior) == "prior", "check1")
    checkTrue(!prior@hier, "check2")
    checkTrue(prior@type, "condconjugate", "check3")
    checkTrue(is.list(prior@par), "check4")
    checkTrue(is.matrix(prior@weight), "check5")
    checkTrue("a" %in% names(prior@par), "check6")
    checkTrue("b" %in% names(prior@par), "check7")
    checkTrue(!("g" %in% names(prior@par)), "check8")
    checkTrue(!("G" %in% names(prior@par)), "check9")
    checkTrue(is.array(prior@par$a), "check10")
    checkTrue(is.array(prior@par$b), "check11")
    checkEquals(dim(prior@par$a)[2], 2)
    checkEquals(dim(prior@par$b)[2], 2)
   
    ##--- Check K = 1 ---##
    model@K <- as.integer(1)
    prior   <- prior(data, model)
}
