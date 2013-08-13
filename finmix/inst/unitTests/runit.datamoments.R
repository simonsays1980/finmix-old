### --- Test Setup --- ###

if(TRUE) {
	## Not really needed, but can be handy 
	## when writing tests 
	library("RUnit")
	library("finmix")
}

".setUp.data" <- function(withInd = FALSE) {
        ## Get path ##
        pkg <- "finmix"
        if (Sys.getenv("RCMDCHECK") == FALSE) {
            data.path <- file.path(getwd(), "..", 
                                   "data", "poisson.data.csv")        
        } else {
            data.path <- system.file(package = pkg, 
                                     'data/poisson.data.csv')
        }
        data <- read.csv(data.path, header = FALSE, sep = ",")
        if(withInd) {
            if (Sys.getenv("RCMDCHECK") == FALSE) {

                ind.path <- file.path(getwd(), "..", 
                                       "data", 
                                       "poisson.ind.csv")        
            } else {
                ind.path <- system.file(package = pkg, 
                                        'data/poisson.ind.csv')
            }               
            ind <- read.csv(ind.path, header = FALSE, sep = ",")
            data <- data(y = as.matrix(data), S = as.matrix(ind), type = "discrete",
                         r = 1, N = nrow(data), sim = TRUE,
                         bycolumn = TRUE)
            return(data)
        } else {
                data <- data(y = as.matrix(data), type = "discrete", r = 1,
                             N = nrow(data), sim = TRUE,
                             bycolumn. = TRUE)
                return(data)
        }
}

"test.groupmoments" <- function() {
    ## Setup test ##
    data <- .setUp.data(withInd = TRUE)
    mom <- groupmoments(data)
    ## --- Check r = 1 --- ##
    checkTrue(class(mom) == "groupmoments", "check1")
    checkTrue("NK" %in% slotNames(mom), "check2")
    checkTrue("mean" %in% slotNames(mom), "check3")
    checkTrue("WK" %in% slotNames(mom), "check4")
    checkTrue("var" %in% slotNames(mom), "check5")
    checkTrue("data" %in% slotNames(mom), "check6")
    checkTrue(!all(is.na(mom@NK)), "check7")
    checkTrue(!all(is.na(mom@mean)), "check8")
    checkTrue(!all(is.na(mom@WK)), "check9")
    checkTrue(!all(is.na(mom@var)), "check10")
    checkEquals(dim(mom@NK)[1], 2)
    checkEquals(nrow(mom@mean), 1)
    checkEquals(ncol(mom@mean), 2)
    checkEquals(dim(mom@WK)[1], 1)
    checkEquals(dim(mom@WK)[2], 1)
    checkEquals(dim(mom@WK)[3], 2)
    checkEquals(dim(mom@var)[1], 1)
    checkEquals(dim(mom@var)[2], 1)
    checkEquals(dim(mom@var)[3], 2)
    ## --- Check r == 2 --- ##
    setY(data) <- cbind(data@y, rev(data@y))
    setR(data) <- 2
    mom <- groupmoments(data)
    checkTrue(class(mom) == "groupmoments", "check11")
    checkTrue("NK" %in% slotNames(mom), "check12")
    checkTrue("mean" %in% slotNames(mom), "check13")
    checkTrue("WK" %in% slotNames(mom), "check14")
    checkTrue("var" %in% slotNames(mom), "check15")
    checkTrue("data" %in% slotNames(mom), "check16")
    checkTrue(!all(is.na(mom@NK)), "check17")
    checkTrue(!all(is.na(mom@mean)), "check18")
    checkTrue(!all(is.na(mom@WK)), "check19")
    checkTrue(!all(is.na(mom@var)), "check20")
    checkEquals(dim(mom@NK)[1], 2)
    checkEquals(nrow(mom@mean), 2)
    checkEquals(ncol(mom@mean), 2)
    checkEquals(dim(mom@WK)[1], 2)
    checkEquals(dim(mom@WK)[2], 2)
    checkEquals(dim(mom@WK)[3], 2)
    checkEquals(dim(mom@var)[1], 2)
    checkEquals(dim(mom@var)[2], 2)
    checkEquals(dim(mom@var)[3], 2)
    ## Check Exceptions ##
    setS(data) <- matrix()
    checkException(groupmoments(data), "check21")
    data <- data()
    checkException(groupmoments(data), "check22")
}


"test.sdatamoments" <- function() {
    ## Setup test ## 
    data <- .setUp.data(withInd = TRUE)
    mom <- sdatamoments(data)
    ## --- Check r = 1 --- ##
    checkTrue(class(mom) == "sdatamoments", "check1")
    checkTrue("data" %in% slotNames(mom), "check2")
    checkTrue("gmoments" %in% slotNames(mom), "check3")
    ## --- Check r = 1 && type = "continuous" --- ##
    setType(data) <- "continuous"
    mom <- sdatamoments(data)
    checkTrue(class(mom) == "csdatamoments", "check4")
    checkTrue("B" %in% slotNames(mom), "check5")
    checkTrue("W" %in% slotNames(mom), "check6")
    checkTrue("T" %in% slotNames(mom), "check7")
    checkTrue("R" %in% slotNames(mom), "check8")
    checkTrue("Rdet" %in% slotNames(mom), "check9")
    checkTrue("Rtr" %in% slotNames(mom), "check10")
    checkTrue("data" %in% slotNames(mom), "check11")
    checkTrue("gmoments" %in% slotNames(mom), "check12")
    checkTrue(class(mom@gmoments) == "groupmoments", "check13")
    checkTrue(!all(is.na(mom@B)), "check14")
    checkTrue(!all(is.na(mom@W)), "check15")
    checkTrue(!all(is.na(mom@T)), "check16")
    checkTrue(!is.na(mom@R), "check17")
    checkTrue(is.na(mom@Rdet), "check18")
    checkTrue(is.na(mom@Rtr), "check19")
    checkEquals(dim(mom@W)[2], 1)
    checkEquals(dim(mom@W)[1], 1)
    checkEquals(dim(mom@B)[2], 1)
    checkEquals(dim(mom@B)[1], 1)
    checkEquals(dim(mom@T)[2], 1)
    checkEquals(dim(mom@T)[1], 1)
    ## --- Check for r = 2 --- ##
    setType(data) <- "discrete"
    setY(data) <- cbind(data@y, rev(data@y))
    setR(data) <- 2
    mom <- sdatamoments(data)
    checkTrue(class(mom) == "sdatamoments", "check20")
    checkTrue("data" %in% slotNames(mom), "check21")
    checkTrue("gmoments" %in% slotNames(mom), "check22")
    ## --- Check r = 2 && type = "continuous" --- ##
    setType(data) <- "continuous"
    mom <- sdatamoments(data)
    checkTrue(class(mom) == "csdatamoments", "check23")
    checkTrue("B" %in% slotNames(mom), "check24")
    checkTrue("W" %in% slotNames(mom), "check25")
    checkTrue("T" %in% slotNames(mom), "check26")
    checkTrue("R" %in% slotNames(mom), "check27")
    checkTrue("Rdet" %in% slotNames(mom), "check28")
    checkTrue("Rtr" %in% slotNames(mom), "check29")
    checkTrue("data" %in% slotNames(mom), "check30")
    checkTrue("gmoments" %in% slotNames(mom), "check31")
    checkTrue(class(mom@gmoments) == "groupmoments", "check32")
    checkTrue(!all(is.na(mom@B)), "check33")
    checkTrue(!all(is.na(mom@W)), "check34")
    checkTrue(!all(is.na(mom@T)), "check35")
    checkTrue(is.na(mom@R), "check36")
    checkTrue(!is.na(mom@Rdet), "check37")
    checkTrue(!is.na(mom@Rtr), "check38")
    checkEquals(dim(mom@W)[2], 2)
    checkEquals(dim(mom@W)[1], 2)
    checkEquals(dim(mom@B)[2], 2)
    checkEquals(dim(mom@B)[1], 2)
    checkEquals(dim(mom@T)[2], 2)
    checkEquals(dim(mom@T)[1], 2)
    ## --- Check r = 2 && type = "discrete" --- ##
    setType(data) <- "discrete"
    mom <- sdatamoments(data)
    checkTrue("data" %in% slotNames(mom), "check39")
    checkTrue("gmoments" %in% slotNames(mom), "check40")
    ## Check exceptions ##
    setS(data) <- matrix()
    checkException(sdatamoments(data), "check41")
    data <- data()
    checkException(sdatamoments(data), "check43")
}

"test.datamoments" <- function() {
    ## Setup test ##
    data <- .setUp.data(withInd = TRUE)
    mom <- datamoments(data)
    ## --- Check for r = 1 --- ##
    checkTrue(class(mom) == "ddatamoments", "check1")
    checkTrue("mean" %in% slotNames(mom), "check2")
    checkTrue("var" %in% slotNames(mom), "check3")
    checkTrue("factorial" %in% slotNames(mom), "check4")
    checkTrue("over" %in% slotNames(mom), "check5")
    checkTrue("zero" %in% slotNames(mom), "check6")
    checkEquals(length(mom@mean), 1)
    checkEquals(length(mom@var), 1)
    checkEquals(dim(mom@factorial)[1], 4)
    checkEquals(dim(mom@factorial)[2], 1)
    checkEquals(length(mom@over), 1)
    checkEquals(length(mom@zero), 1)
    checkTrue(class(mom@data) == "data", "check7")
    checkTrue(class(mom@smoments) == "sdatamoments", "check8")
    checkTrue(!all(is.na(mom@mean)), "check9")
    checkTrue(!all(is.na(mom@var)), "check10")
    checkTrue(!all(is.na(mom@factorial)), "check11")
    checkTrue(!all(is.na(mom@over)), "check12")
    checkTrue(!all(is.na(mom@zero)), "check13")
    ## --- Check for r = 2 --- ##
    setY(data) <- cbind(data@y, rev(data@y))
    setR(data) <- 2
    mom <- datamoments(data)
    checkTrue(class(mom) == "ddatamoments", "check14")
    checkTrue("mean" %in% slotNames(mom), "check15")
    checkTrue("var" %in% slotNames(mom), "check16")
    checkTrue("factorial" %in% slotNames(mom), "check17")
    checkTrue("over" %in% slotNames(mom), "check18")
    checkTrue("zero" %in% slotNames(mom), "check19")
    checkEquals(length(mom@mean), 2)
    checkEquals(length(mom@var), 4)
    checkEquals(dim(mom@factorial)[1], 4)
    checkEquals(dim(mom@factorial)[2], 2)
    checkEquals(length(mom@over), 2)
    checkEquals(length(mom@zero), 2)
    checkTrue(class(mom@data) == "data", "check20")
    checkTrue(class(mom@smoments) == "sdatamoments", "check21")
    checkTrue(!all(is.na(mom@mean)), "check22")
    checkTrue(!all(is.na(mom@var)), "check23")
    checkTrue(!all(is.na(mom@factorial)), "check23")
    checkTrue(!all(is.na(mom@over)), "check24")
    checkTrue(!all(is.na(mom@zero)), "check25")
    ## --- Check r = 1 && type = continuous --- ##
    data <- .setUp.data(withInd = TRUE)
    setType(data) <- "continuous" 
    mom <- datamoments(data)
    checkTrue(class(mom) == "cdatamoments", "check26")
    checkTrue("mean" %in% slotNames(mom), "check27")
    checkTrue("var" %in% slotNames(mom), "check28")
    checkTrue("higher" %in% slotNames(mom), "check29")
    checkTrue("skewness" %in% slotNames(mom), "check30")
    checkTrue("kurtosis" %in% slotNames(mom), "check31")
    checkTrue("corr" %in% slotNames(mom), "check32")
    checkEquals(length(mom@mean), 1)
    checkEquals(length(mom@var), 1)
    checkEquals(dim(mom@higher)[1], 4)
    checkEquals(dim(mom@higher)[2], 1)
    checkEquals(length(mom@skewness), 1)
    checkEquals(length(mom@kurtosis), 1)
    checkEquals(length(mom@corr), 1)
    checkTrue(class(mom@data) == "data", "check33")
    checkTrue(class(mom@smoments) == "csdatamoments", "check34")
    checkTrue(!all(is.na(mom@mean)), "check35")
    checkTrue(!all(is.na(mom@var)), "check36")
    checkTrue(!all(is.na(mom@higher)), "check37")
    checkTrue(!all(is.na(mom@skewness)), "check38")
    checkTrue(!all(is.na(mom@kurtosis)), "check39")
    checkTrue(all(is.na(mom@corr)), "check40")
    ## --- Check r = 2 && type = continuous --- ##
    setY(data) <- cbind(data@y, rev(data@y))
    setR(data) <- 2
    mom <- datamoments(data)
    checkTrue(class(mom) == "cdatamoments", "check41")
    checkTrue("mean" %in% slotNames(mom), "check42")
    checkTrue("var" %in% slotNames(mom), "check43")
    checkTrue("higher" %in% slotNames(mom), "check44")
    checkTrue("skewness" %in% slotNames(mom), "check45")
    checkTrue("kurtosis" %in% slotNames(mom), "check46")
    checkEquals(length(mom@mean), 2)
    checkEquals(length(mom@var), 4)
    checkEquals(dim(mom@higher)[1], 4)
    checkEquals(dim(mom@higher)[2], 2)
    checkEquals(length(mom@skewness), 2)
    checkEquals(length(mom@skewness), 2)
    checkTrue(class(mom@data) == "data", "check47")
    checkTrue(class(mom@smoments) == "csdatamoments", "check48")
    checkTrue(!all(is.na(mom@mean)), "check49")
    checkTrue(!all(is.na(mom@var)), "check50")
    checkTrue(!all(is.na(mom@higher)), "check51")
    checkTrue(!all(is.na(mom@skewness)), "check52")
    checkTrue(!all(is.na(mom@kurtosis)), "check53")
    ## --- Check Exceptions --- ##
    data <- data()
    checkException(datamoments(data), "check54")
    
}
