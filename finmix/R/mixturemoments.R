"mixturemoments.normal" <- function(model, J, meanm) {
	
	zm <- matrix(0, ncol = J, nrow = 1)
	zm[seq(2,J, by = 2)] <- exp(cumsum(log(seq(1,(J-1), by = 2))))
	moments <- matrix(0, nrow = 1, ncol = J)

	for(m in 2:J) { ## first higher moment is always zero
		diff <- apply(model@par$mu, 2, "-", meanm)
		moments[m] <- sum(model@weight * diff^m) 
		for(n in 1:m) {
			cm = diff^(m - n) * model@par$sigma^(n/2)*zm[n]
			moments[m] = moments[m] + choose(m, n) * sum(model@weight * cm)
		}
	}
	
	return(moments) 
}
