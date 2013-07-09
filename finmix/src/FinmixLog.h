#ifndef FINMIXLOG_H
#define FINMIXLOG_H

#include <RcppArmadillo.h>

class FinmixLog {
	public:
		arma::vec mixlik;
		arma::vec mixprior;
		arma::vec cdpost;
		
		FinmixLog(const Rcpp::List& log);
		~FinmixLog();
};

#endif
