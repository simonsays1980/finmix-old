#ifndef FINMIXPOST_H
#define FINMIXPOST_H

#include <RcppArmadillo.h> 

template <typename PostParType>
class FinmixPost {
	public:
		PostParType postPar;
		arma::mat weight;
		FinmixPost(const Rcpp::List &post);
		~FinmixPost();
};

template <typename PostParType> 
FinmixPost<PostParType>::FinmixPost (const Rcpp::List& post) 
{
	Rcpp::List tmpPP(post["par"]);
	Rcpp::NumericMatrix tmpW(post["weight"]);
	unsigned int M = tmpW.rows();
	unsigned int K = tmpW.cols();
	postPar(tmpPP);
	weight(tmpW.begin(), M, K, false, true);
};

#endif
