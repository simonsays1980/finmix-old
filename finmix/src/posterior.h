/**
 * computes the moments of a posterior distribution under a 
 * conditionally conjugate prior
 *  author: Lars Simon Zehnder
 * package: 'finmix'
 * created: 31 Jan. 2013
 */

#include <RcppArmadillo.h>

/**
 * posterior for multinomial distribution
 *
 * as this is used for the weights of any mixture 
 * no specific struct is input argument 
 */

inline arma::rowvec 
posterior_multinomial(const unsigned int &K, const arma::uvec &S, 
			const arma::rowvec &weight) 
{
	arma::umat repS = arma::repmat(S, 1, K);
	arma::umat compM = arma::ones<arma::umat>(S.n_rows, K);
	arma::rowvec par_post(K);
 
	/* create sequence */
	for(unsigned int k = 0; k < K; ++k) {
		compM.col(k) = compM.col(k) * (k + 1);		
	}
	arma::umat ind = (repS == compM);
	arma::mat indDouble = arma::conv_to<arma::mat>::from(ind);
	par_post = sum(indDouble); 
	par_post = par_post + weight;
	
	return par_post;
}

/**
 * posterior for poisson distribution 
 *
 * This posterior is specific for a poisson 
 * mixture, therefore a parprior_poisson object
 * is expected as input argument and only 
 * modified within. 
 */
inline arma::mat 
posterior_poisson (const unsigned int &K, const arma::vec &y, 
			const arma::uvec &S, const arma::mat &hyper_start_par) 
{
	/* output matrix has in first column 'a' and in second 'b' */
	if(K == 1) {
		arma::mat par_out(1, 2);
		par_out(0, 0) = hyper_start_par(0, 0) + arma::accu(y);
		par_out(0, 1) = hyper_start_par(0, 1) + y.n_elem;

		return par_out;
	}
	else { // K > 1
		arma::mat repY = arma::repmat(y, 1, K);
		arma::umat repS = arma::repmat(S, 1, K);
		// check if cumsum() is faster
		arma::umat compM = arma::ones<arma::umat>(S.n_elem, K);
		arma::mat par_out(K, 2);

	  	/* create a sequence */ 
		for(unsigned int k = 0; k < K; ++k) {
			compM.col(k) = compM.col(k) * (k + 1);
		}
		arma::umat ind = (repS == compM);
		arma::mat indDouble = arma::conv_to<arma::mat>::from(ind);
		repY %= indDouble;
		arma::rowvec prod = sum(repY, 0);
		arma::rowvec sprod = sum(indDouble, 0);
		par_out.col(0) = hyper_start_par.col(0) + sprod.t();
		par_out.col(1) = hyper_start_par.col(1) + prod.t(); 
		
		return par_out;
	}
}
