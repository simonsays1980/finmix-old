/**
 * methods to sample classifications of finite mixtures 
 *  author: Lars Simon Zehnder  
 * package: finmix
 * created: 02/01/2013
 */

#include <RcppArmadillo.h>
//#include "liklist.h"    // structure liklist
#include "likelihood.h" // methods for calculating likelihoods

/* structure to contain:
 * logpy    : log likelihood values for 1 to K
 * prob     : classification probability matrix
 * newS     : sampled classifications from classification probability matrix
 * loglikcd : posterior complete data likelihood, knowing classifications
 * mixlik   : posterior mixture likelihood
 * entropy  : entropy of the posterior classification probability distribution
 * postS    : posterior of sampled classifications
 */
struct DataClass {

	arma::mat logPy; //TODO: Check if this is needed! 
	arma::mat prob;  //TODO: Check if this is needed!
	arma::uvec newS;
	arma::vec logLikCd;
	double mixLik;
	double entropy;
	double postS;
	
};

DataClass classification(const bool&, const arma::uvec&, const liklist&, const arma::rowvec&);

inline DataClass 
dataclass_poisson (const bool& INDICFIX, const arma::mat& yM, 
			const arma::uvec& sV, const arma::rowvec& weightV, 
			const arma::rowvec& parV) 
{

	liklist lik_l = likelihood_poisson(yM, parV);
	DataClass dataC = classification(INDICFIX, sV, lik_l, weightV);        
	return dataC;	
}

inline DataClass 
classification (const bool &INDICFIX, const arma::uvec &S, 
			const liklist &lik, const arma::rowvec &weight) 
{

	const unsigned int N = S.n_elem;
	const unsigned int K = weight.n_elem;
	
	/* if indicators are not fix, they are simulated */
	if(!INDICFIX) {
		arma::mat p_m(N, K); 
		arma::uvec newS(N);
		//TODO: Check why it is unused
		double postS = 0.0;
		/* only multinomial incidator model implemented */
		for(unsigned int i = 0; i < N; ++i) {
			p_m.row(i) = lik.lh.row(i) % weight;
		}
		arma::vec sump_v = sum(lik.lh, 1);        	// N x 1 matrix
		arma::vec lsump = arma::log(sump_v) + lik.maxl;   // N x 1 matrix
		double mixlik = sum(lsump);               	// mixture likelihood
		p_m.each_col() /= sump_v;	                  	// classification probability matrix
						          	// N x K matrix

		
		/* simulate only if true mixture */
		if(K > 1) {
			/* simulate classifications from probability matrix p */
			arma::vec rnd = arma::randu<arma::vec>(N); 	// draw N uniform random numbers 
							     	// TODO: check rng from armadillo!!
			arma::mat rndM = arma::repmat(rnd, 1, K);  	// N x K matrix
			arma::mat cumSP = arma::cumsum(p_m, 1);      	// cumulate along rows, N x K matrix
			arma::umat ind = (cumSP < rndM);
			rndM = arma::conv_to<arma::mat>::from(ind);    	// logical N x K matrix
			newS = arma::conv_to<arma::uvec>::from(arma::sum(rndM, 1));        	// new classifications
			
			/* compute posterior log likelihood of S */
			// TODO: check if this is necessary! 
			arma::umat Sm = arma::repmat(newS, 1, K);      // N x K matrix of S
			arma::umat compM = arma::ones<arma::umat>(N, K);
			for(unsigned int k = 0; k < K; ++k) {
				compM.col(k) * (k + 1);
			}
			ind = (Sm == compM);      	     // logical N x K matrix
			arma::mat indDouble = arma::conv_to<arma::mat>::from(ind);
			arma::vec postSm = arma::sum(indDouble % p_m, 1);  // sum along rows
			postSm = arma::log(postSm);
			double postS = arma::sum(postSm);
		}	
		
		
		/* calculate entropy */
		arma::mat logp(N, K);
		arma::uvec col_index(1);
		for(unsigned int k = 0; k < K; ++k) {
			col_index(0) = k;
			arma::uvec zero_index = arma::find(p_m.col(k) == 0);
			arma::uvec index = arma::find(p_m.col(k));
			logp.submat(zero_index, col_index).fill(-99.00);
			logp.submat(index, col_index) = arma::log(p_m.submat(index, col_index));
			
		} 
		double entropy = (-1.0) * arma::accu(logp % p_m);
	
		DataClass dataC = DataClass();
		dataC.logPy = lik.llh;
		dataC.prob = p_m;
		dataC.newS = newS;
		dataC.mixLik = mixlik;
		dataC.entropy = entropy;
		dataC.postS = postS;

		return dataC;
	}
	else { /* indicators fixed */
		arma::vec loglikcd(1, K);
		arma::uvec col_index(1);
		if(!S.is_empty() && K > 1) {
		 	for(unsigned int k = 0; k < K; ++k) {
				col_index(0) = k;
				arma::uvec index = arma::find(S == k);
                                loglikcd(k) = arma::accu(lik.llh.submat(index, col_index));
			}
		}
		else { /* no indicators or no true mixture */
			arma::vec sump_v = sum(lik.lh, 1);        	// N x 1 matrix
			arma::vec lsump = arma::log(sump_v) + lik.maxl;   // N x 1 matrix
			double mixlik = sum(lsump);               	// mixture likelihood
			//TODO: check if better to fill all K entries with mixlik
			loglikcd(0, 0) = mixlik;
		}
		DataClass dataC;
		dataC.logPy = lik.llh;
		dataC.logLikCd = loglikcd;

		return dataC;
	} 
}
