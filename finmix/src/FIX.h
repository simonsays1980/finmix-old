#ifndef FIX_H
#define FIX_H

#include <RcppArmadillo.h>
#include "FinmixData.h"
#include "FinmixModel.h"
#include "FinmixPrior.h"
#include "FinmixMCMC.h"


template <typename PriorType, typename ParType, typename LogType, 
	typename ParOutType> 
class FIX {
	public:
		class Node {
			public:
				const unsigned int K;
				const unsigned int N;
				const unsigned int M;
				const unsigned int BURNIN;
				const unsigned int STORES;
				const bool INDICFIX;
				const bool STARTPAR;
				const bool HIER;
				const bool RANPERM;
				const bool STOREPOST;
				PriorType hyperPar;
				ParType par;
				LogType log;
				const arma::mat y;
				arma::ivec S;
				const arma::mat expos;
				arma::urowvec compIndex;
				arma::urowvec permIndex;
				arma::urowvec compIndex2;

				Node (const FinmixData&, const FinmixModel&,
					const FinmixPrior&, const FinmixMCMC&);
				virtual void update ();
		};
		class Output {
			public: 
				const unsigned int M;
				const bool RANPERM;
				ParOutType par;
				arma::vec* mixlik;
				arma::vec* mixprior;
		
				Output (const Rcpp::S4&);
				virtual void store (const unsigned int&, Node&);
		};
		Node node;
		Output output;

		FIX (const FinmixData&, const FinmixModel&, const FinmixPrior&,
			const FinmixMCMC&, const Rcpp::S4&);
		virtual ~FIX () {}
		virtual void update ();
		virtual void store (const unsigned int&);
};

template <typename PriorType, typename ParType, typename LogType,
	typename ParOutType>
FIX <PriorType, ParType, LogType, ParOutType>::Node::Node (const FinmixData& data,
		const FinmixModel& model, const FinmixPrior& prior,
		const FinmixMCMC& mcmc) : 
		K(model.K), N(data.N), M(mcmc.M), BURNIN(mcmc.burnIn),
		STORES(mcmc.storeS), INDICFIX(model.indicFix), 
		STARTPAR(mcmc.startPar), HIER(prior.hier), 
		RANPERM(mcmc.ranPerm), STOREPOST(mcmc.storePost),
		hyperPar(prior), par(mcmc.startPar, model), log(), 
		y(data.y), S(data.S), expos(data.expos), compIndex(model.K), 
		permIndex(model.K), compIndex2(model.K)
{
        for(unsigned int k = 0; k < K; ++k) {
		compIndex(k) = k;
	}
}

template <typename PriorType, typename ParType, typename LogType,
	typename ParOutType>
void FIX <PriorType, ParType, LogType, ParOutType>::Node::update () 
{
	hyperPar.update(K, y, S, par);
	par.update(hyperPar);
	hyperPar.updateHier(par);
	log.update(K, y, S, expos, par, hyperPar);
	if(RANPERM && K > 1) {
		permIndex 	= arma::shuffle(compIndex, 1);
		compIndex2  	= (permIndex == compIndex);
		if(arma::sum(compIndex) != K) {
			par.lambda(compIndex) = par.lambda(permIndex);
		} 	
	}
}

template <typename PriorType, typename ParType, typename LogType, typename ParOutType>
FIX <PriorType, ParType, LogType, ParOutType>::Output::Output (const Rcpp::S4& classS4) : 
		M(Rcpp::as<unsigned int>((SEXP) classS4.slot("M"))),
		RANPERM(Rcpp::as<bool>((SEXP) classS4.slot("ranperm"))),
		par(Rcpp::as<Rcpp::List>((SEXP) classS4.slot("par"))) 
{
	Rcpp::List tmpLog((SEXP) classS4.slot("log"));
	Rcpp::NumericVector tmpMixLik((SEXP) tmpLog["mixlik"]);
	Rcpp::NumericVector tmpMixPrior((SEXP) tmpLog["mixprior"]);
	mixlik = new arma::vec(tmpMixLik.begin(), M, false, true);
	mixprior = new arma::vec(tmpMixPrior.begin(), M, false, true);
}

template <typename PriorType, typename ParType, typename LogType,
	typename ParOutType>
void FIX <PriorType, ParType, LogType, ParOutType>::Output::store (const unsigned
	int& m, Node& node)
{
	if (m >= node.BURNIN) {
		const unsigned int index = m - node.BURNIN;
		(*mixlik)(index) = node.log.mixlik;
		(*mixprior)(index) = node.log.mixprior;
		par.store(index, node.par);
	}
}

template <typename PriorType, typename ParType, typename LogType,
	typename ParOutType>
FIX <PriorType, ParType, LogType, ParOutType>::FIX (const FinmixData& data,
	const FinmixModel& model, const FinmixPrior& prior, const FinmixMCMC&
	mcmc, const Rcpp::S4& classS4) : 
		node(data, model, prior, mcmc), output(classS4) {}

template <typename PriorType, typename ParType, typename LogType, 
	typename ParOutType>
void FIX <PriorType, ParType, LogType, ParOutType>::update () 
{
	node.update();
}

template <typename PriorType, typename ParType, typename LogType,
	typename ParOutType>
void FIX <PriorType, ParType, LogType, ParOutType>::store (const unsigned int& m) 
{
	output.store(m, node);
}
#endif
