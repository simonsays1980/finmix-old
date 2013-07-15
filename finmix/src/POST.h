#ifndef POST_H
#define POST_H

#include <RcppArmadillo.h>

template <typename Super, typename PostOutType>
class POST : public Super {
	public:
		class Node : public Super::Node {
			public:
				Node (const FinmixData&,
					const FinmixModel&,
					const FinmixPrior&,
					const FinmixMCMC&);
				virtual ~Node () {}
		};
		class Output : public Super::Output {
			public:
				PostOutType post;
				
				Output (const Rcpp::S4&);
				virtual ~Output () {}	
				virtual void store (const
						unsigned int&,
						Node&);
		};	
		Node node;
		Output output;
		
		POST (const FinmixData&, const FinmixModel&,
			const FinmixPrior&, const FinmixMCMC&,
			const Rcpp::S4&);
		virtual ~POST () {}
		virtual void update ();
		virtual void store (const unsigned int&);
};

template <typename Super, typename PostOutType>
POST <Super, PostOutType>::Node::Node (const FinmixData& data,
	const FinmixModel& model, const FinmixPrior& prior,
	const FinmixMCMC& mcmc) :
		Super::Node(data, model, prior, mcmc) {}

template <typename Super, typename PostOutType> 
POST <Super, PostOutType>::Output::Output (const Rcpp::S4& classS4) :
	Super::Output(classS4),
	post(Rcpp::as<Rcpp::List>((SEXP) classS4.slot("post"))) {}

template <typename Super, typename PostOutType>
void POST <Super, PostOutType>::Output::store (const unsigned int& m, 
	Node& node)
{
	Super::Output::store(m, node);
	if(m >= node.BURNIN) {
		const unsigned int index = m - node.BURNIN;
		post.store(index, node.hyperPar);
	}
}

template <typename Super, typename PostOutType>
POST <Super, PostOutType>::POST (const FinmixData& data, 
	const FinmixModel& model, const FinmixPrior& prior, 
	const FinmixMCMC& mcmc, const Rcpp::S4& classS4) :
		Super(data, model, prior, mcmc, classS4),
		node(data, model, prior, mcmc),
		output(classS4) {}

template <typename Super, typename PostOutType>
void POST <Super, PostOutType>::update () 
{
	node.update();
} 

template <typename Super, typename PostOutType>
void POST <Super, PostOutType>::store (const unsigned int& m)
{
	output.store(m, node);
}
#endif
