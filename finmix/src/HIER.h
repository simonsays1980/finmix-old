#ifndef HIER_H
#define HIER_H

#include <RcppArmadillo.h>

template <typename Super, typename HierOutType>
class HIER : public Super {
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
				HierOutType hyper;

				Output (const Rcpp::S4&);
				virtual ~Output () {}
				virtual void store (const 
						unsigned int&, 
						Node&);
		};
		Node node;
		Output output;
		
		HIER (const FinmixData&, const FinmixModel&,
			const FinmixPrior&, const FinmixMCMC&,
			const Rcpp::S4&);
		virtual ~HIER () {}
		virtual void update ();
		virtual void store (const unsigned int&);
};

template <typename Super, typename HierOutType>
HIER <Super, HierOutType>::Node::Node (const FinmixData& data, 
	const FinmixModel& model, const FinmixPrior& prior,
	const FinmixMCMC& mcmc) : 
		Super::Node(data, model, prior, mcmc) {}

template <typename Super, typename HierOutType> 
HIER <Super, HierOutType>::Output::Output (const Rcpp::S4& classS4) :
	Super::Output(classS4),
	hyper(Rcpp::as<Rcpp::List>((SEXP) classS4.slot("hyper"))) {}

template <typename Super, typename HierOutType>
void HIER <Super, HierOutType>::Output::store (const unsigned int& m, 
	Node& node)
{
	Super::Output::store(m, node);
	if (m >= node.BURNIN) {
		const unsigned int index = m - node.BURNIN;
		hyper.store(index, node.hyperPar);
	}
}

template <typename Super, typename HierOutType>
HIER <Super, HierOutType>::HIER (const FinmixData& data, 
	const FinmixModel& model, const FinmixPrior& prior,
	const FinmixMCMC& mcmc, const Rcpp::S4& classS4) :
		Super(data, model, prior, mcmc, classS4),
		node(data, model, prior, mcmc),
		output(classS4) {}

template <typename Super, typename HierOutType> 
void HIER <Super, HierOutType>::update () 
{
	node.update();
}

template <typename Super, typename HierOutType>
void HIER <Super, HierOutType>::store (const unsigned int& m) 
{
	output.store (m, node);
}
#endif
