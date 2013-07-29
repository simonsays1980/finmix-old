/******************************************************************************
 *
 * Copyright (C) 2013 Lars Simon Zehnder. All Rights Reserved.
 *
 * Author: Lars Simon Zehnder <simon.zehnder@gmail.com>
 *
 * This file is part of the R package finmix.
 *
 * finmix is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundatio, either version 3 of the License, or
 * any later version.
 *
 * finmix is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with 'finmix'. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef IND_H
#define IND_H

#include <RcppArmadillo.h>

template <typename Super> 
class IND : public Super {
	public:
		class Node : public Super::Node {
			public: 
				arma::urowvec swapIndex;

				Node (const FinmixData&,
					const FinmixModel&,
					const FinmixPrior&,
					const FinmixMCMC&);
				virtual ~Node () {}
				virtual void update (); 
		};
		class Output : public Super::Output {
			public:
				arma::mat* weight;
				arma::vec* cdpost;
				arma::vec* entropy;
				arma::ivec* ST;
				arma::imat* S;
				arma::imat* NK;
				arma::ivec* clust;
				Output (const Rcpp::S4&);
				virtual ~Output () {}
				virtual void store (const unsigned int&,
					Node&);
		};
		Node node;
		Output output;

		IND (const FinmixData&, const FinmixModel&,
			const FinmixPrior&, const FinmixMCMC&,
			const Rcpp::S4&);
		virtual ~IND () {}
		virtual void update ();
		virtual void store (const unsigned int&);
};

template <typename Super>
IND <Super>::Node::Node (const FinmixData& data, 
	const FinmixModel& model, const FinmixPrior& prior,
	const FinmixMCMC& mcmc) :
		Super::Node(data, model, prior, mcmc),
		swapIndex(model.K) {}

template <typename Super>
void IND <Super>::Node::update () 
{	
	Super::Node::update();
	if (Super::Node::RANPERM && arma::sum(Super::Node::compIndex2) !=
		Super::Node::K) {
		Super::Node::par.weight(Super::Node::compIndex) = 
			Super::Node::par.weight(Super::Node::permIndex);
		swapIndex = arma::sort_index(Super::Node::permIndex);
		for(unsigned int i = 0; i < Super::Node::N; ++i) {
			Super::Node::S(i) = (int) swapIndex((unsigned int) 
				(Super::Node::S(i) - 1)) + 1;
		}	
	}	
}

template <typename Super>
IND <Super>::Output::Output (const Rcpp::S4& classS4) : 
	Super::Output(classS4) 
{
	Rcpp::NumericMatrix tmpWeight((SEXP) classS4.slot("weight"));
	Rcpp::List tmpList((SEXP) classS4.slot("log"));
	Rcpp::NumericVector tmpCDPost((SEXP) tmpList["cdpost"]);
	Rcpp::NumericVector tmpEntropy((SEXP) classS4.slot("entropy"));
	Rcpp::IntegerVector tmpST((SEXP) classS4.slot("ST"));
	Rcpp::IntegerMatrix tmpS((SEXP) classS4.slot("S"));
	Rcpp::IntegerMatrix tmpNK((SEXP) classS4.slot("NK"));
	Rcpp::IntegerVector tmpClust((SEXP) classS4.slot("clust"));
	const unsigned int tmpM = tmpWeight.nrow();
	const unsigned int K = tmpWeight.ncol();
	const unsigned int N = tmpS.nrow();
	const unsigned int STORES = tmpS.ncol();
	weight = new arma::mat(tmpWeight.begin(), tmpM, K, false, true);
	cdpost = new arma::vec(tmpCDPost.begin(), tmpM, false, true);
	entropy = new arma::vec(tmpEntropy.begin(), tmpM, false, true);
	ST = new arma::ivec(tmpST.begin(), tmpM, false, true);
	S = new arma::imat(tmpS.begin(), N, STORES, false, true);
	NK = new arma::imat(tmpNK.begin(), tmpM, K, false, true);
	clust = new arma::ivec(tmpClust.begin(), N, false, true);
}

template <typename Super> 
void IND <Super>::Output::store (const unsigned int& m,
	Node& node)
{
	Super::Output::store(m,node);
	if(m >= node.BURNIN) {
		const unsigned int index = m - node.BURNIN;
		(*weight).row(index) = node.par.weight;
		(*cdpost)(index) = node.log.cdpost;
		(*entropy)(index) = node.log.entropy;
		(*ST)(index) = node.S(node.N - 1);
		if(index >= node.M - node.STORES) {
			if(node.STARTPAR && index != node.M - 1) {
				(*S).col(index - (node.M - node.STORES) + 1) = node.S;
			}
			if (!node.STARTPAR){
				(*S).col(index - (node.M - node.STORES)) = node.S;
			}
		}
		(*NK).row(index) = arma::conv_to<arma::irowvec>::from
			(node.hyperPar.weightPost - node.hyperPar.weightStart);
		if(m == node.BURNIN) {
			node.log.maxcdpost = node.log.cdpost - 1;		
		}
		if(node.log.cdpost > node.log.maxcdpost) {
			(*clust) = node.S;
		}
	}
}

template <typename Super>
IND <Super>::IND (const FinmixData& data, const FinmixModel& model,
	const FinmixPrior& prior, const FinmixMCMC& mcmc,
	const Rcpp::S4& classS4) :
		Super(data, model, prior, mcmc, classS4),
		node(data, model, prior, mcmc),
		output(classS4) {}

template <typename Super>
void IND <Super>::update ()
{
	node.update();
}

template <typename Super>
void IND <Super>::store (const unsigned int& m)
{
	output.store(m, node);
}  
#endif
