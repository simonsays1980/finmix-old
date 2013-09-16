/******************************************************************************
 *
 * TODO: Project Title
 *
 * Copyright (C) 2012-2013 Lars Simon Zehnder. All Rights Reserved.
 * Web: -
 *
 * Author: Lars Simon Zehnder <simon.zehnder@gmail.com>
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 ******************************************************************************/

#ifndef __OPTIMIZE_H__
#define __OPTIMIZE_H__

#include <RcppArmadillo.h>
#include <vector>
#include "distributions.h"

inline
double obj_stephens1997a_poisson (const std::vector<double> &x, 
        std::vector<double> &grad, void *f_data)
{
    std::vector<arma::mat*> *arma_data = static_cast<std::vector<arma::mat*>* >(f_data);
    const unsigned int M = (*arma_data)[0]->n_rows; 
    const unsigned int K = (*arma_data)[0]->n_cols;
    arma::vec rvalues(M);
    arma::vec arma_x(x);
    arma::vec dirich(&x[0], K);
    arma::vec shape(&x[0] + K, K);
    arma::vec rate(&x[0] + 2 * K, K);
    rvalues = arma::exp(lddirichlet((*(*arma_data)[1]), dirich));
    rvalues = rvalues % arma::prod(arma::exp(ldgamma((*(*arma_data)[0]), 
                shape, rate)), 1);
    return arma::sum(arma::log(rvalues));
}

#endif /* __OPTIMIZE_H__ */



