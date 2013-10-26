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
 * by the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * finmix is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with finmix. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef __FINMIX_PAROUTBINOMIAL_H__
#define __FINMIX_PAROUTBINOMIAL_H__

#include "ParBinomialFix.h"

class ParOutBinomial {
    public:
        arma::mat* p;

        ParOutBinomial () {}
        ParOutBinomial (const Rcpp::List&);
        void store(const unsigned int&, const ParBinomialFix&);
};
#endif /* __FINMIX_PAROUTBINOMIAL_H__ */



