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
#ifndef __FINMIX_LOGBINOMIALIND_H__
#define __FINMIX_LOGBINOMIALIND_H__

#include "LogBinomialFix.h"
#include "ParBinomialInd.h"
#include "PriorBinomialInd.h"

class LogBinomialInd : public LogBinomialFix {
    public:
        double cdpost;
        double entropy;
        double maxcdpost;

        LogBinomialInd ();
        virtual ~LogBinomialInd () {}
        void update (const unsigned int&, const arma::mat&,
                arma::ivec&, const arma::mat&, const arma::vec&,
                const ParBinomialInd&, const PriorBinomialInd&);
};
#endif /* __FINMIX_LOGBINOMIALIND_H__ */



