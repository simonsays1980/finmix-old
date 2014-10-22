/******************************************************************************
 *
 * TODO: Project Title
 *
 * Copyright (C) 2003-2009 ascolab GmbH. All Rights Reserved.
 * Web: http://www.ascolab.com
 *
 * Author: Gerhard Gappmeier <gerhard.gappmeier@ascolab.com>
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 ******************************************************************************/

#ifndef __FINMIX_PRIORSTUDMULTIND_H__
#define __FINMIX_PRIORSTUDMULTIND_H__

#include "PriorStudmultFix.h"

/* Forward declaration */
class ParStudmultInd;

class PriorStudmultInd : virtual public PriorStudmultFix {
    public: 
        arma::rowvec weightStart;
        arma::rowvec weightPost;

        PriorStudmultInd (const FinmixPrior&);
        virtual ~PriorStudmultInd () {}
        virtual void update (const unsigned int&, 
                const arma::mat&, arma::ivec&,
                const arma::vec&, ParStudmultInd&);
        virtual void updateDf(const unsigned int&, 
                const arma::mat&, const arma::ivec&,
                ParStudmultInd&);
};
#endif /* __FINMIX_PRIORSTUDMULTIND_H__ */



