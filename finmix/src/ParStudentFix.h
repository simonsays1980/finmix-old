#ifndef __FINMIX_PARSTUDENTFIX_H__
#define __FINMIX_PARSTUDENTFIX_H__

#include "FinmixModel.h"
#include "PriorStudentFix.h"

class ParStudentFix {
    public:
        arma::rowvec mu;
        arma::rowvec sigma;
        arma::rowvec df;
        arma::rowvec acc;
        bool INDEPENDENT;

        ParStudentFix (const bool&, const FinmixModel&);
        virtual ~ParStudentFix () {}
        virtual void update (const PriorStudentFix&);
        virtual void permute (const arma::urowvec&, 
                const arma::urowvec&);
};
#endif /* __FINMIX_PARSTUDENTFIX_H__ */



