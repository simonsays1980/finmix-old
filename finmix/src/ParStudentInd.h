#ifndef __FINMIX_PARSTUDENTIND_H__
#define __FINMIX_PARSTUDENTIND_H__

#include "ParStudentFix.h"
#include "PriorStudentInd.h"

class ParStudentInd : virtual public ParStudentFix {
    public:
        arma::rowvec weight;

        ParStudentInd (const bool&, 
                const FinmixModel&);
        virtual ~ParStudentInd () {}
        virtual void update (const PriorStudentInd&);
};
#endif /* __FINMIX_PARSTUDENTIND_H__ */



