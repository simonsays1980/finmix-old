#ifndef __FINMIX_PRIORSTUDENTIND_H__
#define __FINMIX_PRIORSTUDENTIND_H__

#include "PriorStudentFix.h"

/* Forward declaration */
class ParStudentInd;
class PriorStudentInd : virtual public PriorStudentFix {
    public:
        arma::rowvec weightStart;
        arma::rowvec weightPost;

        PriorStudentInd (const FinmixPrior&);
        virtual ~PriorStudentInd () {}
        virtual void update (const unsigned int&, 
                const arma::mat&, arma::ivec&,
                const arma::vec&, ParStudentInd&);
        virtual void updateDf (const unsigned int&, const arma::mat&, 
                const arma::ivec&, ParStudentInd&);
};
#endif /* __FINMIX_PRIORSTUDENTIND_H__ */



