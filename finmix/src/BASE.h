#ifndef BASE_H
#define BASE_H

#include <RcppArmadillo.h>
 
class BASE {
	public:
		BASE () {}
		virtual ~BASE () {}
		virtual void update () {}
		virtual void store (const unsigned int&) {}
};
#endif
