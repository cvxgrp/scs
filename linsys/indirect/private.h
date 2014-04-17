#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#include "glbopts.h"
#include "scs.h"
#include "cs.h"
#include <math.h>

struct PRIVATE_DATA {
	pfloat *p; /* cg iterate  */
	pfloat *r; /* cg residual */
	pfloat *x;
	pfloat * Ap;
	pfloat * tmp;
	pfloat * Atx;
	idxint * Ati;
	idxint * Atp;
	/* preconditioning */
	pfloat * z;
	pfloat * M;
};

#endif
