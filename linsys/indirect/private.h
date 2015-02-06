#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#include "glbopts.h"
#include "scs.h"
#include "../common.h"

struct PRIVATE_DATA{
	scs_float * p;
	scs_float * r;
	scs_float * Gp;
	scs_float * M; /* pre-conditioner */
	scs_float * z;
   	scs_float * G; /* Gram matrix */
};

#endif
