#ifndef CONES_H_GUARD                                                              
#define CONES_H_GUARD

#include <string.h>
#include <math.h>
#include "glbopts.h"
#include "scs.h"

struct CONE {
    idxint f;           /* number of linear equality constraints */
    idxint l;           /* length of LP cone */
    idxint *q;   	    /* array of second-order cone constraints */
    idxint qsize;       /* length of SOC array */
	idxint *s;			/* array of SD constraints */
	idxint ssize;		/* length of SD array */
    idxint ep;          /* number of primal exponential cone triples */
    idxint ed;          /* number of dual exponential cone triples */
};

idxint initCone(Cone * k);
void finishCone();
void projCone(pfloat *x, Cone *k, idxint iter);
idxint getFullConeDims(Cone * k);
char * getConeHeader(Cone * k);
idxint validateCones(Cone * k);
#endif
