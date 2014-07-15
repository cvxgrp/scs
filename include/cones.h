#ifndef CONES_H_GUARD
#define CONES_H_GUARD

#include "scs.h"

struct CONE {
	idxint f; /* number of linear equality constraints */
	idxint l; /* length of LP cone */
	idxint *q; /* array of second-order cone constraints */
	idxint qsize; /* length of SOC array */
	idxint *s; /* array of SD constraints */
	idxint ssize; /* length of SD array */
	idxint ep; /* number of primal exponential cone triples */
	idxint ed; /* number of dual exponential cone triples */
};

/*
 * boundaries will contain array of indices of rows of A corresponding to
 * cone boundaries, boundaries[0] is starting index for cones of size larger than 1
 * returns length of boundaries array, boundaries malloc-ed here so should be freed
 */
idxint getConeBoundaries(Cone * k, idxint ** boundaries);

idxint initCone(Cone * k);
char * getConeHeader(Cone * k);
idxint validateCones(Data * d, Cone * k);
/* pass in iter to control how accurate the cone projection
 with iteration, set iter < 0 for exact projection, warm_start contains guess
 of solution, can be NULL*/
idxint projDualCone(pfloat *x, Cone *k, const pfloat * warm_start, idxint iter);
void finishCone(void);
char * getConeSummary(Info * info);

#endif
