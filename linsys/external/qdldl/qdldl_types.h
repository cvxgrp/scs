#ifndef QDLDL_TYPES_H
# define QDLDL_TYPES_H

# ifdef __cplusplus
extern "C" {
# endif /* ifdef __cplusplus */

#include "glbopts.h"
#include <limits.h> //for the QDLDL_INT_TYPE_MAX

/* QDLDL integer and float types */

#define QDLDL_int scs_int
#define QDLDL_float scs_float
#define QDLDL_bool scs_int

/* Maximum value of the signed type QDLDL_int */
#ifdef DLONG
#define QDLDL_INT_MAX LLONG_MAX
#else
#define QDLDL_INT_MAX INT_MAX 
#endif

# ifdef __cplusplus
}
# endif /* ifdef __cplusplus */

#endif /* ifndef QDLDL_TYPES_H */

