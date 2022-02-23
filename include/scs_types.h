/*
 * Definitions of primitive types used in SCS.
 */

#ifndef SCS_TYPES_H_GUARD
#define SCS_TYPES_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#ifdef DLONG
/*#ifdef _WIN64
#include <stdint.h>
typedef int64_t scs_int;
#else
typedef long scs_int;
#endif
*/
typedef long long scs_int;
#else
typedef int scs_int;
#endif

#ifndef SFLOAT
typedef double scs_float;
#else
typedef float scs_float;
#endif

#ifdef __cplusplus
}
#endif
#endif
