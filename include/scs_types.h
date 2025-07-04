/*
 * Public header including definitions of primitive types used in SCS.
 * Make sure this file and `scs.h` are somewhere appropriate and then use
 * `#include "scs.h"` to access the SCS public API.
 */

#ifndef SCS_TYPES_H_GUARD
#define SCS_TYPES_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include <complex.h>

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
typedef double _Complex scs_complex_float;
#else
typedef float scs_float;
typedef float _Complex scs_complex_float;
#endif

#ifdef __cplusplus
}
#endif
#endif
