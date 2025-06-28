#ifndef SCS_TYPES_H_GUARD
#define SCS_TYPES_H_GUARD

#ifdef __cplusplus
#  include <complex>   // Use C++ complex type
extern "C" {
#else
#  include <complex.h> // Use C complex.h
#endif

#ifdef DLONG
typedef long long scs_int;
#else
typedef int scs_int;
#endif

#ifndef SFLOAT
typedef double scs_float;
#  ifdef __cplusplus
typedef std::complex<scs_float> scs_complex_float;
#  else
typedef double _Complex scs_complex_float;
#  endif
#else
typedef float scs_float;
#  ifdef __cplusplus
typedef std::complex<scs_float> scs_complex_float;
#  else
typedef float _Complex scs_complex_float;
#  endif
#endif

#ifdef __cplusplus
}
#endif
#endif
