#ifndef UTIL_H_GUARD
#define UTIL_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "cones.h"
#include "scs.h"

/* timing code courtesy of A. Domahidi */
#if (defined NOTIMER)
typedef void *SCS(timer);
#elif (defined _WIN32 || defined _WIN64 || defined _WINDLL)
/* Use Windows QueryPerformanceCounter for timing */
#include <windows.h>
typedef struct SCS(timer) {
  LARGE_INTEGER tic;
  LARGE_INTEGER toc;
  LARGE_INTEGER freq;
} SCS(timer);

#elif (defined __APPLE__)
/* Use MAC OSX mach_time for timing */
#include <mach/mach_time.h>
typedef struct SCS(timer) {
  uint64_t tic;
  uint64_t toc;
  mach_timebase_info_data_t tinfo;
} SCS(timer);

#else
/* Use POSIX clock_gettime() for timing on other machines */
#include <time.h>
typedef struct SCS(timer) {
  struct timespec tic;
  struct timespec toc;
} SCS(timer);

#endif

#if EXTRA_VERBOSE > 1
extern SCS(timer) global_timer;
#endif

/* these all return milli-seconds */
void SCS(tic)(SCS(timer) * t);
scs_float SCS(toc)(SCS(timer) * t);
scs_float SCS(str_toc)(char *str, SCS(timer) * t);
scs_float SCS(tocq)(SCS(timer) * t);

void SCS(print_cone_data)(const ScsCone *k);
void SCS(print_data)(const ScsData *d);
void SCS(print_work)(const ScsWork *w);
void SCS(print_array)(const scs_float *arr, scs_int n, const char *name);
void SCS(set_default_settings)(ScsData *d);
void SCS(free_sol)(ScsSolution *sol);
void SCS(free_data)(ScsData *d, ScsCone *k);

#ifdef __cplusplus
}
#endif
#endif
