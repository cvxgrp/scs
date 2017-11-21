#ifndef UTIL_H_GUARD
#define UTIL_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "scs.h"
#include <stdlib.h>
#include <stdio.h>
#include "cones.h"

/* timing code courtesty of A. Domahidi */
#if (defined NOTIMER)
typedef void *timer;
#elif(defined _WIN32 || defined _WIN64 || defined _WINDLL)
/* Use Windows QueryPerformanceCounter for timing */
#include <windows.h>
typedef struct timer {
  LARGE_INTEGER tic;
  LARGE_INTEGER toc;
  LARGE_INTEGER freq;
} timer;

#elif(defined __APPLE__)
/* Use MAC OSX mach_time for timing */
#include <mach/mach_time.h>
typedef struct timer {
  uint64_t tic;
  uint64_t toc;
  mach_timebase_info_data_t tinfo;
} timer;

#else
/* Use POSIX clock_gettime() for timing on other machines */
#include <time.h>
typedef struct timer {
  struct timespec tic;
  struct timespec toc;
} timer;

#endif

#if EXTRA_VERBOSE > 1
extern timer global_timer;
#endif

/* these all return milli-seconds */
void scs_tic(timer *t);
scs_float scs_toc(timer *t);
scs_float scs_str_toc(char *str, timer *t);
scs_float tocq(timer *t);

void print_cone_data(const ScsCone *k);
void print_data(const ScsData *d);
void print_work(const ScsWork *w);
void print_array(const scs_float *arr, scs_int n, const char *name);
void set_default_scs_settings(ScsData *d);
void free_sol(ScsSolution *sol);
void free_data(ScsData *d, ScsCone *k);

#ifdef __cplusplus
}
#endif
#endif
