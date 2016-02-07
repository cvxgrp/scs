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

#if EXTRAVERBOSE > 1
extern timer globalTimer;
#endif

/* these all return milli-seconds */
void tic(timer *t);
scs_float toc(timer *t);
scs_float strtoc(char *str, timer *t);
scs_float tocq(timer *t);

void printConeData(const Cone *k);
void printData(const Data *d);
void printWork(const Work *w);
void printArray(const scs_float *arr, scs_int n, const char *name);
void setDefaultSettings(Data *d);
void freeSol(Sol *sol);
void freeData(Data *d, Cone *k);

#ifdef __cplusplus
}
#endif
#endif
