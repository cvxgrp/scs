#ifndef UTIL_H_GUARD
#define UTIL_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include "scs.h"
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

/** \brief SCS timer timer timer */
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
    void freeInfo(Info *info);

    /**
     * \brief Custom print function for SCS.
     * 
     * This functions allows to print in different streams. The argument \c print_mode
     * specifies whether it is allowed to override the default stream and use a 
     * print function other than \c printf. 
     * 
     * For example, if SCS is interfaced via a MEX function, MATLAB expects to 
     * use \c printf exclusively which it then delegated to \c mexPrintf; a function
     * that prints the program's output to the console of MATLAB.
     * 
     * When SCS is called from other software, it is likely that \c print_mode 
     * has to be set to 0.
     * 
     * @param print_mode whether to override the default behavior (using \c printf)
     * @param __stream an output stream
     * @param __format string format 
     * @param ... arguments specifying data to print
     * @return return value of \c print or \c vfprintf
     */
    int scs_special_print(scs_int print_mode,
            FILE *__restrict __stream,
            const char *__restrict __format, ...);

#ifdef __cplusplus
}
#endif
#endif
