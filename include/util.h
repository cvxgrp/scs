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

    /**
     * \brief Starts the timer
     * @param t timer
     */
    void tic(timer *t);
    /**
     * \brief Stops the timer 
     * 
     * This function prints a message of the form: 
     * 
     * ~~~~~
     * time: %8.4f milli-seconds.\n
     * ~~~~~
     * 
     * to the standard output using \c printf.
     * 
     * @param t timer 
     * @return elapsed time in milliseconds
     */
    scs_float toc(timer *t);
    /**
     * \brief Stops the timer and prints a custom message
     * 
     * 
     * @param str string
     * @param t timer 
     * @return elapsed time
     */
    scs_float strtoc(char *str, timer *t);
    /**
     * \brief Stops the timer 
     * @param t timer 
     * @return elapsed time in milliseconds
     */
    scs_float tocq(timer *t);

    /**
     * \brief Prints the content of a Cone object
     * @param k pointer to cone
     */
    void printConeData(const Cone *k);
    /**
     * \brief Prints the content of a Data object
     * @param d pointer to data
     */
    void printData(const Data *d);
    /**
     * \brief Prints the content of a Work object
     * @param w pointer to work
     */
    void printWork(const Work *w);

    /**
     * \brief Prints an array
     * 
     * @param arr pointer to array
     * @param n length of array
     * @param name name of the array
     */
    void printArray(const scs_float *arr, scs_int n, const char *name);

    /**
     * \brief Sets the settings to certain default values
     * 
     * <table>
     * <tr><th>Parameter<th>Default value
     * <tr><td>\ref Settings::normalize "normalize"<td>1
     * <tr><td>\ref Settings::scale "scale"<td>1.0
     * <tr><td>\ref Settings::rho_x "rho_x"<td>1.0
     * <tr><td>\ref Settings::max_iters "max_iters"<td>2500
     * <tr><td>\ref Settings::previous_max_iters "previous_max_iters"<td>-1
     * <tr><td>\ref Settings::eps "eps"<td>1e-3
     * <tr><td>\ref Settings::alpha "alpha"<td>1.5
     * <tr><td>\ref Settings::cg_rate "cg_rate"<td>2.0
     * <tr><td>\ref Settings::verbose "verbose"<td>1
     * <tr><td>\ref Settings::warm_start "warm_start"<td>0
     * <tr><td>\ref Settings::do_super_scs "do_super_scs"<td>1
     * <tr><td>\ref Settings::k0 "k0"<td>0
     * <tr><td>\ref Settings::k1 "k1"<td>1
     * <tr><td>\ref Settings::k2 "k2"<td>1
     * <tr><td>\ref Settings::c_bl "c_bl"<td>0.999
     * <tr><td>\ref Settings::c_bl "c1"<td>0.9999
     * <tr><td>\ref Settings::ls "ls"<td>10
     * <tr><td>\ref Settings::beta "beta"<td>0.5
     * <tr><td>\ref Settings::sigma "sigma"<td>0.01
     * <tr><td>\ref Settings::direction "direction"<td>\ref restarted_broyden "restarted_broyden"
     * <tr><td>\ref Settings::thetabar "thetabar"<td>0.1
     * <tr><td>\ref Settings::memory "memory"<td>10
     * <tr><td>\ref Settings::tRule "tRule"<td>1
     * <tr><td>\ref Settings::broyden_init_scaling "broyden_init_scaling"<td>1
     * <tr><td>\ref Settings::do_record_progress "do_record_progress"<td>0
     * <tr><td>\ref Settings::do_override_streams "do_override_streams"<td>0
     * <tr><td>\ref Settings::output_stream "output_stream"<td>\c stdout
     * </table>
     * 
     * @param d Pointer to data
     * 
     * \warning If you want to increase the maximum number of iteration with respect
     * to the previous run and you have set \ref Settings::do_record_progress "do_record_progress"
     * to \c 1, then you should not use this function. If you really want to use it, however,
     * you should set the parameter \ref Settings::previous_max_iters "previous_max_iters"
     * to the maximum number of iterations you used in the previous run. This is in 
     * order to avoid memory management errors. 
     * 
     * \warning Alternatively, a simple solution is to invoke ::freeInfo after you 
     * call ::scs and then again ::initInfo. Then it is safe to call this function and
     * run ::scs again.
     * 
     * \note If you have set \ref Settings::do_record_progress "do_record_progress" to \c 0,
     * you may ignote this warning.
     */
    void setDefaultSettings(Data *d);

    /**
     * \brief Frees the memory allocated for a Sol object
     * 
     * @param sol
     * 
     * \sa initSol
     */
    void freeSol(Sol *sol);
    /**
     * \brief Frees the memory allocate of a Data and a Cone object
     * @param d
     * @param k
     * 
     * \sa setDefaultSettings
     * \sa initData
     */
    void freeData(Data *d, Cone *k);
    /**
     * \brief Frees the memory allocated for an Info object
     * @param info
     * 
     * \sa ::initInfo
     */
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
