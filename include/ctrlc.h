/*
 * Interface for SCS signal handling.
 */

#ifndef __CTRLC_H__
#define __CTRLC_H__

#if CTRLC > 0

#if defined MATLAB_MEX_FILE

/* No header file available here; define the prototypes ourselves */
extern int utIsInterruptPending();
extern int utSetInterruptEnabled(int);

#elif (defined _WIN32 || defined _WIN64 || defined _WINDLL )

/* Use Windows SetConsoleCtrlHandler for signal handling */
#include <Windows.h>

#else

/* Use POSIX clocl_gettime() for timing on non-Windows machines */
#include <signal.h>

#endif

/* METHODS are the same for both */
void init_ctrlc(void);
void remove_ctrlc(void);
int check_ctrlc(void);

#else /* CTRLC = 0 */

/* No signal handling. */
#define init_ctrlc()
#define remove_ctrlc()
#define check_ctrlc() 0

#endif /* END IF CTRLC > 0 */

#endif /* END IFDEF __TIMER_H__ */

