/*
 * Interface for SCS signal handling.
 */

#ifndef CTRLC_H_GUARD
#define CTRLC_H_GUARD

#if CTRLC > 0

#if defined MATLAB_MEX_FILE

/* No header file available here; define the prototypes ourselves */
extern int utIsInterruptPending();
extern int utSetInterruptEnabled(int);

#elif(defined _WIN32 || defined _WIN64 || defined _WINDLL)

/* Use Windows SetConsoleCtrlHandler for signal handling */
#include <windows.h>

#else

/* Use POSIX clocl_gettime() for timing on non-Windows machines */
#include <signal.h>

#endif

/* METHODS are the same for both */
void startInterruptListener(void);
void endInterruptListener(void);
int isInterrupted(void);

#else /* CTRLC = 0 */

/* No signal handling. */
#define startInterruptListener()
#define endInterruptListener()
#define isInterrupted() 0

#endif /* END IF CTRLC > 0 */

#endif /* END IFDEF __TIMER_H__ */
