/*
 * Interface for SCS signal handling.
 */

#ifndef CTRLC_H_GUARD
#define CTRLC_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#if CTRLC > 0

#if defined MATLAB_MEX_FILE

/* No header file available here; define the prototypes ourselves */
extern int utIsInterruptPending();
extern int utSetInterruptEnabled(int);

#elif (defined _WIN32 || defined _WIN64 || defined _WINDLL)

/* Use Windows set_console_ctrl_handler for signal handling */
#include <windows.h>

#else

/* Use POSIX clocl_gettime() for timing on non-Windows machines */
#include <signal.h>

#endif

/* METHODS are the same for both */
void scs_start_interrupt_listener(void);
void scs_end_interrupt_listener(void);
int scs_is_interrupted(void);

#else /* CTRLC = 0 */

/* Simply to suppress empty translation unit warnings. */
typedef int scs_make_iso_compilers_happy;

/* No signal handling. */
#define scs_start_interrupt_listener()
#define scs_end_interrupt_listener()
#define scs_is_interrupted() 0

#endif /* END IF CTRLC > 0 */

#ifdef __cplusplus
}
#endif
#endif
