/*
 * Implements signal handling (ctrl-c) for SCS.
 *
 * Under Windows, we use SetConsoleCtrlHandler.
 * Under Unix systems, we use sigaction.
 * For Mex files, we use utSetInterruptEnabled/utIsInterruptPending.
 *
 */

#include "ctrlc.h"

#if CTRLC > 0

#ifdef MATLAB_MEX_FILE

static int istate;
void startInterruptListener(void) {
    istate = utSetInterruptEnabled(1);
}

void endInterruptListener(void) {
    utSetInterruptEnabled(istate);
}

int isInterrupted(void) {
    return utIsInterruptPending();
}

#elif defined _WIN32 || defined _WIN64

static int int_detected;
BOOL WINAPI handle_ctrlc(DWORD dwCtrlType) {
    if (dwCtrlType != CTRL_C_EVENT)
        return FALSE;
    int_detected = 1;
    return TRUE;
}

void startInterruptListener(void) {
    int_detected = 0;
    SetConsoleCtrlHandler(handle_ctrlc, TRUE);
}

void endInterruptListener(void) {
    SetConsoleCtrlHandler(handle_ctrlc, FALSE);
}

int isInterrupted(void) {
    return int_detected;
}

#else /* Unix */

#include <signal.h>
static int int_detected;
struct sigaction oact;
void handle_ctrlc(int dummy) {
    int_detected = dummy ? dummy : -1;
}

void startInterruptListener(void) {
    struct sigaction act;
    int_detected = 0;
    act.sa_flags = 0;
    sigemptyset(&act.sa_mask);
    act.sa_handler = handle_ctrlc;
    sigaction(SIGINT, &act, &oact);
}

void endInterruptListener(void) {
    struct sigaction act;
    sigaction(SIGINT, &oact, &act);
}

int isInterrupted(void) {
    return int_detected;
}

#endif /* END IF MATLAB_MEX_FILE / WIN32 */

#endif /* END IF CTRLC > 0 */
