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
void start_interrupt_listener(void) {
  istate = utSetInterruptEnabled(1);
}

void end_interrupt_listener(void) {
  utSetInterruptEnabled(istate);
}

int is_interrupted(void) {
  return utIsInterruptPending();
}

#elif(defined _WIN32 || _WIN64 || defined _WINDLL)

static int int_detected;
BOOL WINAPI scs_handle_ctrlc(DWORD dwCtrlType) {
  if (dwCtrlType != CTRL_C_EVENT) {
    return FALSE;
  }
  int_detected = 1;
  return TRUE;
}

void start_interrupt_listener(void) {
  int_detected = 0;
  SetConsoleCtrlHandler(scs_handle_ctrlc, TRUE);
}

void end_interrupt_listener(void) {
  SetConsoleCtrlHandler(scs_handle_ctrlc, FALSE);
}

int is_interrupted(void) {
  return int_detected;
}

#else /* Unix */

#include <signal.h>
static int int_detected;
struct sigaction oact;
void scs_handle_ctrlc(int dummy) {
  int_detected = dummy ? dummy : -1;
}

void start_interrupt_listener(void) {
  struct sigaction act;
  int_detected = 0;
  act.sa_flags = 0;
  sigemptyset(&act.sa_mask);
  act.sa_handler = scs_handle_ctrlc;
  sigaction(SIGINT, &act, &oact);
}

void end_interrupt_listener(void) {
  struct sigaction act;
  sigaction(SIGINT, &oact, &act);
}

int is_interrupted(void) {
  return int_detected;
}

#endif /* END IF MATLAB_MEX_FILE / WIN32 */

#endif /* END IF CTRLC > 0 */
