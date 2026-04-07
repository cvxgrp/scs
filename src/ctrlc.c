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
#include <stdbool.h>

extern bool utIsInterruptPending(void);
extern bool utSetInterruptEnabled(bool);

static int istate;
void scs_start_interrupt_listener(void) {
  istate = (int)utSetInterruptEnabled(true);
}

void scs_end_interrupt_listener(void) {
  utSetInterruptEnabled((bool)istate);
}

int scs_is_interrupted(void) {
  return (int)utIsInterruptPending();
}

#elif (defined _WIN32 || _WIN64 || defined _WINDLL)
#include <windows.h>

static volatile LONG int_detected;
static INIT_ONCE ctrlc_init_once = INIT_ONCE_STATIC_INIT;
static CRITICAL_SECTION ctrlc_cs;
static int listener_count = 0;

static BOOL CALLBACK init_ctrlc_cs(PINIT_ONCE once, PVOID param, PVOID *ctx) {
  (void)once;
  (void)param;
  (void)ctx;
  InitializeCriticalSection(&ctrlc_cs);
  return TRUE;
}

static BOOL WINAPI scs_handle_ctrlc(DWORD dwCtrlType) {
  if (dwCtrlType != CTRL_C_EVENT) {
    return FALSE;
  }
  InterlockedExchange(&int_detected, 1);
  return TRUE;
}

void scs_start_interrupt_listener(void) {
  InitOnceExecuteOnce(&ctrlc_init_once, init_ctrlc_cs, NULL, NULL);
  EnterCriticalSection(&ctrlc_cs);
  if (listener_count == 0) {
    InterlockedExchange(&int_detected, 0);
    SetConsoleCtrlHandler(scs_handle_ctrlc, TRUE);
  }
  listener_count++;
  LeaveCriticalSection(&ctrlc_cs);
}

void scs_end_interrupt_listener(void) {
  InitOnceExecuteOnce(&ctrlc_init_once, init_ctrlc_cs, NULL, NULL);
  EnterCriticalSection(&ctrlc_cs);
  if (listener_count > 0) {
    listener_count--;
    if (listener_count == 0) {
      SetConsoleCtrlHandler(scs_handle_ctrlc, FALSE);
    }
  }
  LeaveCriticalSection(&ctrlc_cs);
}

int scs_is_interrupted(void) {
  return (int)InterlockedCompareExchange(&int_detected, 0, 0);
}

#else /* Unix */

#include <signal.h>
#include <pthread.h>

static volatile sig_atomic_t int_detected;
static struct sigaction oact;
static pthread_mutex_t ctrlc_mutex = PTHREAD_MUTEX_INITIALIZER;
static int listener_count = 0;

static void scs_handle_ctrlc(int dummy) {
  int_detected = dummy ? dummy : -1;
}

void scs_start_interrupt_listener(void) {
  pthread_mutex_lock(&ctrlc_mutex);
  if (listener_count == 0) {
    struct sigaction act;
    int_detected = 0;
    act.sa_flags = 0;
    sigemptyset(&act.sa_mask);
    act.sa_handler = scs_handle_ctrlc;
    sigaction(SIGINT, &act, &oact);
  }
  listener_count++;
  pthread_mutex_unlock(&ctrlc_mutex);
}

void scs_end_interrupt_listener(void) {
  pthread_mutex_lock(&ctrlc_mutex);
  if (listener_count > 0) {
    listener_count--;
    if (listener_count == 0) {
      struct sigaction act;
      sigaction(SIGINT, &oact, &act);
    }
  }
  pthread_mutex_unlock(&ctrlc_mutex);
}

int scs_is_interrupted(void) {
  return int_detected;
}

#endif /* END IF MATLAB_MEX_FILE / WIN32 */

#endif /* END IF CTRLC > 0 */
