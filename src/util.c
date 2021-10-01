#include "util.h"

#include "glbopts.h"
#include "linsys.h"
#include "scs_matrix.h"

/* return milli-seconds */
#if (defined NOTIMER)

void SCS(tic)(SCS(timer) * t) {
}
scs_float SCS(tocq)(SCS(timer) * t) {
  return NAN;
}

#elif (defined _WIN32 || _WIN64 || defined _WINDLL)

void SCS(tic)(SCS(timer) * t) {
  QueryPerformanceFrequency(&t->freq);
  QueryPerformanceCounter(&t->tic);
}

scs_float SCS(tocq)(SCS(timer) * t) {
  QueryPerformanceCounter(&t->toc);
  return (1e3 * (t->toc.QuadPart - t->tic.QuadPart) /
          (scs_float)t->freq.QuadPart);
}

#elif (defined __APPLE__)
#include <stdint.h>

void SCS(tic)(SCS(timer) * t) {
  /* read current clock cycles */
  t->tic = mach_absolute_time();
}

scs_float SCS(tocq)(SCS(timer) * t) {
  uint64_t duration; /* elapsed time in clock cycles*/

  t->toc = mach_absolute_time();
  duration = t->toc - t->tic;

  /*conversion from clock cycles to nanoseconds*/
  mach_timebase_info(&(t->tinfo));
  duration *= t->tinfo.numer;
  duration /= t->tinfo.denom;

  return (scs_float)duration / 1e6;
}

#else

void SCS(tic)(SCS(timer) * t) {
  clock_gettime(CLOCK_MONOTONIC, &t->tic);
}

scs_float SCS(tocq)(SCS(timer) * t) {
  struct timespec temp;

  clock_gettime(CLOCK_MONOTONIC, &t->toc);

  if ((t->toc.tv_nsec - t->tic.tv_nsec) < 0) {
    temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec - 1;
    temp.tv_nsec = 1e9 + t->toc.tv_nsec - t->tic.tv_nsec;
  } else {
    temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec;
    temp.tv_nsec = t->toc.tv_nsec - t->tic.tv_nsec;
  }
  return (scs_float)temp.tv_sec * 1e3 + (scs_float)temp.tv_nsec / 1e6;
}

#endif

void SCS(free_data)(ScsData *d, ScsCone *k, ScsSettings *stgs) {
  if (d) {
    scs_free(d->b);
    scs_free(d->c);
    if (d->A) {
      SCS(free_scs_matrix)(d->A);
    }
    if (d->P) {
      SCS(free_scs_matrix)(d->P);
    }
    scs_free(d);
  }
  if (k) {
    scs_free(k->bu);
    scs_free(k->bl);
    scs_free(k->q);
    scs_free(k->s);
    scs_free(k->p);
    scs_free(k);
  }
  if (stgs) {
    scs_free(stgs);
  }
}

void SCS(free_sol)(ScsSolution *sol) {
  if (sol) {
    scs_free(sol->x);
    scs_free(sol->y);
    scs_free(sol->s);
    scs_free(sol);
  }
}

/* assumes stgs already allocated memory */
void SCS(set_default_settings)(ScsSettings *stgs) {
  /* These constants are defined in include/glbopts.h */
  stgs->max_iters = MAX_ITERS;
  stgs->eps_abs = EPS_ABS;
  stgs->eps_rel = EPS_REL;
  stgs->eps_infeas = EPS_INFEAS;
  stgs->alpha = ALPHA;
  stgs->rho_x = RHO_X;
  stgs->scale = SCALE;
  stgs->verbose = VERBOSE;
  stgs->normalize = NORMALIZE;
  stgs->warm_start = WARM_START;
  stgs->acceleration_lookback = ACCELERATION_LOOKBACK;
  stgs->acceleration_interval = ACCELERATION_INTERVAL;
  stgs->adaptive_scale = ADAPTIVE_SCALE;
  stgs->write_data_filename = WRITE_DATA_FILENAME;
  stgs->log_csv_filename = LOG_CSV_FILENAME;
  stgs->time_limit_secs = TIME_LIMIT_SECS;
}
