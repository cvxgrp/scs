#include "util.h"

#include "glbopts.h"
#include "linsys.h"

/* return milli-seconds */
#if (defined NOTIMER)

void SCS(tic)(SCS(timer) * t) {}
scs_float SCS(tocq)(SCS(timer) * t) { return NAN; }

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

void SCS(tic)(SCS(timer) * t) { clock_gettime(CLOCK_MONOTONIC, &t->tic); }

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

scs_float SCS(toc)(SCS(timer) * t) {
  scs_float time = SCS(tocq)(t);
  scs_printf("time: %8.4f milli-seconds.\n", time);
  return time;
}

scs_float SCS(str_toc)(char *str, SCS(timer) * t) {
  scs_float time = SCS(tocq)(t);
  scs_printf("%s - time: %8.4f milli-seconds.\n", str, time);
  return time;
}

void SCS(print_cone_data)(const ScsCone *k) {
  scs_int i;
  scs_printf("num zeros = %i\n", (int)k->f);
  scs_printf("num LP = %i\n", (int)k->l);
  scs_printf("num SOCs = %i\n", (int)k->qsize);
  scs_printf("soc array:\n");
  for (i = 0; i < k->qsize; i++) {
    scs_printf("%i\n", (int)k->q[i]);
  }
  scs_printf("num SDCs = %i\n", (int)k->ssize);
  scs_printf("sdc array:\n");
  for (i = 0; i < k->ssize; i++) {
    scs_printf("%i\n", (int)k->s[i]);
  }
  scs_printf("num ep = %i\n", (int)k->ep);
  scs_printf("num ed = %i\n", (int)k->ed);
  scs_printf("num PCs = %i\n", (int)k->psize);
  scs_printf("pow array:\n");
  for (i = 0; i < k->psize; i++) {
    scs_printf("%4f\n", (double)k->p[i]);
  }
}

void SCS(print_work)(const ScsWork *w) {
  scs_int i, l = w->n + w->m;
  scs_printf("\n u_t is \n");
  for (i = 0; i < l; i++) {
    scs_printf("%f\n", w->u_t[i]);
  }
  scs_printf("\n u is \n");
  for (i = 0; i < l; i++) {
    scs_printf("%f\n", w->u[i]);
  }
  scs_printf("\n v is \n");
  for (i = 0; i < l; i++) {
    scs_printf("%f\n", w->v[i]);
  }
}

void SCS(print_data)(const ScsData *d) {
  scs_printf("m = %i\n", (int)d->m);
  scs_printf("n = %i\n", (int)d->n);

  scs_printf("max_iters = %i\n", (int)d->stgs->max_iters);
  scs_printf("verbose = %i\n", (int)d->stgs->verbose);
  scs_printf("normalize = %i\n", (int)d->stgs->normalize);
  scs_printf("warm_start = %i\n", (int)d->stgs->warm_start);
  scs_printf("acceleration_lookback = %i\n",
             (int)d->stgs->acceleration_lookback);
  scs_printf("eps = %4f\n", d->stgs->eps);
  scs_printf("alpha = %4f\n", d->stgs->alpha);
  scs_printf("rho_x = %4f\n", d->stgs->rho_x);
  scs_printf("cg_rate = %4f\n", d->stgs->cg_rate);
  scs_printf("scale = %4f\n", d->stgs->scale);
  scs_printf("write_data_filename = %s\n",
             (char *)d->stgs->write_data_filename);
}

void SCS(print_array)(const scs_float *arr, scs_int n, const char *name) {
  scs_int i, j, k = 0;
  scs_int num_on_one_line = 10;
  scs_printf("\n");
  for (i = 0; i < n / num_on_one_line; ++i) {
    for (j = 0; j < num_on_one_line; ++j) {
      scs_printf("%s[%li] = %4f, ", name, (long)k, arr[k]);
      k++;
    }
    scs_printf("\n");
  }
  for (j = k; j < n; ++j) {
    scs_printf("%s[%li] = %4f, ", name, (long)j, arr[j]);
  }
  scs_printf("\n");
}

void SCS(free_data)(ScsData *d, ScsCone *k) {
  if (d) {
    scs_free(d->b);
    scs_free(d->c);
    scs_free(d->stgs);
    if (d->A) {
      SCS(free_a_matrix)(d->A);
    }
    scs_free(d);
  }
  if (k) {
    scs_free(k->q);
    scs_free(k->s);
    scs_free(k->p);
    scs_free(k);
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

/* assumes d->stgs already allocated memory */
void SCS(set_default_settings)(ScsData *d) {
  /* These constants are defined in include/glbopts.h */
  d->stgs->max_iters = MAX_ITERS;
  d->stgs->eps = EPS;
  d->stgs->alpha = ALPHA;
  d->stgs->rho_x = RHO_X;
  d->stgs->scale = SCALE;
  d->stgs->cg_rate = CG_RATE;
  d->stgs->verbose = VERBOSE;
  d->stgs->normalize = NORMALIZE;
  d->stgs->warm_start = WARM_START;
  d->stgs->acceleration_lookback = ACCELERATION_LOOKBACK;
  d->stgs->write_data_filename = WRITE_DATA_FILENAME;
}
