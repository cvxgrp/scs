#include "cones.h"
#include "glbopts.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"

#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* With OpenMP the PSD blocks are projected in parallel (cones.c). Each block
 * goes through the same LAPACK calls with the same workspace sizes as in the
 * serial projection, so with single-threaded BLAS the result is the same bit
 * for bit. Project one vector with 1 and with 4 threads and compare.
 *
 * 40 real blocks (four above LAPACK's blocking crossover, ~32, of different
 * sizes) and 4 complex ones: with workspaces sized per thread instead of the
 * serial sizes, the large blocks would round differently, and this mix shows
 * it. The largest block holds less than half of the work, so the parallel
 * path is used. Other cones before and after the PSD blocks check the
 * offsets. */

#ifdef _OPENMP
#if defined(__GNUC__) && defined(__ELF__)
/* threaded BLAS would round differently in the serial projection; pin the
 * libraries that can be pinned at run time (weak references: ELF only) */
#define PSD_PAR_PIN_BLAS (1)
#ifdef __cplusplus
extern "C" {
#endif
void openblas_set_num_threads(int) __attribute__((weak));
int openblas_get_num_threads(void) __attribute__((weak));
void MKL_Set_Num_Threads(int) __attribute__((weak));
int MKL_Get_Max_Threads(void) __attribute__((weak));
#ifdef __cplusplus
}
#endif
#else
#define PSD_PAR_PIN_BLAS (0) /* e.g. Accelerate: VECLIB_MAXIMUM_THREADS=1 */
#endif

static scs_int psd_par_s[] = {1,  2,  3, 5,  8, 13, 21, 1, 2,  3,  48, 5, 8, 13,
                              21, 1,  2, 3,  5, 8,  64, 21, 1, 2,  3,  5, 8, 13,
                              21, 1, 40, 3,  5, 8,  13, 36, 1, 2,  3,  5};
static scs_int psd_par_cs[] = {2, 3, 6, 10};
static scs_int psd_par_q[] = {4, 6};
#endif

static const char *test_psd_parallel(void) {
#ifndef _OPENMP
  return mu_skipped; /* the projection is serial without OpenMP */
#else
  ScsCone k = {0};
  ScsConeWork *cs, *cp;
  scs_float *v, *vs, *vp, *r_y;
  scs_int i, m, nt = 0;
  int nt0 = omp_get_max_threads();
  const char *fail = SCS_NULL;
#if PSD_PAR_PIN_BLAS
  int ob0 = openblas_get_num_threads ? openblas_get_num_threads() : 0;
  int mkl0 = MKL_Get_Max_Threads ? MKL_Get_Max_Threads() : 0;
  if (openblas_set_num_threads)
    openblas_set_num_threads(1);
  if (MKL_Set_Num_Threads)
    MKL_Set_Num_Threads(1);
#endif

  k.z = 5;
  k.l = 10;
  k.qsize = sizeof(psd_par_q) / sizeof(scs_int);
  k.q = psd_par_q;
  k.ssize = sizeof(psd_par_s) / sizeof(scs_int);
  k.s = psd_par_s;
  k.cssize = sizeof(psd_par_cs) / sizeof(scs_int);
  k.cs = psd_par_cs;
  k.ep = 2;
  m = k.z + k.l + 3 * k.ep;
  for (i = 0; i < k.qsize; ++i)
    m += k.q[i];
  for (i = 0; i < k.ssize; ++i)
    m += (k.s[i] * (k.s[i] + 1)) / 2;
  for (i = 0; i < k.cssize; ++i)
    m += k.cs[i] * k.cs[i];

  v = (scs_float *)scs_calloc(m, sizeof(scs_float));
  vs = (scs_float *)scs_calloc(m, sizeof(scs_float));
  vp = (scs_float *)scs_calloc(m, sizeof(scs_float));
  r_y = (scs_float *)scs_calloc(m, sizeof(scs_float));
  ran_start(99);
  for (i = 0; i < m; ++i) {
    v[i] = rand_scs_float();
    r_y[i] = 1.5 + rand_scs_float(); /* non-uniform metric */
  }
  memcpy(vs, v, m * sizeof(scs_float));
  memcpy(vp, v, m * sizeof(scs_float));

  omp_set_num_threads(1); /* an OpenMP BLAS then runs serially too */
  cs = SCS(init_cone)(&k, m);
  if (cs)
    SCS(proj_dual_cone)(vs, cs, SCS_NULL, r_y);
  omp_set_num_threads(4);
  cp = SCS(init_cone)(&k, m);
  if (cp)
    SCS(proj_dual_cone)(vp, cp, SCS_NULL, r_y);

  if (!cs || !cp) {
    fail = "test_psd_parallel: cone setup failed";
  } else {
    nt = cs->psd_nthreads == 1 ? cp->psd_nthreads : -1;
    scs_printf("psd blocks projected by %li threads\n", (long)nt);
    if (nt < 2)
      fail = "test_psd_parallel: parallel PSD projection not set up";
    else if (memcmp(vs, vp, m * sizeof(scs_float)))
      fail = "test_psd_parallel: parallel projection differs from serial";
  }

  SCS(finish_cone)(cs);
  SCS(finish_cone)(cp);
  scs_free(v);
  scs_free(vs);
  scs_free(vp);
  scs_free(r_y);
  omp_set_num_threads(nt0);
#if PSD_PAR_PIN_BLAS
  if (openblas_set_num_threads)
    openblas_set_num_threads(ob0);
  if (MKL_Set_Num_Threads)
    MKL_Set_Num_Threads(mkl0);
#endif
  return fail;
#endif
}
