#include "scs.h"

#include <dlfcn.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Checks, with a real MKL, that the parallel PSD projection keeps MKL at one
 * thread in each of its threads: the calls of MKL_Set_Num_Threads_Local around
 * the projection in src/cones.c. Built and run only by
 * .github/workflows/build_openmp_mkl.yml.
 *
 * That guard matters when MKL cannot see the parallel region of SCS, because
 * the two run on different OpenMP runtimes: a GNU-compiled SCS on libgomp and
 * MKL on libiomp5, its default threading layer. Wrappers get them that way (a
 * Python extension module, a Julia JLL), each library with its own symbols,
 * and this program loads SCS as they do: it links neither runtime, and opens
 * the SCS library given as its argument with a local scope. All symbols are
 * bound at that point (RTLD_NOW), before MKL loads libiomp5, so the regions of
 * SCS stay on libgomp.
 *
 * It solves an SDP with four PSD blocks, large enough for MKL to use threads,
 * and prints the number of threads of the process after the solve. The
 * workflow compares that number between three runs:
 *   - SCS_PSD_NUM_THREADS=1: the blocks are projected in order,
 *   - the parallel projection,
 *   - the parallel projection of a library built without the guard.
 * With the guard the parallel projection starts no MKL thread, so it ends with
 * the threads of the serial run. Without it every thread of SCS gets MKL
 * threads of its own, which shows that the comparison can tell.
 *
 * usage: mkl_psd_guard <SCS shared library> [iterations]
 * With a number of iterations the run stops there and only counts threads
 * (they all start in the first iterations); without, the solve must succeed.
 *
 * The problem is: minimize x subject to I - x * M_k >= 0 (PSD) for each block
 * k, with M_k = I - (1 + a_k) v_k v_k' for a unit vector v_k. M_k has the
 * eigenvalues 1 and -a_k, so the constraints are -1 / a_k <= x <= 1. The
 * largest a_k is 1, so the solution is x = -1.
 */

#define PSD_N 200
#define PSD_BLOCKS 4
#define PSD_LEN (PSD_N * (PSD_N + 1) / 2)

typedef ScsWork *(*init_fn)(const ScsData *, const ScsCone *,
                            const ScsSettings *);
typedef scs_int (*solve_fn)(ScsWork *, ScsSolution *, ScsInfo *, scs_int);
typedef void (*finish_fn)(ScsWork *);
typedef void (*defaults_fn)(ScsSettings *);
typedef int (*mkl_set_local_fn)(int);
typedef int (*mkl_get_max_fn)(void);

/* uniform in (-1, 1), from a linear congruential generator */
static double next_rand(unsigned long *state) {
  *state = (*state * 1103515245UL + 12345UL) & 0x7fffffffUL;
  return 2. * ((double)*state + 0.5) / 2147483648. - 1.;
}

/* Number of threads of this process, or -1 where /proc does not give it */
static int os_threads(void) {
  char line[256];
  int n = -1;
  FILE *f = fopen("/proc/self/status", "r");
  if (!f) {
    return -1;
  }
  while (fgets(line, sizeof(line), f)) {
    if (sscanf(line, "Threads: %d", &n) == 1) {
      break;
    }
  }
  fclose(f);
  return n;
}

int main(int argc, char **argv) {
  scs_int m = PSD_BLOCKS * PSD_LEN, blocks[PSD_BLOCKS];
  scs_int i, j, k, idx = 0, status = -1;
  scs_float c[1] = {1.};
  scs_int A_p[2] = {0, 0};
  scs_float *A_x, *b;
  scs_int *A_i;
  double v[PSD_N], nrm, a, x = 0.;
  unsigned long seed = 1234;
  int mkl_nt = 0, mkl_now = 0, threads, fail = 0;
  int short_run = argc > 2;
  ScsMatrix A = {0};
  ScsData d = {0};
  ScsCone cone = {0};
  ScsSettings stgs;
  ScsSolution sol = {0};
  ScsInfo info;
  ScsWork *w;
  void *lib;
  init_fn init = 0;
  solve_fn solve = 0;
  finish_fn finish = 0;
  defaults_fn defaults = 0;
  mkl_set_local_fn mkl_set_local = 0;
  mkl_get_max_fn mkl_get_max = 0;

  if (argc < 2) {
    fprintf(stderr, "usage: %s <SCS shared library> [iterations]\n", argv[0]);
    return 1;
  }
  setvbuf(stdout, NULL, _IOLBF, 0); /* keep the output if the solve crashes */
  lib = dlopen(argv[1], RTLD_NOW | RTLD_LOCAL);
  if (!lib) {
    fprintf(stderr, "cannot load %s: %s\n", argv[1], dlerror());
    return 1;
  }
  /* the assignment through void ** is the conversion POSIX asks for */
  *(void **)(&init) = dlsym(lib, "scs_init");
  *(void **)(&solve) = dlsym(lib, "scs_solve");
  *(void **)(&finish) = dlsym(lib, "scs_finish");
  *(void **)(&defaults) = dlsym(lib, "scs_set_default_settings");
  if (!init || !solve || !finish || !defaults) {
    fprintf(stderr, "%s is not an SCS library\n", argv[1]);
    return 1;
  }
  /* MKL's own functions, found among the libraries that SCS brings in */
  *(void **)(&mkl_set_local) = dlsym(lib, "MKL_Set_Num_Threads_Local");
  *(void **)(&mkl_get_max) = dlsym(lib, "MKL_Get_Max_Threads");

  A_x = (scs_float *)calloc(m, sizeof(scs_float));
  A_i = (scs_int *)calloc(m, sizeof(scs_int));
  b = (scs_float *)calloc(m, sizeof(scs_float));
  if (!A_x || !A_i || !b) {
    fprintf(stderr, "out of memory\n");
    return 1;
  }
  /* column 0 of A is M_k and b is the identity, block after block: the lower
   * triangle by columns, off-diagonal entries scaled by sqrt(2) */
  for (k = 0; k < PSD_BLOCKS; ++k) {
    blocks[k] = PSD_N;
    a = 1. - 0.1 * (double)k;
    nrm = 0.;
    for (i = 0; i < PSD_N; ++i) {
      v[i] = next_rand(&seed);
      nrm += v[i] * v[i];
    }
    nrm = sqrt(nrm);
    for (i = 0; i < PSD_N; ++i) {
      v[i] /= nrm;
    }
    for (j = 0; j < PSD_N; ++j) {
      for (i = j; i < PSD_N; ++i) {
        A_i[idx] = idx;
        if (i == j) {
          A_x[idx] = 1. - (1. + a) * v[i] * v[i];
          b[idx] = 1.;
        } else {
          A_x[idx] = -sqrt(2.) * (1. + a) * v[i] * v[j];
        }
        ++idx;
      }
    }
  }
  A_p[1] = m;
  A.x = A_x;
  A.i = A_i;
  A.p = A_p;
  A.m = m;
  A.n = 1;
  d.m = m;
  d.n = 1;
  d.A = &A;
  d.b = b;
  d.c = c;
  cone.s = blocks;
  cone.ssize = PSD_BLOCKS;
  defaults(&stgs);
  stgs.verbose = 1; /* the header says whether the projection is parallel */
  if (short_run) {
    stgs.max_iters = (scs_int)atoi(argv[2]);
  }

  w = init(&d, &cone, &stgs);
  if (w) {
    if (mkl_set_local && mkl_get_max) {
      /* A thread-local setting for this thread, equal to the global one so
       * that it changes nothing: the projection must leave it as it is. */
      mkl_nt = mkl_get_max();
      mkl_set_local(mkl_nt);
    }
    status = solve(w, &sol, &info, 0);
    x = sol.x ? (double)sol.x[0] : 0.;
  }
  threads = os_threads();
  if (mkl_set_local && mkl_get_max) {
    mkl_now = mkl_set_local(0); /* returns the setting and removes it */
  }

  printf("mkl_psd_guard: status=%d iters=%ld x=%.6f os_threads=%d "
         "mkl_threads=%d mkl_threads_after=%d\n",
         (int)status, w ? (long)info.iter : 0L, x, threads, mkl_nt, mkl_now);
  if (!w) {
    fprintf(stderr, "mkl_psd_guard: scs_init failed\n");
    fail = 1;
  } else if (!short_run && (status != SCS_SOLVED || fabs(x + 1.) > 1e-2)) {
    fprintf(stderr, "mkl_psd_guard: wrong solve, expected status %d, x = -1\n",
            (int)SCS_SOLVED);
    fail = 1;
  }
  if (mkl_now != mkl_nt) {
    fprintf(stderr,
            "mkl_psd_guard: the solve changed the MKL thread setting of its "
            "caller from %d to %d\n",
            mkl_nt, mkl_now);
    fail = 1;
  }
  if (w) {
    finish(w);
  }
  free(sol.x);
  free(sol.y);
  free(sol.s);
  free(A_x);
  free(A_i);
  free(b);
  /* the library stays loaded: its OpenMP threads are still there */
  return fail;
}
