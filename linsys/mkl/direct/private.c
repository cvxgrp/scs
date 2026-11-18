#include "private.h"
#include <string.h>

#define PARDISO_SYMBOLIC (11)
#define PARDISO_NUMERIC (22)
#define PARDISO_SOLVE (33)
#define PARDISO_CLEANUP (-1)

/*
 * MKL interface layer constants. MKL has two integer interfaces:
 *
 *   LP64  (MKL_INTERFACE_LP64  = 0): BLAS/LAPACK use 32-bit integers (int).
 *   ILP64 (MKL_INTERFACE_ILP64 = 1): BLAS/LAPACK use 64-bit integers (long long).
 *
 * These affect the standard BLAS/LAPACK symbols (dgemm, dpotrf, etc.).
 * PARDISO has separate entry points for each integer width:
 *
 *   pardiso    — 32-bit integer indices (used when !DLONG)
 *   pardiso_64 — 64-bit integer indices (used when DLONG)
 *
 * The pardiso/pardiso_64 choice is independent of the interface layer; each is
 * a distinct symbol that always uses its own integer width regardless of what
 * MKL_Set_Interface_Layer says.
 *
 * The BLAS integer width is controlled by the use_blas64 meson option:
 *
 *   use_blas64=false (default) -> links mkl-dynamic-lp64-*,  expects LP64
 *   use_blas64=true            -> links mkl-dynamic-ilp64-*, expects ILP64
 *
 * See meson.build for the linkage logic. The MKL-specific initialization in
 * src/scs.c uses BLAS64 to pick the right expected interface layer.
 *
 * PARDISO is independent: pardiso_64 always uses 64-bit ints regardless of
 * the interface layer or BLAS64 setting.
 */
#define MKL_INTERFACE_LP64 0
#define MKL_INTERFACE_ILP64 1

#ifdef DLONG
#define _PARDISO pardiso_64
#else
#define _PARDISO pardiso
#endif

/* Prototypes for Pardiso and MKL service functions. */
void _PARDISO(void **pt, const scs_int *maxfct, const scs_int *mnum,
              const scs_int *mtype, const scs_int *phase, const scs_int *n,
              const scs_float *a, const scs_int *ia, const scs_int *ja,
              scs_int *perm, const scs_int *nrhs, scs_int *iparm,
              const scs_int *msglvl, scs_float *b, scs_float *x,
              scs_int *error);
/* BLAS64 requires DLONG for MKL builds: the interface layer (ILP64) set by
 * MKL_Set_Interface_Layer must match the PARDISO entry point (pardiso_64).
 * Without DLONG, pardiso (32-bit) is used but ILP64 makes its internal BLAS
 * calls expect 64-bit integers, causing hangs or memory corruption. */
#if defined(BLAS64) && !defined(DLONG)
#error "MKL PARDISO requires DLONG when BLAS64 is set (pardiso_64 needs 64-bit ints)"
#endif

/* Small linear systems: Pardiso on one thread.
 *
 * What it does. With a threaded MKL, every Pardiso call (the symbolic phase,
 * each numerical factorization, the solve of every iteration) runs on MKL's
 * threads (MKL_NUM_THREADS, else OMP_NUM_THREADS). For a small system this
 * workspace instead sets MKL's thread count to one for the calling thread
 * around each of its Pardiso calls (MKL_Set_Num_Threads_Local) and restores
 * the previous setting after the call. "Small" is decided once, at setup:
 *   - when the KKT matrix has fewer than PARDISO_MIN_THREADED_NNZ nonzeros
 *     (upper triangle, as Pardiso gets it), the symbolic phase runs on one
 *     thread;
 *   - if the factor it computes (iparm[17]) has fewer than
 *     PARDISO_MIN_THREADED_NNZ nonzeros too, every later call (the
 *     factorizations, the solves) runs on one thread as well. If not, the
 *     handle is released and the symbolic phase is redone on MKL's threads,
 *     and everything runs on them.
 * One handle never mixes thread counts between its phases: Pardiso keeps the
 * thread count of the symbolic phase in what it sets up there, and a
 * symbolic phase on one thread followed by solves on 4 threads made each
 * solve 6.6x slower (CONT-050). The factor is never smaller than the KKT
 * matrix, so a system that is large by the first test is never redone.
 * Pardiso's control parameters are the same either way, and MKL's BLAS and
 * LAPACK (the PSD projection, Anderson acceleration) keep their threads.
 *
 * Why the factor nonzeros: the solve of every iteration is a forward and a
 * backward substitution through the factor, so its work is proportional to
 * them, while the KKT nonzeros can be 25 times fewer (CVXQP1_L: 8.0e4 in the
 * KKT matrix, 2.1e6 in the factor, 1.9x faster per iteration on 4 threads).
 *
 * Measured (November 2026, oneMKL 2026.1, static gnu_thread layer, 4
 * threads, AMD Zen 3 on Modal, 34 LP, QP and SDP problems, 3 repeats):
 *   - one Pardiso solve costs about 25-30 microseconds more on 4 threads
 *     than on one, so tiny problems iterated about 2x slower on 4 threads;
 *   - the break-even point, all of Pardiso on one thread against all on 4,
 *     lies between about 3e4 and 6e4 factor nonzeros and depends on the
 *     problem (PRIMAL3 3.4e4: one thread 1.12x slower; QSHIP08S 4.4e4: 0.81x;
 *     QSHIP12S 5.4e4: 1.32x), and above it threads win clearly (DTOC3 1.1e5:
 *     one thread 1.8x slower, CONT-201 5.0e6: 2.3x);
 *   - with the threshold at 5e4, ms per iteration with the rule against
 *     without it: geometric mean 0.845 below the threshold (15 problems) with
 *     MKL_DYNAMIC=FALSE and 0.968 with MKL's default MKL_DYNAMIC=TRUE (which
 *     already runs the tiniest systems on one thread by itself), and 1.005
 *     and 1.025 above it (19 problems, where the rule changes nothing and the
 *     ratio is run-to-run noise). No problem was slower by more than that
 *     noise. A threshold of 1e5 made problems between 6e4 and 1e5 6-30%
 *     slower.
 *   - When the symbolic phase is redone, setup takes 25-48 ms longer
 *     (CONT-050, biggs, QSHIP12S, AUG3DQP), once per solve.
 *
 * Tradeoff. One fixed threshold for every machine: a system just below it
 * might have been a little faster on threads, one just above it a little
 * slower, and the break-even point moves with the CPU, the thread count and
 * the MKL version.
 *
 * Deterministic: the choice depends only on the problem data, once at setup,
 * never on timings or load. The factor size it reads always comes from a
 * symbolic phase run on one thread (the factor is only tested when the KKT
 * matrix is small), so the decision is the same whatever the thread count.
 *
 * Where it does nothing: a sequential MKL (the Python wheels) and MKL on one
 * thread. It does not remove the per-call overhead of a threaded MKL on one
 * thread, which is higher than that of the sequential MKL, so small problems
 * stay faster with a sequential MKL.
 *
 * MKL_Set_Num_Threads_Local exists in every MKL threading layer (in the
 * sequential one it does nothing) and in mkl_rt, so the reference is strong.
 * A weak reference, as in src/cones.c, does not pull the function out of a
 * static MKL archive and stays NULL there. This strong reference also
 * resolves the weak one in src/cones.c, so the MKL guard of the parallel PSD
 * projection is active in static MKL builds too. */
#ifndef PARDISO_MIN_THREADED_NNZ
#define PARDISO_MIN_THREADED_NNZ (50000)
#endif

int MKL_Set_Num_Threads_Local(int nt);

/* One Pardiso call of workspace p, on one MKL thread when p->one_thread */
static void pardiso_call(ScsLinSysWork *p, scs_int phase, const scs_float *a,
                         scs_float *b) {
  int prev = 0;
  if (p->one_thread) {
    prev = MKL_Set_Num_Threads_Local(1);
  }
  _PARDISO(p->pt, &(p->maxfct), &(p->mnum), &(p->mtype), &phase,
           &(p->n_plus_m), a, p->kkt->p, p->kkt->i, SCS_NULL, &(p->nrhs),
           p->iparm, &(p->msglvl), b, SCS_NULL, &(p->error));
  if (p->one_thread) {
    MKL_Set_Num_Threads_Local(prev);
  }
}

static scs_int prod_setup(ScsLinSysWork *p);
static scs_int prod_make_handle(ScsLinSysWork *p);
static void prod_free_handle(ScsLinSysWork *p);
static void prod_free(ScsLinSysWork *p);

const char *scs_get_lin_sys_method() {
  return "sparse-direct-mkl-pardiso";
}

void scs_free_lin_sys_work(ScsLinSysWork *p) {
  if (p) {
    /* Pardiso cleanup needs the KKT index arrays; on early allocation
     * failures the workspace may be freed before the KKT matrix (or
     * the factorization) ever existed, so gate the solver teardown on
     * what was actually initialized. The sparse BLAS handle of the
     * products refers to the KKT arrays, so it goes first. */
    prod_free(p);
    if (p->kkt) {
      /* release on the handle's own thread count, like every other call */
      p->phase = PARDISO_CLEANUP;
      pardiso_call(p, PARDISO_CLEANUP, SCS_NULL, SCS_NULL);
      if (p->error != 0) {
        scs_printf("Error during MKL Pardiso cleanup: %d", (int)p->error);
      }
      SCS(cs_spfree)(p->kkt);
    }
    if (p->diag_r_idxs)
      scs_free(p->diag_r_idxs);
    if (p->diag_p)
      scs_free(p->diag_p);
    scs_free(p);
  }
}

ScsLinSysWork *scs_init_lin_sys_work(const ScsMatrix *A, const ScsMatrix *P,
                                     const scs_float *diag_r) {
  scs_int i;
  ScsLinSysWork *p = scs_calloc(1, sizeof(ScsLinSysWork));
  if (!p)
    return SCS_NULL;
  p->n = A->n;
  p->m = A->m;
  p->n_plus_m = p->n + p->m;

  p->diag_r_idxs = (scs_int *)scs_calloc(p->n_plus_m, sizeof(scs_int));
  p->diag_p = (scs_float *)scs_calloc(p->n, sizeof(scs_float));
  if (!p->diag_r_idxs || !p->diag_p) {
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  /* MKL pardiso requires upper triangular CSR matrices. The KKT matrix stuffed
   * as CSC lower triangular is equivalent. Pass upper=0. */
  p->kkt = SCS(form_kkt)(A, P, p->diag_p, diag_r, p->diag_r_idxs, 0);
  if (!(p->kkt)) {
    scs_printf("Error in forming KKT matrix");
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }
  /* small systems: Pardiso on one thread, see PARDISO_MIN_THREADED_NNZ */
  p->one_thread = p->kkt->p[p->n_plus_m] < PARDISO_MIN_THREADED_NNZ;

  for (i = 0; i < 64; i++) {
    p->iparm[i] = 0; /* Setup Pardiso control parameters */
    p->pt[i] = 0;    /* Initialize the internal solver memory pointer */
  }

  /* Set Pardiso variables */
  p->mtype = -2;         /* Real symmetric indefinite matrix */
  p->nrhs = 1;           /* Number of right hand sides */
  p->maxfct = 1;         /* Maximum number of numerical factorizations */
  p->mnum = 1;           /* Which factorization to use */
  p->error = 0;          /* Initialize error flag */
  p->msglvl = VERBOSITY; /* Printing information */

  /* For all iparm vars see the oneMKL "pardiso iparm Parameter" table, which
   * uses C, zero-based indices (iparm[24] here is iparm(25) in the Fortran
   * docs; iparm[34] below only sets the indexing of the matrix arrays, not of
   * iparm). With iparm[0] = 1 every entry not set here is 0 rather than
   * Intel's default. That differs from the default for iparm[20]: 0 is 1x1
   * diagonal pivoting (the default would be Bunch-Kaufman), enough for SCS's
   * quasi-definite KKT matrix. iparm[10] = iparm[12] = 0 (no scaling or
   * matching) is the default for symmetric indefinite matrices and what the
   * two-level factorization below requires. */
  p->iparm[0] = 1;  /* Pardiso must inspect iparm */
  p->iparm[1] = 3;  /* Parallel (OpenMP) nested dissection reordering */
  p->iparm[5] = 1;  /* Write solution into b */
  p->iparm[7] = 0;  /* Automatic iterative refinement calculation */
  p->iparm[9] = 13; /* Perturb small pivots with 1e-13 (Intel's default for
                       symmetric indefinite matrices is 8, i.e. 1e-8) */
  /* Two-level scheduling for the parallel factorization. Measured on oneMKL
   * 2026.1 (2026-11): 0 (the classic algorithm) segfaulted in the solve phase
   * on every run, and 10 (the improved two-level algorithm) broke the inertia
   * check below (iparm[21]), so keep 1. */
  p->iparm[23] = 1;
  /* Forward/backward solve, run in every iteration. Intel's table: 0 (the
   * default) is the parallel solve, which for one right-hand side, as here,
   * partitions the matrix; 1 is the sequential solve; 2 always partitions the
   * matrix. SCS set 1 up to and including 3.3.2, with a comment calling it
   * parallel. Measured on oneMKL 2026.1 (2026-11, threaded static MKL,
   * AMD Zen 3 (Milan) on Modal, 4 threads): 0, 1 and 2 all gave a parallel
   * solve, 2-2.7x faster than on one thread, with no difference between them
   * beyond run-to-run noise; the value was traced to reach Pardiso as
   * element 24, so this was not an indexing mix-up. Why 1 is not sequential
   * there is unknown; two untested guesses: only the traversal of the
   * elimination tree is sequential while the dense kernels inside it still
   * use MKL's threads, or the two-level factorization (iparm[23] = 1) has a
   * solve path of its own that ignores iparm[24]. 0 is set because it is the
   * documented parallel solve, so SCS does not depend on that undocumented
   * behaviour if a later MKL honours 1 = sequential. With a sequential MKL
   * (as in the PyPI wheel), or for a small system (PARDISO_MIN_THREADED_NNZ
   * above), every value runs on one thread. */
  p->iparm[24] = 0;
  p->iparm[17] = -1;        /* Report the factor nonzeros after phase 11 */
  p->iparm[34] = 1;         /* Use C-style indexing for indices */
  /* p->iparm[36] = -80; */ /* Form block sparse matrices */

#ifdef SFLOAT
  p->iparm[27] = 1; /* 1 is single precision, 0 is double */
#endif

  /* Permutation and symbolic factorization */
  pardiso_call(p, PARDISO_SYMBOLIC, p->kkt->x, SCS_NULL);

  if (p->error != 0) {
    scs_printf("Error during symbolic factorization: %d", (int)p->error);
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  /* small systems, see PARDISO_MIN_THREADED_NNZ: the symbolic phase ran on
   * one thread because the KKT matrix is small, but the factor is not. One
   * handle must keep one thread count in all its phases, so release it and
   * redo the symbolic phase on MKL's threads. */
  if (p->one_thread && p->iparm[17] >= PARDISO_MIN_THREADED_NNZ) {
    pardiso_call(p, PARDISO_CLEANUP, SCS_NULL, SCS_NULL);
    for (i = 0; i < 64; i++) {
      p->pt[i] = 0;
    }
    p->one_thread = 0;
    p->iparm[17] = -1; /* report the factor nonzeros again */
    pardiso_call(p, PARDISO_SYMBOLIC, p->kkt->x, SCS_NULL);
    if (p->error != 0) {
      scs_printf("Error during symbolic factorization: %d", (int)p->error);
      scs_free_lin_sys_work(p);
      return SCS_NULL;
    }
  }

  /* Numerical factorization */
  p->phase = PARDISO_NUMERIC;
  pardiso_call(p, PARDISO_NUMERIC, p->kkt->x, SCS_NULL);

  if (p->error) {
    scs_printf("Error during numerical factorization: %d", (int)p->error);
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  if (p->iparm[21] < p->n) {
    scs_printf("KKT matrix has < n positive eigenvalues. P not PSD.\n");
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  /* the refinement's products with K in MKL's sparse BLAS, except for a K too
   * large for MKL's 32-bit sparse integers (prod_setup's size rule, then SCS
   * computes them): if they cannot be set up, setup fails (prod_setup says
   * why) rather than leaving them to SCS */
  if (prod_setup(p) != 0) {
    scs_printf("MKL Pardiso: the refinement's products with K could not be "
               "set up\n");
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }
  p->prod_ok = !p->prod_none;

  return p;
}

/* Returns solution to linear system Ax = b with solution stored in b */
scs_int scs_solve_lin_sys(ScsLinSysWork *p, scs_float *b, const scs_float *ws,
                          scs_float tol) {
  /* Back substitution and iterative refinement */
  p->phase = PARDISO_SOLVE;
  pardiso_call(p, PARDISO_SOLVE, p->kkt->x, b);
  if (p->error != 0) {
    scs_printf("Error during linear system solution: %d", (int)p->error);
  }
  return p->error;
}

/* Update factorization when R changes */
scs_int scs_update_lin_sys_diag_r(ScsLinSysWork *p, const scs_float *diag_r) {
  scs_int i;

  for (i = 0; i < p->n; ++i) {
    /* top left is R_x + P, bottom right is -R_y */
    p->kkt->x[p->diag_r_idxs[i]] = p->diag_p[i] + diag_r[i];
  }
  for (i = p->n; i < p->n + p->m; ++i) {
    /* top left is R_x + P, bottom right is -R_y */
    p->kkt->x[p->diag_r_idxs[i]] = -diag_r[i];
  }

  /* Perform numerical factorization */
  p->phase = PARDISO_NUMERIC;
  pardiso_call(p, PARDISO_NUMERIC, p->kkt->x, SCS_NULL);

  if (p->error != 0) {
    scs_printf("Error in PARDISO factorization when updating: %d.\n",
               (int)p->error);
    return (scs_int)p->error;
  }

  /* MKL may keep its own copy of the values after mkl_sparse_optimize, so
   * the product's handle is rebuilt from the updated kkt->x */
  if (p->prod_ok) {
    prod_free_handle(p);
    if (prod_make_handle(p) != 0) {
      scs_printf("Error in MKL sparse BLAS when updating.\n");
      return -1;
    }
  }
  return 0;
}

/* ------------------------------------------------------------------------
 * Products with K for SCS's per-iteration refinement (optional, linsys.h),
 * in MKL's inspector-executor sparse BLAS: K is the upper triangle Pardiso
 * gets (kkt, CSR), passed as a symmetric matrix, analysed once by
 * mkl_sparse_optimize, and r = b - K x is one mkl_sparse_?_mv. The guess x0
 * and r are host vectors of this workspace.
 *
 * Threads. The product runs on MKL's threads, except for a small system,
 * where Pardiso runs on one thread (one_thread, see PARDISO_MIN_THREADED_NNZ
 * above): there the product runs on one thread too, around every MKL sparse
 * call on this handle (create, optimize, mv), so that the handle never
 * mixes thread counts either. Measured per product against SCS's own
 * single-threaded one (November 2026, oneMKL 2026.1, AMD Zen 3 on Modal, 4
 * threads, 46 LPs and QPs): 0.2 to 0.5 times its time on MKL's threads
 * (co-100: 0.76 against 2.39 ms) and 0.45 to 0.93 times on one thread, so
 * MKL's product is the faster one at every size; forcing 4 threads on the
 * smallest systems made the product about 2 times slower than on one, which
 * is why it follows Pardiso's thread rule.
 *
 * Reproducibility: on a fixed number of threads the product gave the same
 * bits in repeated runs (same measurements); a different thread count can
 * change its last bits.
 *
 * Setup cost: the handle is created and analysed (mkl_sparse_optimize) at
 * setup and again at every scs_update_lin_sys_diag_r (each scale update),
 * whether or not the refinement turns on, since the every-25-iterations
 * check uses the product too. That cost has not been measured on its own.
 *
 * The functions are declared here, like Pardiso's, so that no MKL header is
 * needed; the types and enum values are those of mkl_spblas.h. MKL's integer
 * is 64-bit with BLAS64 (ILP64, as for Pardiso) and in MATLAB, else 32-bit;
 * when that differs from scs_int (DLONG without BLAS64) the indices are
 * copied once; a K too large for 32-bit indices gets SCS's own products
 * (not provided, by size), while a BLAS64=1 build provides them at any
 * size. scs_mkl_int has the width of scs_blas.h's blas_int, which is defined
 * only with USE_LAPACK; an MKL build need not use LAPACK, so it is defined
 * here.
 * ---------------------------------------------------------------------- */
#if defined(MATLAB_MEX_FILE)
#include <stddef.h>
typedef ptrdiff_t scs_mkl_int;
#elif defined(BLAS64)
typedef long long scs_mkl_int;
#else
typedef int scs_mkl_int;
#endif

typedef struct scs_mkl_sparse_matrix *scs_mkl_sparse_t; /* sparse_matrix_t */
struct scs_mkl_matrix_descr { /* struct matrix_descr */
  int type;                   /* sparse_matrix_type_t */
  int mode;                   /* sparse_fill_mode_t */
  int diag;                   /* sparse_diag_type_t */
};
#define SCS_MKL_SPARSE_STATUS_SUCCESS (0)
#define SCS_MKL_SPARSE_OPERATION_NON_TRANSPOSE (10)
#define SCS_MKL_SPARSE_MATRIX_TYPE_SYMMETRIC (21)
#define SCS_MKL_SPARSE_INDEX_BASE_ZERO (0)
#define SCS_MKL_SPARSE_FILL_MODE_UPPER (41)
#define SCS_MKL_SPARSE_DIAG_NON_UNIT (50)

#ifndef SFLOAT
#define SCS_MKL_SPARSE_CREATE_CSR mkl_sparse_d_create_csr
#define SCS_MKL_SPARSE_MV mkl_sparse_d_mv
#else
#define SCS_MKL_SPARSE_CREATE_CSR mkl_sparse_s_create_csr
#define SCS_MKL_SPARSE_MV mkl_sparse_s_mv
#endif
int SCS_MKL_SPARSE_CREATE_CSR(scs_mkl_sparse_t *A, int indexing,
                              scs_mkl_int rows, scs_mkl_int cols,
                              scs_mkl_int *rows_start, scs_mkl_int *rows_end,
                              scs_mkl_int *col_indx, scs_float *values);
int SCS_MKL_SPARSE_MV(int operation, scs_float alpha,
                      const scs_mkl_sparse_t A,
                      struct scs_mkl_matrix_descr descr, const scs_float *x,
                      scs_float beta, scs_float *y);
int mkl_sparse_set_mv_hint(const scs_mkl_sparse_t A, int operation,
                           struct scs_mkl_matrix_descr descr,
                           scs_mkl_int expected_calls);
int mkl_sparse_optimize(scs_mkl_sparse_t A);
int mkl_sparse_destroy(scs_mkl_sparse_t A);

static struct scs_mkl_matrix_descr prod_descr(void) {
  struct scs_mkl_matrix_descr d;
  d.type = SCS_MKL_SPARSE_MATRIX_TYPE_SYMMETRIC;
  d.mode = SCS_MKL_SPARSE_FILL_MODE_UPPER;
  d.diag = SCS_MKL_SPARSE_DIAG_NON_UNIT;
  return d;
}

/* the MKL thread count of every sparse BLAS call on this handle: one for a
 * small system, as for Pardiso (see the threads note above) */
static int prod_threads_begin(const ScsLinSysWork *p) {
  return p->one_thread ? MKL_Set_Num_Threads_Local(1) : 0;
}

static void prod_threads_end(const ScsLinSysWork *p, int prev) {
  if (p->one_thread) {
    MKL_Set_Num_Threads_Local(prev);
  }
}

static void prod_free_handle(ScsLinSysWork *p) {
  if (p->prod_k) {
    int prev = prod_threads_begin(p);
    mkl_sparse_destroy((scs_mkl_sparse_t)p->prod_k);
    prod_threads_end(p, prev);
    p->prod_k = SCS_NULL;
  }
}

/* the handle of K (upper triangle of kkt as a symmetric matrix), analysed for
 * many products; uses kkt->x as it is now */
static scs_int prod_make_handle(ScsLinSysWork *p) {
  scs_mkl_sparse_t A = SCS_NULL;
  scs_mkl_int nm = (scs_mkl_int)p->n_plus_m;
  scs_mkl_int *ptr = p->prod_ptr ? (scs_mkl_int *)p->prod_ptr
                                 : (scs_mkl_int *)p->kkt->p;
  scs_mkl_int *ind = p->prod_ind ? (scs_mkl_int *)p->prod_ind
                                 : (scs_mkl_int *)p->kkt->i;
  struct scs_mkl_matrix_descr descr = prod_descr();
  int prev = prod_threads_begin(p), st;
  st = SCS_MKL_SPARSE_CREATE_CSR(&A, SCS_MKL_SPARSE_INDEX_BASE_ZERO, nm, nm,
                                 ptr, ptr + 1, ind, p->kkt->x);
  if (st == SCS_MKL_SPARSE_STATUS_SUCCESS) {
    p->prod_k = (void *)A;
    /* expected calls: one per iteration while refinement is on, which is
     * many; the hint only tells mkl_sparse_optimize that analysis pays */
    st = mkl_sparse_set_mv_hint(A, SCS_MKL_SPARSE_OPERATION_NON_TRANSPOSE,
                                descr, 1000000);
    if (st == SCS_MKL_SPARSE_STATUS_SUCCESS) {
      st = mkl_sparse_optimize(A);
    }
  }
  prod_threads_end(p, prev);
  return st == SCS_MKL_SPARSE_STATUS_SUCCESS ? 0 : -1;
}

static void prod_free(ScsLinSysWork *p) {
  prod_free_handle(p);
  scs_free(p->prod_ind);
  scs_free(p->prod_ptr);
  scs_free(p->prod_r);
  scs_free(p->prod_x0);
  p->prod_ind = p->prod_ptr = SCS_NULL;
  p->prod_r = p->prod_x0 = SCS_NULL;
}

static scs_int prod_setup(ScsLinSysWork *p) {
  scs_int k, nm = p->n_plus_m, nnz = p->kkt->p[nm];
  p->prod_r = (scs_float *)scs_calloc(nm, sizeof(scs_float));
  p->prod_x0 = (scs_float *)scs_calloc(nm, sizeof(scs_float));
  if (!p->prod_r || !p->prod_x0) {
    scs_printf("MKL sparse BLAS: out of memory\n");
    return -1;
  }
  if (sizeof(scs_mkl_int) != sizeof(scs_int)) {
    scs_mkl_int *ptr, *ind;
    if ((scs_int)(scs_mkl_int)nnz != nnz || (scs_int)(scs_mkl_int)nm != nm) {
      /* A size rule, fixed by the problem data: K does not fit MKL's 32-bit
       * sparse BLAS integers (DLONG without BLAS64), while Pardiso, with
       * 64-bit indices, still solves it. So this workspace does not provide
       * the products and SCS computes them, as on master, instead of failing
       * setup. A BLAS64 (ILP64) build provides them at any size. */
      scs_free(p->prod_r);
      scs_free(p->prod_x0);
      p->prod_r = p->prod_x0 = SCS_NULL;
      p->prod_none = 1;
      return 0;
    }
    ptr = (scs_mkl_int *)scs_calloc(nm + 1, sizeof(scs_mkl_int));
    ind = (scs_mkl_int *)scs_calloc(nnz > 0 ? nnz : 1, sizeof(scs_mkl_int));
    p->prod_ptr = ptr;
    p->prod_ind = ind;
    if (!ptr || !ind) {
      scs_printf("MKL sparse BLAS: out of memory\n");
      return -1;
    }
    for (k = 0; k <= nm; ++k) {
      ptr[k] = (scs_mkl_int)p->kkt->p[k];
    }
    for (k = 0; k < nnz; ++k) {
      ind[k] = (scs_mkl_int)p->kkt->i[k];
    }
  }
  if (prod_make_handle(p) != 0) {
    scs_printf("MKL sparse BLAS: creating or analysing the handle of K "
               "failed\n");
    return -1;
  }
  return 0;
}

/* r = b - K x (x = x0 when x is SCS_NULL), r kept in prod_r. Returns
 * ||r||_inf or SCS_LIN_SYS_PRODUCTS_FAILED. */
scs_float scs_lin_sys_residual(ScsLinSysWork *p, const scs_float *b,
                               const scs_float *x) {
  scs_int i, nm = p->n_plus_m;
  scs_float res = 0.;
  int prev, st;
  if (!p->prod_ok) {
    /* the size rule in prod_setup; otherwise never (setup fails without
     * them) */
    return p->prod_none ? SCS_LIN_SYS_PRODUCTS_NOT_PROVIDED
                        : SCS_LIN_SYS_PRODUCTS_FAILED;
  }
  memcpy(p->prod_r, b, nm * sizeof(scs_float));
  prev = prod_threads_begin(p);
  st = SCS_MKL_SPARSE_MV(SCS_MKL_SPARSE_OPERATION_NON_TRANSPOSE, -1.,
                         (scs_mkl_sparse_t)p->prod_k, prod_descr(),
                         x ? x : p->prod_x0, 1., p->prod_r);
  prod_threads_end(p, prev);
  if (st != SCS_MKL_SPARSE_STATUS_SUCCESS) {
    return SCS_LIN_SYS_PRODUCTS_FAILED;
  }
  for (i = 0; i < nm; ++i) {
    res = MAX(res, ABS(p->prod_r[i]));
  }
  return res;
}

void scs_lin_sys_keep(ScsLinSysWork *p, const scs_float *x) {
  if (p->prod_ok) {
    memcpy(p->prod_x0, x, p->n_plus_m * sizeof(scs_float));
  }
}

/* x = x0 + K^{-1} r */
scs_int scs_lin_sys_correct(ScsLinSysWork *p, scs_float *x) {
  scs_int i, status;
  if (!p->prod_ok) {
    return -1;
  }
  memcpy(x, p->prod_r, p->n_plus_m * sizeof(scs_float));
  status = scs_solve_lin_sys(p, x, SCS_NULL, -1.);
  if (status == 0) {
    for (i = 0; i < p->n_plus_m; ++i) {
      x[i] += p->prod_x0[i];
    }
  }
  return status;
}
