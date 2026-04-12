#include "private.h"

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
 *   use_blas64=false (default) -> links mkl-dynamic-lp64-seq,  expects LP64
 *   use_blas64=true            -> links mkl-dynamic-ilp64-seq, expects ILP64
 *
 * See meson.build for the linkage logic. The runtime check in
 * scs_init_lin_sys_work uses BLAS64 to pick the right expected layer.
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

/* MKL_Set_Interface_Layer and MKL_Set_Threading_Layer are exported by mkl_rt
 * (the single dynamic library). All build systems (Makefile, CMake, Meson)
 * should link against mkl_rt so these functions are available. */
int MKL_Set_Interface_Layer(int);
int MKL_Set_Threading_Layer(int);
#define MKL_THREADING_SEQUENTIAL 2

const char *scs_get_lin_sys_method() {
  return "sparse-direct-mkl-pardiso";
}

void scs_free_lin_sys_work(ScsLinSysWork *p) {
  if (p) {
    p->phase = PARDISO_CLEANUP;
    _PARDISO(p->pt, &(p->maxfct), &(p->mnum), &(p->mtype), &(p->phase),
             &(p->n_plus_m), SCS_NULL, p->kkt->p, p->kkt->i, SCS_NULL,
             &(p->nrhs), p->iparm, &(p->msglvl), SCS_NULL, SCS_NULL,
             &(p->error));
    if (p->error != 0) {
      scs_printf("Error during MKL Pardiso cleanup: %d", (int)p->error);
    }
    if (p->kkt)
      SCS(cs_spfree)(p->kkt);
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

  /* If SCS was not compiled with OpenMP, force MKL to use sequential
   * threading. Without this, mkl_rt may try to load Intel OpenMP (libiomp5)
   * which can deadlock if the binary was not linked against it. */
#ifndef _OPENMP
  MKL_Set_Threading_Layer(MKL_THREADING_SEQUENTIAL);
#endif

  /* Enforce the correct MKL interface layer for BLAS/LAPACK calls.
   *
   * The interface layer must match what we linked against:
   *   BLAS64 defined   -> ILP64 (64-bit BLAS integers, mkl-dynamic-ilp64-seq)
   *   BLAS64 undefined -> LP64  (32-bit BLAS integers, mkl-dynamic-lp64-seq)
   *
   * If another library in the process set the wrong layer, our BLAS calls
   * would silently receive the wrong integer width, causing memory corruption.
   *
   * This only protects the BLAS layer. The pardiso_64 entry point is
   * unaffected by the interface layer — it always uses 64-bit integers.
   *
   * All build systems must link against mkl_rt so that this function is
   * available. Do NOT link component libraries (mkl_intel_lp64 etc.)
   * directly — mkl_rt dispatches to the correct component at runtime. */
  {
#ifdef BLAS64
    int expected = MKL_INTERFACE_ILP64;
#else
    int expected = MKL_INTERFACE_LP64;
#endif
    MKL_Set_Interface_Layer(expected);
  }
  p->n = A->n;
  p->m = A->m;
  p->n_plus_m = p->n + p->m;

  p->diag_r_idxs = (scs_int *)scs_calloc(p->n_plus_m, sizeof(scs_int));
  p->diag_p = (scs_float *)scs_calloc(p->n, sizeof(scs_float));

  /* MKL pardiso requires upper triangular CSR matrices. The KKT matrix stuffed
   * as CSC lower triangular is equivalent. Pass upper=0. */
  p->kkt = SCS(form_kkt)(A, P, p->diag_p, diag_r, p->diag_r_idxs, 0);
  if (!(p->kkt)) {
    scs_printf("Error in forming KKT matrix");
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

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

  /* For all iparm vars see MKL documentation */
  p->iparm[0] = 1;          /* Parsido must inspect iparm */
  p->iparm[1] = 3;          /* Fill-in reordering from OpenMP */
  p->iparm[5] = 1;          /* Write solution into b */
  p->iparm[7] = 0;          /* Automatic iterative refinement calculation */
  p->iparm[9] = 13;         /* Perturb the pivot elements with 1E-13 (default) */
  p->iparm[23] = 1;         /* Two-level scheduling for parallel factorization */
  p->iparm[24] = 1;         /* Parallel forward/backward solve */
  p->iparm[34] = 1;         /* Use C-style indexing for indices */
  /* p->iparm[36] = -80; */ /* Form block sparse matrices */

#ifdef SFLOAT
  p->iparm[27] = 1; /* 1 is single precision, 0 is double */
#endif

  /* Permutation and symbolic factorization */
  scs_int phase = PARDISO_SYMBOLIC;
  _PARDISO(p->pt, &(p->maxfct), &(p->mnum), &(p->mtype), &phase, &(p->n_plus_m),
           p->kkt->x, p->kkt->p, p->kkt->i, SCS_NULL, &(p->nrhs), p->iparm,
           &(p->msglvl), SCS_NULL, SCS_NULL, &(p->error));

  if (p->error != 0) {
    scs_printf("Error during symbolic factorization: %d", (int)p->error);
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  /* Numerical factorization */
  p->phase = PARDISO_NUMERIC;
  _PARDISO(p->pt, &(p->maxfct), &(p->mnum), &(p->mtype), &(p->phase),
           &(p->n_plus_m), p->kkt->x, p->kkt->p, p->kkt->i, SCS_NULL,
           &(p->nrhs), p->iparm, &(p->msglvl), SCS_NULL, SCS_NULL, &(p->error));

  if (p->error) {
    scs_printf("Error during numerical factorization: %d", (int)p->error);
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  if (p->iparm[21] < p->n) {
    scs_printf("KKT matrix has < n positive eigenvalues. P not PSD.");
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  return p;
}

/* Returns solution to linear system Ax = b with solution stored in b */
scs_int scs_solve_lin_sys(ScsLinSysWork *p, scs_float *b, const scs_float *ws,
                          scs_float tol) {
  /* Back substitution and iterative refinement */
  p->phase = PARDISO_SOLVE;
  _PARDISO(p->pt, &(p->maxfct), &(p->mnum), &(p->mtype), &(p->phase),
           &(p->n_plus_m), p->kkt->x, p->kkt->p, p->kkt->i, SCS_NULL,
           &(p->nrhs), p->iparm, &(p->msglvl), b, SCS_NULL, &(p->error));
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
  _PARDISO(p->pt, &(p->maxfct), &(p->mnum), &(p->mtype), &(p->phase),
           &(p->n_plus_m), p->kkt->x, p->kkt->p, p->kkt->i, SCS_NULL,
           &(p->nrhs), p->iparm, &(p->msglvl), SCS_NULL, SCS_NULL, &(p->error));

  if (p->error != 0) {
    scs_printf("Error in PARDISO factorization when updating: %d.\n",
               (int)p->error);
    return (scs_int)p->error;
  }
  return 0;
}
