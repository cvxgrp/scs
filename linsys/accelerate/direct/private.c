#include "private.h"

#ifdef DLONG
#error "Apple Accelerate backend does not support DLONG (64-bit integers)"
#endif

const char *scs_get_lin_sys_method(void) {
  return "sparse-direct-apple-accelerate";
}

static scs_int plain_solve(ScsLinSysWork *p, scs_float *b, const scs_float *s,
                           scs_float tol);

/* The inertia check of every factorization: at setup and after each
 * refactorization at an update of R (scs_update_lin_sys_diag_r), as QDLDL,
 * the dense solver, MKL Pardiso and cuDSS check at every factorization.
 *
 * What it checks. K = [R_x + P, A'; A, -R_y] with R > 0 has n positive
 * eigenvalues exactly when R_x + P + A' R_y^{-1} A is positive definite,
 * which holds whenever P is PSD; then the other m are negative. Fewer than n
 * positive pivots means that the factorization is not of a quasi-definite K,
 * so it fails: at setup (P not PSD) SCS stops, and at a scale update SCS
 * goes back to the previous R (scale_update_rollback in src/scs.c). Same
 * rule as QDLDL's: too many positive pivots are not refused (that needs a
 * negative R_y entry, which SCS never makes), and a P that is only mildly
 * indefinite can pass.
 *
 * How the counts are read, without changing the factorization. SCS uses
 * Accelerate's unpivoted LDL^T (SparseFactorizationLDLTUnpivoted), which
 * factorizes S K S = Q L D L' Q' (S a diagonal scaling, Q a permutation)
 * with 1x1 pivots only: D is diagonal, its entries are the pivots, and by
 * Sylvester's law of inertia their signs are those of K's eigenvalues.
 * Accelerate's SparseGetInertia is documented for the pivoted LDL^T
 * (SparseFactorizationLDLTTPP) only, and on the unpivoted one it does not
 * return: it stops the process (SIGTRAP, macOS 26.5). So D's diagonal is
 * read through its subfactor: SparseCreateSubfactor(SparseSubfactorD) wraps
 * the factorization (it extracts nothing), and multiplying it in place with
 * a vector of ones gives the pivots. A pivot that is not a number counts as
 * neither positive nor negative, so it fails too.
 *
 * Cost: n + m numbers of memory (pivots) and, per factorization (so at most
 * once per scale update), the multiply and one pass over the pivots.
 * Measured (Apple M4, one thread, K = [I, A'; A, -I] with A random, about
 * 5 nonzeros per column): 3 us at n + m = 300, where a refactorization took
 * 25 us; 28 us at 6,000 (1.9% of a refactorization); 0.36 ms at 60,000
 * (0.07%).
 *
 * Tested by test_non_psd_p (at setup) and test_factor_rollback (a scale
 * update whose K has n - 1 positive eigenvalues fails, and the restore to the
 * previous R, with the right inertia, passes). Returns 0, or -1 after
 * printing the counts. */
static scs_int check_inertia(ScsLinSysWork *p) {
  scs_int i, pos = 0, neg = 0, n_plus_m = p->n + p->m;
#ifdef SFLOAT
  SparseOpaqueSubfactor_Float d_sub;
  DenseVector_Float dvec = {.count = n_plus_m, .data = p->pivots};
#else
  SparseOpaqueSubfactor_Double d_sub;
  DenseVector_Double dvec = {.count = n_plus_m, .data = p->pivots};
#endif
  for (i = 0; i < n_plus_m; ++i) {
    p->pivots[i] = 1.;
  }
  d_sub = SparseCreateSubfactor(SparseSubfactorD, p->factorization);
  SparseMultiply(d_sub, dvec);
  SparseCleanup(d_sub);
  for (i = 0; i < n_plus_m; ++i) {
    if (p->pivots[i] > 0) {
      ++pos;
    } else if (p->pivots[i] < 0) {
      ++neg;
    }
  }
  (void)neg; /* read only by scs_printf, which NO_PRINTING removes */
  if (pos < p->n) {
    scs_printf("Apple Accelerate: K has %li positive and %li negative "
               "pivots, fewer than n = %li positive (m = %li negative "
               "expected): R_x + P + A' R_y^-1 A is not positive definite "
               "(P not PSD?)\n",
               (long)pos, (long)neg, (long)p->n, (long)p->m);
    return -1;
  }
  return 0;
}

ScsLinSysWork *scs_init_lin_sys_work(const ScsMatrix *A, const ScsMatrix *P,
                                     const scs_float *diag_r) {
  scs_int i, n_plus_m;
  ScsLinSysWork *p = (ScsLinSysWork *)scs_calloc(1, sizeof(ScsLinSysWork));
  if (!p) {
    return SCS_NULL;
  }

  n_plus_m = A->n + A->m;
  p->m = A->m;
  p->n = A->n;

  p->diag_p = (scs_float *)scs_calloc(A->n, sizeof(scs_float));
  p->diag_r_idxs = (scs_int *)scs_calloc(n_plus_m, sizeof(scs_int));
  p->pivots = (scs_float *)scs_calloc(n_plus_m, sizeof(scs_float));
  if (!p->diag_p || !p->diag_r_idxs || !p->pivots) {
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  /* Form upper triangular KKT matrix in CSC format */
  p->kkt = SCS(form_kkt)(A, P, p->diag_p, diag_r, p->diag_r_idxs, 1);
  if (!p->kkt) {
    scs_printf("Error in forming KKT matrix.\n");
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  /* Convert column starts from scs_int (int) to long for Accelerate API */
  p->col_starts = (long *)scs_calloc(n_plus_m + 1, sizeof(long));
  if (!p->col_starts) {
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }
  for (i = 0; i <= n_plus_m; i++) {
    p->col_starts[i] = (long)p->kkt->p[i];
  }

  /* Set up Accelerate sparse matrix wrapper pointing to KKT data */
#ifdef SFLOAT
  p->accel_mat = (SparseMatrix_Float){
#else
  p->accel_mat = (SparseMatrix_Double){
#endif
      .structure = {.rowCount = n_plus_m,
                    .columnCount = n_plus_m,
                    .columnStarts = p->col_starts,
                    .rowIndices = p->kkt->i,
                    .attributes =
                        {
                            .kind = SparseSymmetric,
                            .triangle = SparseUpperTriangle,
                        },
                    .blockSize = 1},
      .data = p->kkt->x};

  /* Perform symbolic and numeric factorization */
  p->factorization =
      SparseFactor(SparseFactorizationLDLTUnpivoted, p->accel_mat);
  if (p->factorization.status != SparseStatusOK) {
    scs_printf("Error in Apple Accelerate LDLt factorization: %d\n",
               (int)p->factorization.status);
    SparseCleanup(p->factorization);
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }
  p->factorizations++;

  /* the inertia check, as after every refactorization */
  if (check_inertia(p) != 0) {
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }

  if (SCS(lin_sys_shared_init)(&p->sh, A, P, diag_r, 0, plain_solve,
                               SCS_NULL, SCS_NULL) != 0) {
    scs_free_lin_sys_work(p);
    return SCS_NULL;
  }
  return p;
}

/* the plain solve K x = b, solution stored in b (s and tol unused) */
static scs_int plain_solve(ScsLinSysWork *p, scs_float *b, const scs_float *s,
                           scs_float tol) {
  scs_int n_plus_m = p->n + p->m;
#ifdef SFLOAT
  DenseVector_Float bvec = {.count = n_plus_m, .data = b};
#else
  DenseVector_Double bvec = {.count = n_plus_m, .data = b};
#endif
  (void)s;
  (void)tol;
  SparseSolve(p->factorization, bvec);
  return 0;
}

/* the warm start, refinement and residual around plain_solve, see linsys.h */
scs_int scs_solve_lin_sys(ScsLinSysWork *p, scs_float *b, const scs_float *s,
                          scs_float tol, scs_int max_refine_steps,
                          scs_float *res) {
  return SCS(lin_sys_shared_solve)(&p->sh, p, b, s, tol, max_refine_steps,
                                   res);
}

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

  /* Refactor with updated values (reuses symbolic factorization) */
  SparseRefactor(p->accel_mat, &(p->factorization));
  if (p->factorization.status != SparseStatusOK) {
    scs_printf("Error in Apple Accelerate LDLt refactorization: %d\n",
               (int)p->factorization.status);
    return -1;
  }
  p->factorizations++;

  /* K with the wrong inertia: fails, and SCS restores the previous R with
   * another call, which refactorizes */
  return check_inertia(p);
}

void scs_free_lin_sys_work(ScsLinSysWork *p) {
  if (p) {
    if (p->factorizations > 0) {
      SparseCleanup(p->factorization);
    }
    SCS(cs_spfree)(p->kkt);
    scs_free(p->diag_p);
    scs_free(p->diag_r_idxs);
    scs_free(p->col_starts);
    scs_free(p->pivots);
    SCS(lin_sys_shared_free)(&p->sh);
    scs_free(p);
  }
}
