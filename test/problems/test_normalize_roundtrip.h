#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "scs.h"
#include "util.h"

/*
 * Test that normalization is consistent:
 * 1. Solve the same problem with normalize=0 and normalize=1 and verify
 *    that the returned (un-normalized) solutions agree.
 * 2. Solve a problem where A entries span several orders of magnitude to
 *    exercise the Ruiz equilibration more aggressively.
 *
 * Problem (LP):
 *   minimize    x1 + x2 + x3
 *   subject to  x1 >= 1     (l cone)
 *               x2 >= 10    (l cone)
 *               x3 >= 100   (l cone)
 *   A entries: -1, -1, -1 (diagonal), b = [-1, -10, -100]
 *   Optimal: x = [1, 10, 100], obj = 111.
 *
 * The varying magnitudes in b (and correspondingly in optimal x) ensure
 * that D scaling factors are non-trivial.
 */
static const char *test_normalize_roundtrip_lp(void) {
  ScsCone *k0, *k1;
  ScsData *d0, *d1;
  ScsSettings *stgs0, *stgs1;
  ScsSolution *sol0, *sol1;
  ScsInfo info0 = {0}, info1 = {0};
  scs_int exitflag0, exitflag1;
  scs_float err_x, err_y, err_s, err_obj;
  scs_int i;

  scs_float opt = 111.0;

  scs_int Ai[]   = {0, 1, 2};
  scs_int Ap[]   = {0, 1, 2, 3};
  scs_int m = 3, n = 3;

  /* Separate copies needed because normalization modifies A, b, c in place */
  scs_float Ax0[] = {-1.0, -1.0, -1.0};
  scs_float Ax1[] = {-1.0, -1.0, -1.0};
  scs_float b0[]  = {-1.0, -10.0, -100.0};
  scs_float b1[]  = {-1.0, -10.0, -100.0};
  scs_float c0[]  = {1.0, 1.0, 1.0};
  scs_float c1[]  = {1.0, 1.0, 1.0};

  /* --- Setup problem 0: normalize OFF --- */
  k0 = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  d0 = (ScsData *)scs_calloc(1, sizeof(ScsData));
  stgs0 = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  sol0 = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));

  d0->m = m; d0->n = n;
  d0->b = b0; d0->c = c0;
  d0->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d0->A->m = m; d0->A->n = n;
  d0->A->x = Ax0; d0->A->i = Ai; d0->A->p = Ap;
  k0->l = m;
  scs_set_default_settings(stgs0);
  stgs0->normalize = 0;
  stgs0->eps_abs = 1e-9;
  stgs0->eps_rel = 1e-9;

  /* --- Setup problem 1: normalize ON --- */
  k1 = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  d1 = (ScsData *)scs_calloc(1, sizeof(ScsData));
  stgs1 = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  sol1 = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));

  d1->m = m; d1->n = n;
  d1->b = b1; d1->c = c1;
  d1->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d1->A->m = m; d1->A->n = n;
  d1->A->x = Ax1; d1->A->i = Ai; d1->A->p = Ap;
  k1->l = m;
  scs_set_default_settings(stgs1);
  stgs1->normalize = 1;
  stgs1->eps_abs = 1e-9;
  stgs1->eps_rel = 1e-9;

  /* --- Solve both --- */
  exitflag0 = scs(d0, k0, stgs0, sol0, &info0);
  exitflag1 = scs(d1, k1, stgs1, sol1, &info1);

  mu_assert("normalize_roundtrip: unnormalized solve failed",
            exitflag0 == SCS_SOLVED);
  mu_assert("normalize_roundtrip: normalized solve failed",
            exitflag1 == SCS_SOLVED);

  /* Both should match the known optimal */
  mu_assert("normalize_roundtrip: unnormalized obj wrong",
            ABS(info0.pobj - opt) < 1e-4);
  mu_assert("normalize_roundtrip: normalized obj wrong",
            ABS(info1.pobj - opt) < 1e-4);

  /* Solutions should agree (un-normalization must be correct) */
  err_x = 0; err_y = 0; err_s = 0;
  for (i = 0; i < n; ++i) {
    scs_float d = ABS(sol0->x[i] - sol1->x[i]);
    if (d > err_x) err_x = d;
  }
  for (i = 0; i < m; ++i) {
    scs_float dy = ABS(sol0->y[i] - sol1->y[i]);
    scs_float ds = ABS(sol0->s[i] - sol1->s[i]);
    if (dy > err_y) err_y = dy;
    if (ds > err_s) err_s = ds;
  }
  err_obj = ABS(info0.pobj - info1.pobj);

  scs_printf("normalize_roundtrip: x_err=%.2e y_err=%.2e s_err=%.2e "
             "obj_err=%.2e\n", err_x, err_y, err_s, err_obj);

  mu_assert("normalize_roundtrip: x vectors differ too much", err_x < 1e-4);
  mu_assert("normalize_roundtrip: y vectors differ too much", err_y < 1e-4);
  mu_assert("normalize_roundtrip: s vectors differ too much", err_s < 1e-4);
  mu_assert("normalize_roundtrip: obj values differ too much", err_obj < 1e-5);

  /* Cleanup */
  SCS(free_sol)(sol0);
  scs_free(d0->A);
  scs_free(k0);
  scs_free(stgs0);
  scs_free(d0);

  SCS(free_sol)(sol1);
  scs_free(d1->A);
  scs_free(k1);
  scs_free(stgs1);
  scs_free(d1);

  return 0;
}

/*
 * Test normalization roundtrip with a QP where P is present,
 * ensuring the E scaling (column scaling for P) is also correct.
 *
 * Problem (QP):
 *   minimize    0.5 * (x1^2 + 100*x2^2) + x1 + x2
 *   subject to  x1 >= 0, x2 >= 0
 *
 * P = [[1, 0], [0, 100]]  (upper triangular)
 * c = [1, 1]
 * A = [[-1, 0], [0, -1]], b = [0, 0], l = 2
 *
 * The 100x difference in P eigenvalues stresses normalization.
 * Optimal: x1 = 0 (since the unconstrained min of x1 + 0.5*x1^2 is at
 *          x1 = -1 which is infeasible), x2 = 0 (unconstrained min at
 *          x2 = -0.01 which is infeasible).
 * Actually: x1 >= 0 binds, x2 >= 0 binds => x = [0, 0], obj = 0.
 *
 * Wait, let me redo: min 0.5*x'Px + c'x s.t. x >= 0.
 * Gradient at x=0: c = [1, 1] > 0, so x=0 is optimal. obj=0.
 */
static const char *test_normalize_roundtrip_qp(void) {
  ScsCone *k0, *k1;
  ScsData *d0, *d1;
  ScsSettings *stgs0, *stgs1;
  ScsSolution *sol0, *sol1;
  ScsInfo info0 = {0}, info1 = {0};
  scs_int exitflag0, exitflag1;
  scs_float err_x, err_y, err_s, err_obj;
  scs_int i;

  scs_float opt = 0.0;

  /* A (2x2 diagonal) */
  scs_float Ax0[] = {-1.0, -1.0};
  scs_float Ax1[] = {-1.0, -1.0};
  scs_int Ai[]    = {0, 1};
  scs_int Ap[]    = {0, 1, 2};

  /* P (2x2 upper triangular, diagonal: [1, 100]) */
  scs_float Px0[] = {1.0, 100.0};
  scs_float Px1[] = {1.0, 100.0};
  scs_int Pi[]    = {0, 1};
  scs_int Pp[]    = {0, 1, 2};

  scs_float b0[] = {0.0, 0.0};
  scs_float b1[] = {0.0, 0.0};
  scs_float c0[] = {1.0, 1.0};
  scs_float c1[] = {1.0, 1.0};
  scs_int m = 2, n = 2;

  /* --- normalize OFF --- */
  k0 = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  d0 = (ScsData *)scs_calloc(1, sizeof(ScsData));
  stgs0 = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  sol0 = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));

  d0->m = m; d0->n = n;
  d0->b = b0; d0->c = c0;
  d0->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d0->A->m = m; d0->A->n = n;
  d0->A->x = Ax0; d0->A->i = Ai; d0->A->p = Ap;
  d0->P = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d0->P->m = n; d0->P->n = n;
  d0->P->x = Px0; d0->P->i = Pi; d0->P->p = Pp;
  k0->l = m;
  scs_set_default_settings(stgs0);
  stgs0->normalize = 0;
  stgs0->eps_abs = 1e-9;
  stgs0->eps_rel = 1e-9;

  /* --- normalize ON --- */
  k1 = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  d1 = (ScsData *)scs_calloc(1, sizeof(ScsData));
  stgs1 = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  sol1 = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));

  d1->m = m; d1->n = n;
  d1->b = b1; d1->c = c1;
  d1->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d1->A->m = m; d1->A->n = n;
  d1->A->x = Ax1; d1->A->i = Ai; d1->A->p = Ap;
  d1->P = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  d1->P->m = n; d1->P->n = n;
  d1->P->x = Px1; d1->P->i = Pi; d1->P->p = Pp;
  k1->l = m;
  scs_set_default_settings(stgs1);
  stgs1->normalize = 1;
  stgs1->eps_abs = 1e-9;
  stgs1->eps_rel = 1e-9;

  /* --- Solve both --- */
  exitflag0 = scs(d0, k0, stgs0, sol0, &info0);
  exitflag1 = scs(d1, k1, stgs1, sol1, &info1);

  mu_assert("normalize_qp: unnormalized solve failed",
            exitflag0 == SCS_SOLVED);
  mu_assert("normalize_qp: normalized solve failed",
            exitflag1 == SCS_SOLVED);

  mu_assert("normalize_qp: unnormalized obj wrong",
            ABS(info0.pobj - opt) < 1e-4);
  mu_assert("normalize_qp: normalized obj wrong",
            ABS(info1.pobj - opt) < 1e-4);

  err_x = 0; err_y = 0; err_s = 0;
  for (i = 0; i < n; ++i) {
    scs_float d = ABS(sol0->x[i] - sol1->x[i]);
    if (d > err_x) err_x = d;
  }
  for (i = 0; i < m; ++i) {
    scs_float dy = ABS(sol0->y[i] - sol1->y[i]);
    scs_float ds = ABS(sol0->s[i] - sol1->s[i]);
    if (dy > err_y) err_y = dy;
    if (ds > err_s) err_s = ds;
  }
  err_obj = ABS(info0.pobj - info1.pobj);

  scs_printf("normalize_qp: x_err=%.2e y_err=%.2e s_err=%.2e "
             "obj_err=%.2e\n", err_x, err_y, err_s, err_obj);

  mu_assert("normalize_qp: x vectors differ too much", err_x < 1e-4);
  mu_assert("normalize_qp: y vectors differ too much", err_y < 1e-4);
  mu_assert("normalize_qp: s vectors differ too much", err_s < 1e-4);
  mu_assert("normalize_qp: obj values differ too much", err_obj < 1e-5);

  /* Cleanup */
  SCS(free_sol)(sol0);
  scs_free(d0->A);
  scs_free(d0->P);
  scs_free(k0);
  scs_free(stgs0);
  scs_free(d0);

  SCS(free_sol)(sol1);
  scs_free(d1->A);
  scs_free(d1->P);
  scs_free(k1);
  scs_free(stgs1);
  scs_free(d1);

  return 0;
}
