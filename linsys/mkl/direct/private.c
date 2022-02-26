#include "private.h"
#include "linsys.h"

#define MKL_INTERFACE_LP64 0
#define MKL_INTERFACE_ILP64 1

#define PARDISO_SYMBOLIC (11)
#define PARDISO_NUMERIC (22)
#define PARDISO_SOLVE (33)
#define PARDISO_CLEANUP (-1)

#ifdef DLONG
#define _PARDISO pardiso_64
#else
#define _PARDISO pardiso
#endif

/* Prototypes for Pardiso functions */
void _PARDISO(void **,           /* pt */
              const scs_int *,   /* maxfct */
              const scs_int *,   /* mnum */
              const scs_int *,   /* mtype */
              const scs_int *,   /* phase */
              const scs_int *,   /* n */
              const scs_float *, /* a */
              const scs_int *,   /* ia */
              const scs_int *,   /* ja */
              scs_int *,         /* perm */
              const scs_int *,   /* nrhs */
              scs_int *,         /* iparam */
              const scs_int *,   /* msglvl */
              scs_float *,       /* b */
              scs_float *,       /* x */
              scs_int *          /* error */
);
scs_int MKL_Set_Interface_Layer(scs_int);

const char *scs_get_lin_sys_method() {
  return "sparse-direct-mkl-pardiso";
}

void scs_free_lin_sys_work(ScsLinSysWork *s) {
  if (s) {
    s->phase = PARDISO_CLEANUP;
    _PARDISO(s->pt, &(s->maxfct), &(s->mnum), &(s->mtype), &(s->phase),
             &(s->n_plus_m), &(s->fdum), s->kkt->p, s->kkt->i, &(s->idum),
             &(s->nrhs), s->iparm, &(s->msglvl), &(s->fdum), &(s->fdum),
             &(s->error));
    if (s->error != 0) {
      scs_printf("Error during MKL Pardiso cleanup: %d", (int)s->error);
    }
    if (s->kkt)
      SCS(cs_spfree)(s->kkt);
    if (s->sol)
      scs_free(s->sol);
    if (s->diag_r_idxs)
      scs_free(s->diag_r_idxs);
    if (s->diag_p)
      scs_free(s->diag_p);
    scs_free(s);
  }
}

ScsLinSysWork *scs_init_lin_sys_work(const ScsMatrix *A, const ScsMatrix *P,
                                     const scs_float *diag_r) {
  scs_int i;
  ScsLinSysWork *s = scs_calloc(1, sizeof(ScsLinSysWork));

  /* TODO: is this necessary with pardiso_64? */
  /* Set MKL interface layer */
#ifdef DLONG
  MKL_Set_Interface_Layer(MKL_INTERFACE_ILP64);
#else
  MKL_Set_Interface_Layer(MKL_INTERFACE_LP64);
#endif

  s->n = A->n;
  s->m = A->m;
  s->n_plus_m = s->n + s->m;

  /* Even though we overwrite rhs with sol pardiso requires the memory */
  s->sol = (scs_float *)scs_malloc(sizeof(scs_float) * s->n_plus_m);
  s->diag_r_idxs = (scs_int *)scs_calloc(s->n_plus_m, sizeof(scs_int));
  s->diag_p = (scs_float *)scs_calloc(s->n, sizeof(scs_float));

  /* MKL pardiso requires upper triangular CSR matrices. The KKT matrix stuffed
   * as CSC lower triangular is equivalent. Pass upper=0. */
  s->kkt = SCS(form_kkt)(A, P, s->diag_p, diag_r, s->diag_r_idxs, 0);
  if (!(s->kkt)) {
    scs_printf("Error in forming KKT matrix");
    scs_free_lin_sys_work(s);
    return SCS_NULL;
  }

  for (i = 0; i < 64; i++) {
    s->iparm[i] = 0; /* Setup Pardiso control parameters */
    s->pt[i] = 0;    /* Initialize the internal solver memory pointer */
  }

  /* Set Pardiso variables */
  s->mtype = -2;         /* Real symmetric indefinite matrix */
  s->nrhs = 1;           /* Number of right hand sides */
  s->maxfct = 1;         /* Maximum number of numerical factorizations */
  s->mnum = 1;           /* Which factorization to use */
  s->error = 0;          /* Initialize error flag */
  s->msglvl = VERBOSITY; /* Printing information */

  /* For all iparm vars see MKL documentation */
  s->iparm[0] = 1;          /* Parsido must inspect iparm */
  s->iparm[1] = 3;          /* Fill-in reordering from OpenMP */
  s->iparm[5] = 1;          /* Write solution into b */
  s->iparm[7] = 0;          /* Automatic iterative refinement calculation */
  s->iparm[9] = 8;          /* Perturb the pivot elements with 1E-8 */
  s->iparm[34] = 1;         /* Use C-style indexing for indices */
  /* s->iparm[36] = -80; */ /* Form block sparse matrices */

#ifdef SFLOAT
  s->iparm[27] = 1; /* 1 is single precision, 0 is double */
#endif

  /* Permutation and symbolic factorization */
  scs_int phase = PARDISO_SYMBOLIC;
  _PARDISO(s->pt, &(s->maxfct), &(s->mnum), &(s->mtype), &phase, &(s->n_plus_m),
           s->kkt->x, s->kkt->p, s->kkt->i, &(s->idum), &(s->nrhs), s->iparm,
           &(s->msglvl), &(s->fdum), &(s->fdum), &(s->error));
  if (s->error != 0) {
    scs_printf("Error during symbolic factorization: %d", (int)s->error);
    scs_free_lin_sys_work(s);
    return SCS_NULL;
  }

  /* Numerical factorization */
  s->phase = PARDISO_NUMERIC;
  _PARDISO(s->pt, &(s->maxfct), &(s->mnum), &(s->mtype), &(s->phase),
           &(s->n_plus_m), s->kkt->x, s->kkt->p, s->kkt->i, &(s->idum),
           &(s->nrhs), s->iparm, &(s->msglvl), &(s->fdum), &(s->fdum),
           &(s->error));

  if (s->error) {
    scs_printf("Error during numerical factorization: %d", (int)s->error);
    scs_free_lin_sys_work(s);
    return SCS_NULL;
  }

  if (s->iparm[21] < s->n) {
    scs_printf("KKT matrix has < n positive eigenvalues. P not PSD.");
    return SCS_NULL;
  }

  return s;
}

/* Returns solution to linear system Ax = b with solution stored in b */
scs_int scs_solve_lin_sys(ScsLinSysWork *s, scs_float *b,
                          const scs_float *warm_start, scs_float tol) {
  /* Back substitution and iterative refinement */
  s->phase = PARDISO_SOLVE;
  _PARDISO(s->pt, &(s->maxfct), &(s->mnum), &(s->mtype), &(s->phase),
           &(s->n_plus_m), s->kkt->x, s->kkt->p, s->kkt->i, &(s->idum),
           &(s->nrhs), s->iparm, &(s->msglvl), b, s->sol, &(s->error));
  if (s->error != 0) {
    scs_printf("Error during linear system solution: %d", (int)s->error);
    return 1;
  }
  return 0;
}

/* Update factorization when R changes */
void scs_update_lin_sys_diag_r(ScsLinSysWork *p, const scs_float *diag_r) {
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
           &(p->n_plus_m), p->kkt->x, p->kkt->p, p->kkt->i, &(p->idum),
           &(p->nrhs), p->iparm, &(p->msglvl), &(p->fdum), &(p->fdum),
           &(p->error));

  if (p->error != 0) {
    scs_printf("Error in PARDISO factorization when updating.\n");
    scs_free_lin_sys_work(p);
    return;
  }
}
