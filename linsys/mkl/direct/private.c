#include "private.h"

#define PARDISO_SYMBOLIC (11)
#define PARDISO_NUMERIC (22)
#define PARDISO_SOLVE (33)
#define PARDISO_CLEANUP (-1)

/* TODO: is it necessary to use pardiso_64 and MKL_Set_Interface_Layer ? */
/*
#define MKL_INTERFACE_LP64 0
#define MKL_INTERFACE_ILP64 1
*/
#ifdef DLONG
#define _PARDISO pardiso_64
#else
#define _PARDISO pardiso
#endif

/* Prototypes for Pardiso functions */
void _PARDISO(void **pt, const scs_int *maxfct, const scs_int *mnum,
              const scs_int *mtype, const scs_int *phase, const scs_int *n,
              const scs_float *a, const scs_int *ia, const scs_int *ja,
              scs_int *perm, const scs_int *nrhs, scs_int *iparm,
              const scs_int *msglvl, scs_float *b, scs_float *x,
              scs_int *error);
/* scs_int MKL_Set_Interface_Layer(scs_int); */

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
    if (p->sol)
      scs_free(p->sol);
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

  /* TODO: is this necessary with pardiso_64? */
  /* Set MKL interface layer */
  /*
#ifdef DLONG
  MKL_Set_Interface_Layer(MKL_INTERFACE_ILP64);
#else
  MKL_Set_Interface_Layer(MKL_INTERFACE_LP64);
#endif
  */
  p->n = A->n;
  p->m = A->m;
  p->n_plus_m = p->n + p->m;

  /* Even though we overwrite rhs with sol pardiso requires the memory */
  p->sol = (scs_float *)scs_malloc(sizeof(scs_float) * p->n_plus_m);
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
  p->iparm[9] = 8;          /* Perturb the pivot elements with 1E-8 */
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
           &(p->nrhs), p->iparm, &(p->msglvl), b, p->sol, &(p->error));
  if (p->error != 0) {
    scs_printf("Error during linear system solution: %d", (int)p->error);
  }
  return p->error;
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
           &(p->n_plus_m), p->kkt->x, p->kkt->p, p->kkt->i, SCS_NULL,
           &(p->nrhs), p->iparm, &(p->msglvl), SCS_NULL, SCS_NULL, &(p->error));

  if (p->error != 0) {
    scs_printf("Error in PARDISO factorization when updating: %d.\n",
               (int)p->error);
    scs_free_lin_sys_work(p);
  }
}
