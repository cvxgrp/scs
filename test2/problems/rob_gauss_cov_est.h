#include "amatrix.h"
#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "scs.h"
#include "util.h"

static scs_float Ax[] = {-1.0, -0.25035003149, -0.214684983694, -0.0013696512491, -1.0, -1.0, 0.0941450010301, -0.0929840802212, 0.0143905024814, -1.0, -1.41421356237, -1.0, 0.0941450010301, -0.0929840802212, 0.0143905024814, 1.0, -1.0, -0.0354035554388, -0.0402731435884, -0.151196563215, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.41421356237, 1.0, -1.41421356237, 1.0, -1.41421356237, 1.0, -1.0, 1.0, -1.41421356237, 1.0, -1.41421356237, 1.0, -1.0, 1.0, 1.0, -1.41421356237, 1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
static scs_int Ai[] = {11, 15, 16, 17, 36, 12, 15, 16, 17, 18, 37, 13, 15, 16, 17, 18, 14, 15, 16, 17, 38, 15, 19, 16, 20, 17, 21, 19, 20, 21, 25, 19, 20, 21, 19, 22, 20, 23, 21, 24, 3, 26, 4, 5, 27, 7, 28, 9, 29, 6, 30, 8, 31, 10, 32, 11, 33, 12, 13, 34, 14, 35, 1, 7, 0, 8, 9, 2, 10, 1, 3, 41, 2, 6, 44, 25, 39, 25, 42};
static scs_int Ap[] = {0, 5, 11, 16, 21, 23, 25, 27, 31, 34, 36, 38, 40, 42, 45, 47, 49, 51, 53, 55, 57, 60, 62, 64, 66, 67, 69, 72, 75, 77, 79};
static scs_float b[] = {-0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 1.0, -0.0, -0.0, 1.0, -0.0};
static scs_float c[] ={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

static scs_int m = 45;
static scs_int n = 30;

static scs_int f = 19;
static scs_int l = 7;
static scs_int *q = SCS_NULL;
static scs_int qsize = 0;
static scs_int s[] = {4, 2};
static scs_int ssize = 2;
static scs_int ep = 2;
static scs_int ed = 0;
static scs_float *p = SCS_NULL;
static scs_int psize = 0;

static scs_float opt = -4.8912;

static const char *rob_gauss_cov_est(void) {
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  ScsInfo info = {0};
  scs_int exitflag;
  scs_float perr, derr;
  scs_int success;

  d->m = m;
  d->n = n;
  d->b = b;
  d->c = c;

  d->stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));

  d->A->m = m;
  d->A->n = n;
  d->A->x = Ax;
  d->A->i = Ai;
  d->A->p = Ap;

  k->f = f;
  k->l = l;
  k->q = q;
  k->qsize = qsize;
  k->s = s;
  k->ssize = ssize;
  k->ep = ep;
  k->ed = ed;
  k->p = p;
  k->psize = psize;

  SCS(set_default_settings)(d);

  exitflag = scs(d, k, sol, &info);

  perr = SCS(dot)(d->c, sol->x, d->n) - opt;
  derr = -SCS(dot)(d->b, sol->y, d->m) - opt;

  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;

  SCS(free_sol)(sol);
  scs_free(d->A);
  scs_free(k);
  scs_free(d->stgs);
  scs_free(d);

  mu_assert("rob_gauss_cov_est: SCS failed to produce outputflag SCS_SOLVED",
            success);
  return 0;
}
