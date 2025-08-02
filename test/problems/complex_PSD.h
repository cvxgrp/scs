#include "glbopts.h"
#include "linalg.h"
#include "minunit.h"
#include "problem_utils.h"
#include "rw.h"
#include "scs.h"
#include "scs_matrix.h"
#include "util.h"

static const char *complex_PSD(void)
{
    ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
    ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
    ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
    ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
    ScsInfo info = {0};
    scs_int exitflag;
    scs_float perr, derr;
    scs_int success;
    const char *fail;

    /* data */
    scs_float Ax[] = {1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, 1.0, -1.0};
    scs_int Ai[] = {0, 1, 2, 3, 4, 5, 0, 6, 7, 8, 0, 9};
    scs_int Ap[] = {0, 2, 3, 4, 5, 6, 8, 9, 10, 12};
    scs_float b[] = {1.0, 0., 0., 0., 0., 0., 0., 0., 0., 0.};
    scs_float c[] = {1, 2 * sqrt(2), 3 * sqrt(2), 4 * sqrt(2), 5 * sqrt(2), 6,
                     7 * sqrt(2), -8 * sqrt(2), 9};

    scs_int m = 10;
    scs_int n = 9;
    scs_int cs[] = {3};
    scs_int cssize = 1;
    k->z = 1;
    k->cs = cs;
    k->cssize = cssize;

    // computed offline
    scs_float opt = -5.228930;
    /* end data */

    d->m = m;
    d->n = n;
    d->b = b;
    d->c = c;

    d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
    d->A->m = m;
    d->A->n = n;
    d->A->x = Ax;
    d->A->i = Ai;
    d->A->p = Ap;

    scs_set_default_settings(stgs);
    stgs->eps_abs = 1e-7;
    stgs->eps_rel = 1e-7;
    stgs->eps_infeas = 1e-9;

    exitflag = scs(d, k, stgs, sol, &info);

    perr = SCS(dot)(d->c, sol->x, d->n) - opt;
    derr = -SCS(dot)(d->b, sol->y, d->m) - opt;

    scs_printf("primal obj error %4e\n", perr);
    scs_printf("dual obj error %4e\n", derr);

    success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;

    mu_assert("complex_PSD: SCS failed to produce outputflag SCS_SOLVED", success);

    fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);
    if (fail)
        return fail;

    /* kill data */
    scs_free(d->A);
    scs_free(k);
    scs_free(stgs);
    scs_free(d);
    SCS(free_sol)(sol);

    return fail;
}
