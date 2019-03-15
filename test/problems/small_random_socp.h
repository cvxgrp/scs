#include "rw.h"
#include "scs.h"
#include "util.h"

#define OPT (1.530897)

static const char *small_random_socp(void) {
  ScsData *d;
  ScsCone *k;
  ScsSolution *sol;
  ScsInfo info = {0};
  scs_int exitflag, success;
  scs_float perr, derr;
  const char *filename = "test/data/small_random_socp";

  SCS(read_data)(filename, &d, &k);
  sol = scs_calloc(1, sizeof(ScsSolution));
  exitflag = scs(d, k, sol, &info);

  perr = SCS(dot)(d->c, sol->x, d->n) - OPT;
  derr = -SCS(dot)(d->b, sol->y, d->m) - OPT;

  scs_printf("primal obj error %4e\n", perr);
  scs_printf("dual obj error %4e\n", derr);

  success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;

  SCS(free_data)(d, k);
  SCS(free_sol)(sol);
  mu_assert("small_random_socp: SCS failed to produce outputflag SCS_SOLVED",
            success);
  return 0;
}
