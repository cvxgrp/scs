#include <time.h> /* to seed random */

#include "amatrix.h"
#include "problem_utils.h"
#include "scs.h"

/*
 create data for problem:

 minimize 	    c'*x
 subject to 	Ax <=_K b

 where K is a product of zero, linear, and second-order cones. A is a sparse
 matrix in
 CSC format. A is 3n by n with about sqrt(n) nonzeros per column.

 Construct data in such a way that the problem is primal and dual
 feasible and thus bounded.
 */

int main(int argc, char **argv) {
  scs_int n, m, col_nnz, nnz, i, q_total, q_num_rows, max_q;
  ScsCone *k;
  ScsData *d;
  ScsSolution *sol, *opt_sol;
  ScsInfo info = {0};
  scs_float p_f, p_l;
  scs_float geo_mean, sum_log_iter = 0.;
  scs_int num_failures = 0, max_seeds = 100;
  int seed = 0;

  /* default parameters */
  p_f = 0.1;
  p_l = 0.3;
  seed = time(SCS_NULL);

  switch (argc) {
    /* no break */
    case 5:
      max_seeds = atoi(argv[4]);
    case 4:
      p_f = atof(argv[2]);
      p_l = atof(argv[3]);
    /* no break */
    case 2:
      n = atoi(argv[1]);
      break;
    default:
      scs_printf(
          "usage:\t%s n p_f p_l max_seeds\n"
          "\tcreates an SOCP with n variables where p_f fraction of "
          "rows correspond\n"
          "\tto equality constraints, p_l fraction of rows correspond "
          "to LP constraints,\n"
          "\tand the remaining percentage of rows are involved in "
          "second-order\n"
          "\tcone constraints.\n",
          argv[0]);
      scs_printf(
          "\nusage:\t%s n p_f p_l\n"
          "\tdefaults the seed to the system time\n",
          argv[0]);
      scs_printf(
          "\nusage:\t%s n\n"
          "\tdefaults to using p_f = 0.1 and p_l = 0.3\n",
          argv[0]);
      return 0;
  }

  k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  d->stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  opt_sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));

  for (seed = 0; seed <= max_seeds; seed++) {
    geo_mean = sum_log_iter / seed;
    geo_mean = exp(geo_mean);

    if (seed * 100 % max_seeds == 0 && (seed * 100 / max_seeds) % 10 == 0) {
      scs_printf("\rdone: %ld%% | iter geomean: %f | %% fail %f\t\t",
                 (long)(seed * 100 / max_seeds), geo_mean,
                 (100. * num_failures) / seed);
      fflush(stdout);
    }
    srand(seed);

    m = 3 * n;
    col_nnz = (int)ceil(sqrt(n));
    nnz = n * col_nnz;

    max_q = (scs_int)ceil(3 * n / log(3 * n));

    if (p_f + p_l > 1.0) {
      printf("error: p_f + p_l > 1.0!\n");
      return 1;
    }

    k->f = (scs_int)floor(3 * n * p_f);
    k->l = (scs_int)floor(3 * n * p_l);

    k->qsize = 0;
    q_num_rows = 3 * n - k->f - k->l;
    k->q = (scs_int *)scs_malloc(q_num_rows * sizeof(scs_int));

    while (q_num_rows > max_q) {
      int size;
      size = (rand() % max_q) + 1;
      k->q[k->qsize] = size;
      k->qsize++;
      q_num_rows -= size;
    }
    if (q_num_rows > 0) {
      k->q[k->qsize] = q_num_rows;
      k->qsize++;
    }

    q_total = 0;
    for (i = 0; i < k->qsize; i++) {
      q_total += k->q[i];
    }

    k->s = SCS_NULL;
    k->ssize = 0;
    k->ep = 0;
    k->ed = 0;

    /* set up SCS structures */
    d->m = m;
    d->n = n;
    gen_random_prob_data(nnz, col_nnz, d, k, opt_sol);
    SCS(set_default_settings)(d);

    d->stgs->verbose = 0;

    /* solve! */
    scs(d, k, sol, &info);

    if (info.status_val == 2 || info.status_val < 0) {
      /* scs_printf("%i\n", seed); */
      num_failures += 1;
    }

    sum_log_iter += log(info.iter);
  }
  SCS(free_data)(d, k);
  SCS(free_sol)(sol);
  SCS(free_sol)(opt_sol);

  geo_mean = sum_log_iter / seed;
  geo_mean = exp(geo_mean);

  scs_printf("\rdone: %i%% | iter geomean: %f | %% fail %f\n",
             100, geo_mean,
             (100. * num_failures) / max_seeds);

  return 0;
}
