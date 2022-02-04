#include "problem_utils.h"
#include "scs.h"
#include "scs_matrix.h"
#include <time.h> /* to seed random */

/*
 create data for problem:

 minimize 	    c'*x
 subject to 	Ax <=_K b

 where K is a product of zero, linear, and second-order cones. A is a sparse
 matrix in
 CSC format. A is factor * n by n with about sqrt(n) nonzeros per column.

 Construct data in such a way that the problem is primal and dual
 feasible and thus bounded.
 */

int main(int argc, char **argv) {
  scs_int n, m, col_nnz, nnz, i, q_total, q_num_rows, max_q;
  ScsCone *k;
  ScsData *d;
  ScsSettings *stgs;
  ScsSolution *sol, *opt_sol;
  ScsInfo info = {0};
  scs_float p_f, p_l;
  scs_int factor = 4;
  int seed = 0;

  /* default parameters */
  p_f = 0.1;
  p_l = 0.3;
  seed = time(SCS_NULL);

  switch (argc) {
  case 5:
    seed = atoi(argv[4]);
  /* no break */
  case 4:
    p_f = atof(argv[2]);
    p_l = atof(argv[3]);
  /* no break */
  case 2:
    n = atoi(argv[1]);
    break;
  default:
    scs_printf("usage:\t%s n p_f p_l s\n"
               "\tcreates an SOCP with n variables where p_f fraction of "
               "rows correspond\n"
               "\tto equality constraints, p_l fraction of rows correspond "
               "to LP constraints,\n"
               "\tand the remaining percentage of rows are involved in "
               "second-order\n"
               "\tcone constraints. the random number generator is seeded "
               "with s.\n"
               "\tnote that p_f + p_l should be less than or equal to 1, "
               "and that\n"
               "\tp_f should be less than .33, since that corresponds to "
               "as many equality\n"
               "\tconstraints as variables.\n",
               argv[0]);
    scs_printf("\nusage:\t%s n p_f p_l\n"
               "\tdefaults the seed to the system time\n",
               argv[0]);
    scs_printf("\nusage:\t%s n\n"
               "\tdefaults to using p_f = 0.1 and p_l = 0.3\n",
               argv[0]);
    return 0;
  }
  scs_printf("seed : %i\n", seed);

  k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  opt_sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));

  m = factor * n;
  col_nnz = (int)ceil(sqrt(n));
  nnz = n * col_nnz;

  max_q = (scs_int)ceil(factor * n / log(factor * n));

  if (p_f + p_l > 1.0) {
    printf("error: p_f + p_l > 1.0!\n");
    return 1;
  }

  k->z = (scs_int)floor(factor * n * p_f);
  k->l = (scs_int)floor(factor * n * p_l);

  k->qsize = 0;
  q_num_rows = factor * n - k->z - k->l;
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
  k->bsize = 0;
  k->bu = SCS_NULL;
  k->bl = SCS_NULL;

  scs_printf("\nA is %ld by %ld, with %ld nonzeros per column.\n", (long)m,
             (long)n, (long)col_nnz);
  scs_printf("A has %ld nonzeros (%f%% dense).\n", (long)nnz,
             100 * (scs_float)col_nnz / m);
  scs_printf("Nonzeros of A take %f GB of storage.\n",
             ((scs_float)nnz * sizeof(scs_float)) / POWF(2, 30));
  scs_printf("Row idxs of A take %f GB of storage.\n",
             ((scs_float)nnz * sizeof(scs_int)) / POWF(2, 30));
  scs_printf("Col ptrs of A take %f GB of storage.\n\n",
             ((scs_float)n * sizeof(scs_int)) / POWF(2, 30));

  printf("ScsCone information:\n");
  printf("Zero cone rows: %ld\n", (long)k->z);
  printf("LP cone rows: %ld\n", (long)k->l);
  printf("Number of second-order cones: %ld, covering %ld rows, with sizes\n[",
         (long)k->qsize, (long)q_total);
  for (i = 0; i < k->qsize; i++) {
    printf("%ld, ", (long)k->q[i]);
  }
  printf("]\n");
  printf("Number of rows covered is %ld out of %ld.\n\n",
         (long)(q_total + k->z + k->l), (long)m);

  /* set up SCS structures */
  d->m = m;
  d->n = n;
  gen_random_prob_data(nnz, col_nnz, d, k, opt_sol, seed);
  scs_set_default_settings(stgs);

  /* stgs->write_data_filename = "random_socp_prob"; */

  scs_printf("true pri opt = %4f\n", SCS(dot)(d->c, opt_sol->x, d->n));
  scs_printf("true dua opt = %4f\n", -SCS(dot)(d->b, opt_sol->y, d->m));
  /* solve! */
  scs(d, k, stgs, sol, &info);
  scs_printf("true pri opt = %4f\n", SCS(dot)(d->c, opt_sol->x, d->n));
  scs_printf("true dua opt = %4f\n", -SCS(dot)(d->b, opt_sol->y, d->m));

  if (sol->x) {
    scs_printf("scs pri obj= %4f\n", SCS(dot)(d->c, sol->x, d->n));
  }
  if (sol->y) {
    scs_printf("scs dua obj = %4f\n", -SCS(dot)(d->b, sol->y, d->m));
  }

  SCS(free_data)(d);
  SCS(free_cone)(k);
  SCS(free_sol)(sol);
  SCS(free_sol)(opt_sol);
  scs_free(stgs);

  return 0;
}
