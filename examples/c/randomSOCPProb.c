#include "scs.h"
#include "../linsys/amatrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h> /* pow function */
#include <time.h> /* to seed random */
#include "../../include/cones.h"
#include "linAlg.h"

pfloat rand_pfloat(void);

/*
 create data for problem:

 minimize 	    c'*x
 subject to 	Ax <=_K b

 where K is a product of zero, linear, and second-order cones. A is a sparse matrix in
 CSC format. A is 3n by n with about sqrt(n) nonzeros per column.

 Construct data in such a way that the problem is primal and dual
 feasible and thus bounded.
 */

int main(int argc, char **argv) {
	idxint n, m, col_nnz, nnz, i, j, r;
	int seed = 0;
	AMatrix A;
	pfloat * c, *b;
	pfloat * z, *s, *y, *x;
	Cone k;
	Data d;
	Sol sol = { 0 };
	Info info = { 0 };
	idxint q_num_rows;
	int q_total = 0;

	pfloat p_f, p_l;
	idxint max_q;

	/* default parameters */
	p_f = 0.1;
	p_l = 0.3;
	seed = time(NULL);

	switch(argc){
		case 5:
			seed = atoi(argv[4]);
		case 4:
			p_f = atof(argv[2]);
			p_l = atof(argv[3]);
		case 2:
			n = atoi(argv[1]);
			break;
		default:
			scs_printf("usage:\t%s n p_f p_l s\n\t"
				"creates an SOCP with n variables where p_f fraction of rows correspond "
				"to equality constraints, p_l fraction of rows correspond to LP constraints, "
				"and the remaining percentage of rows are involved in second-order "
				"cone constraints. the random number generator is seeded with s. "
				"note that p_f + p_l should be less than or equal to 1, and that "
				"p_f should be less than .33, since that corresponds to as many equality "
				"constraints as variables, which would not be an interesting problem.\n",argv[0]);
			scs_printf("usage:\t%s n p_f p_l\n\t"
				"defaults the seed to the system time\n",argv[0]);
			scs_printf("usage:\t%s n\n\t"
				"defaults p_f = 0.1 and p_l = 0.3\n",argv[0]);
			return 0;
	}
	srand(seed);
	scs_printf("seed : %i\n", seed);

	m = 3 * n;
	col_nnz = (int) ceil(sqrt(n));
	nnz = n * col_nnz;


	max_q = (idxint) ceil(3*n/log(3*n)); 

	if(p_f + p_l > 1.0){
		printf("error: p_f + p_l > 1.0!\n");
		return 1;
	}

	k.f = (idxint)floor(3*n*p_f);
	k.l = (idxint)floor(3*n*p_l);

	k.qsize = 0;
	q_num_rows = 3*n - k.f - k.l;
	k.q = scs_malloc(q_num_rows*sizeof(idxint));

	while(q_num_rows > max_q){
		int size;
		size = (rand() % max_q ) + 1;
		k.q[k.qsize] = size;
		k.qsize++;
		q_num_rows -= size;
	}
	if(q_num_rows > 0){
		k.q[k.qsize] = q_num_rows;
		k.qsize++;
	}

	for(i=0;i<k.qsize;i++){
		q_total += k.q[i];
	}

	k.s = NULL;
	k.ssize = 0;
	k.ep = 0;
	k.ed = 0;

	scs_printf("\nA is %d by %d, with %d nonzeros per column.\n", m, n, col_nnz);
	scs_printf("A has %d nonzeros (%f%% dense).\n", nnz, 100 * (pfloat) col_nnz / m);
	scs_printf("Nonzeros of A take %f GB of storage.\n\n", ((pfloat) nnz * sizeof(pfloat)) / pow(2, 30));

	printf("Cone information:\n");
	printf("Zero cone rows: %d\n", k.f);
	printf("LP cone rows: %d\n", k.l);
	printf("Number of second-order cones: %d, covering %d rows, with sizes\n[",  k.qsize, q_total);
	for(i=0;i<k.qsize;i++){
		printf("%d, ", k.q[i]);
	}
	printf("]\n");
	printf("Number of rows covered is %d out of %d.\n\n", q_total+k.f+k.l, m);

	A.i = scs_malloc(nnz * sizeof(idxint));
	A.p = scs_malloc((n + 1) * sizeof(idxint));
	A.x = scs_malloc(nnz * sizeof(pfloat));
	c = scs_malloc(n * sizeof(pfloat));
	b = scs_malloc(m * sizeof(pfloat));
	z = scs_malloc(m * sizeof(pfloat));
	y = scs_malloc(m * sizeof(pfloat));
	s = scs_malloc(m * sizeof(pfloat));
	x = scs_malloc(n * sizeof(pfloat));

	/* y, s >= 0 and y'*s = 0 */
	for (i = 0; i < m; i++) {
		z[i] = rand_pfloat();
		y[i] = z[i];
	}

	projDualCone(y, &k, -1);

	for (i = 0; i < m; i++) {
		s[i] = y[i] - z[i];
		b[i] = s[i];
	}

	for (i = 0; i < n; i++) {
		x[i] = rand_pfloat();
		c[i] = 0.0;
	}

	/* 	c = -A'*y
	 b = A*x + s
	 */
	A.p[0] = 0;
	for (j = 0; j < n; j++) { /* column */
		for (r = 0; r < col_nnz; r++) { /* row index */
			i = rand() % m; /* row */
			A.x[r + j * col_nnz] = rand_pfloat();
			A.i[r + j * col_nnz] = i;

			b[i] += A.x[r + j * col_nnz] * x[j];

			c[j] -= A.x[r + j * col_nnz] * y[i];
		}
		A.p[j + 1] = (j + 1) * col_nnz;
	}

	/* set up SCS structures */
	d.m = m;
	d.n = n;
	d.A = &A;
	d.b = b;
	d.c = c;
	d.MAX_ITERS = 2500; /* maximum iterations to take: 2500 */
	d.EPS = 1e-3; /* convergence tolerance: 1e-3 */
	d.ALPHA = 1.8; /* relaxation parameter: 1.8 */
	d.RHO_X = 1e-3; /* x equality constraint scaling: 1e-3 */
	d.SCALE = 1; /* if normalized, rescales by this factor: 1 */
	d.CG_RATE = 1.5; /* for indirect, tolerance goes down like (1/iter)^CG_RATE: 1.5 */
	d.VERBOSE = 1; /* boolean, write out progress: 1 */
	d.NORMALIZE = 1; /* boolean, heuristic data rescaling: 1 */
	d.WARM_START = 0;

	scs(&d, &k, &sol, &info);

	scs_printf("true pri opt = %4f\n", innerProd(c, x, d.n));
	scs_printf("true dua opt = %4f\n", -innerProd(b, y, d.m));

	scs_free(A.i);
	scs_free(A.p);
	scs_free(A.x);
	scs_free(c);
	scs_free(b);
	scs_free(z);
	scs_free(y);
	scs_free(s);
	scs_free(x);

	scs_free(k.q);

	scs_free(sol.x);
	scs_free(sol.y);
	scs_free(sol.s);

	return 0;
}

/* uniform random number in [-1,1] */
pfloat rand_pfloat(void) {
	return 2 * (((pfloat) rand()) / RAND_MAX) - 1;
}
