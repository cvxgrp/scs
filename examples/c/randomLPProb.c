#include "scs.h"
#include "../linsys/amatrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h> /* pow function */
#include <time.h> /* to seed random */

pfloat rand_pfloat(void);
pfloat pos(pfloat x);

/*
create data for problem:

minimize 	c'*x
subject to 	Ax <= b

where the constraints are linear inequalities. A is a sparse matrix in
CSC format. A is 3n by n with about sqrt(n) nonzeros per column.

Construct data in such a way that the problem is primal and dual
feasible and thus bounded.
*/

int main(int argc, char **argv) {
	idxint n, m, col_nnz, nnz, i, j, r;
	AMatrix A;
	pfloat * c, * b;
	pfloat * z, * s, * y, * x;
	Cone k;
	Data d;
	Sol sol = { 0 };
	Info info = { 0 };

	if(argc < 2){
		printf("usage: \"%s [n]\", where n is the number of columns in the matrix A."\
			"Will seed the random number generator with the system time.\n",
			argv[0]);
		printf("\"%s [n] [s]\" will form A with n columns and seed the RNG with s.\n",argv[0]);

		return 0;
	}else if(argc == 2){
		srand(time(NULL));
	}else if(argc == 3){
		srand(atoi(argv[2]));
	}

	n = atoi(argv[1]);
	m = 3*n;
	col_nnz = (int)ceil(sqrt(n));
	nnz = n*col_nnz;


	printf("\n\nA is %d by %d, with %d nonzeros per column.\n", m,n, col_nnz);
	printf("A has %d nonzeros (%f %% dense).\n", nnz, 100*(pfloat)col_nnz/m);
	printf("Nonzeros of A take %f GB of storage.\n\n", ((pfloat)nnz*sizeof(pfloat))/pow(2,30));

	A.i = scs_malloc( nnz * sizeof(idxint));
	A.p = scs_malloc( (n+1)* sizeof(idxint));
	A.x = scs_malloc( nnz * sizeof(pfloat));
	c = scs_malloc( n * sizeof(pfloat));
	b = scs_malloc( m * sizeof(pfloat));
	z = scs_malloc( m * sizeof(pfloat));
	y = scs_malloc( m * sizeof(pfloat));
	s = scs_malloc( m * sizeof(pfloat));
	x = scs_malloc( n * sizeof(pfloat));


	/* y, s >= 0 and y'*s = 0 */
	for(i=0;i<m;i++){
		z[i] = rand_pfloat();
		y[i] = pos(z[i]);
		s[i] = y[i] - z[i];
		b[i] = 0.0;
	}
	for(i=0;i<n;i++){
		x[i] = rand_pfloat();
		c[i] = 0.0;
	}


	/* 	c = -A'*y
		b = A*x + s
	*/
	A.p[0] = 0;
	for(j=0;j<n;j++){
		for(r=0;r<col_nnz;r++){
			i = rand() % m;
			A.x[r+j*col_nnz] = rand_pfloat();
			A.i[r+j*col_nnz] = i;

			b[i] += A.x[r+j*col_nnz]*x[j] + s[i];

			c[j] -= A.x[r+j*col_nnz]*y[i];
		}
		A.p[j+1] = (j+1)*col_nnz;
	}


	/* set up SCS structures */
	k.f = 0;
	k.l = m;
	k.q = NULL;
	k.qsize = 0;
	k.s = NULL;
	k.ssize = 0;
	k.ep = 0;
	k.ed = 0;

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

	scs_free(A.i);
	scs_free(A.p);
	scs_free(A.x);
	scs_free(c);
	scs_free(b);
	scs_free(z);
	scs_free(y);
	scs_free(s);
	scs_free(x);

	return 0;
}

/* uniform random number in [-1,1] */
pfloat rand_pfloat(void){
	return 2*(((pfloat)rand())/RAND_MAX) - 1;
}

pfloat pos(pfloat x){
	return x > 0 ? x : 0;
}