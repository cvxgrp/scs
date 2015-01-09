#include "scs.h"
#include "linsys/amatrix.h"
#include "problemUtils.h"

#ifndef DEMO_PATH
#define DEMO_PATH "examples/raw/demo_data"
#endif

#define NUM_TRIALS (5)
#define RHOX (1e-3)
#define TEST_WARM_START (1)

scs_int read_in_data(FILE * fp, Data * d, Cone * k);
scs_int open_file(scs_int argc, char ** argv, scs_int idx, char * default_file, FILE ** fb);
/* void printSol(Data * d, Sol * sol, Info * info); */

int main(int argc, char **argv) {
	FILE * fp;
	Cone * k;
	Data * d;
	Work * w;
	Sol * sol;
	Info info = { 0 };
	scs_int i;

	if (open_file(argc, argv, 1, DEMO_PATH, &fp) < 0)
		return -1;

	k = scs_calloc(1, sizeof(Cone));
	d = scs_calloc(1, sizeof(Data));
	Settings * stgs = scs_calloc(1, sizeof(Settings));
	sol = scs_calloc(1, sizeof(Sol));

	if (read_in_data(fp, d, k) == -1) {
		printf("Error reading in data, aborting.\n");
		return -1;
	}
	fclose(fp);
	scs_printf("solve once using scs\n");
	stgs->cg_rate = 2;
	d->stgs = stgs;
	scs(d, k, sol, &info);
	if (TEST_WARM_START) {
		scs_printf("solve %i times with warm-start and (if applicable) factorization caching.\n", NUM_TRIALS);
		/* warm starts stored in Sol */
		w = scs_init(d, k, &info);
		if (w) {
			for (i = 0; i < NUM_TRIALS; i++) {
				/* perturb b and c */
				perturbVector(d->b, d->m);
				perturbVector(d->c, d->n);
				stgs->warm_start = 1;
				stgs->cg_rate = 4;
				scs_solve(w, d, k, sol, &info);
				stgs->warm_start = 0;
				stgs->cg_rate = 2;
				scs_solve(w, d, k, sol, &info);
			}
		}
		scs_printf("finished\n");
		scs_finish(w);
	}
	freeData(d, k);
	freeSol(sol);
	return 0;
}

scs_int read_in_data(FILE * fp, Data * d, Cone * k) {
	/* MATRIX IN DATA FILE MUST BE IN COLUMN COMPRESSED FORMAT */
#define LEN64 64 /* variable-size arrays not allowed in ansi */
	char s[LEN64], *token;
	scs_int i, Anz;
	AMatrix * A;
	Settings * stgs = (Settings *) d->stgs;
	stgs->rho_x = RHOX;
	stgs->warm_start = 0;
	stgs->scale = 1;
	if (fscanf(fp, INTRW, &(d->n)) != 1)
		return -1;
	if (fscanf(fp, INTRW, &(d->m)) != 1)
		return -1;
	if (fscanf(fp, INTRW, &(k->f)) != 1)
		return -1;
	if (fscanf(fp, INTRW, &(k->l)) != 1)
		return -1;
	if (fscanf(fp, INTRW, &(k->qsize)) != 1)
		return -1;
	if (fscanf(fp, INTRW, &(k->ssize)) != 1)
		return -1;

	/* allow arbitrary additional cones, simply add to below: */
	if (fgets(s, LEN64, fp) == NULL)
		return -1;
	token = strtok(s, " ");
	if (token)
		k->ep = atoi(token);
	token = strtok(NULL, " ");
	if (token)
		k->ed = atoi(token);
	token = strtok(NULL, " ");

	/*if(fscanf(fp, INTRW, &(k->ep))!= 1) return -1; */
	/*if(fscanf(fp, INTRW, &(k->ed))!= 1) return -1; */

	if (fscanf(fp, INTRW, &(stgs->max_iters)) != 1)
		return -1;
	if (fscanf(fp, INTRW, &(stgs->verbose)) != 1)
		return -1;
	if (fscanf(fp, INTRW, &(stgs->normalize)) != 1)
		return -1;
	if (fscanf(fp, FLOATRW, &(stgs->alpha)) != 1)
		return -1;
	if (fscanf(fp, FLOATRW, &(stgs->eps)) != 1)
		return -1;
	k->q = malloc(sizeof(scs_int) * k->qsize);
	for (i = 0; i < k->qsize; i++) {
		if (fscanf(fp, INTRW, &k->q[i]) != 1)
			return -1;
	}
	k->s = malloc(sizeof(scs_int) * k->ssize);
	for (i = 0; i < k->ssize; i++) {
		if (fscanf(fp, INTRW, &k->s[i]) != 1)
			return -1;
	}
	d->b = malloc(sizeof(scs_float) * d->m);
	for (i = 0; i < d->m; i++) {
		if (fscanf(fp, FLOATRW, &d->b[i]) != 1)
			return -1;
	}
	d->c = malloc(sizeof(scs_float) * d->n);
	for (i = 0; i < d->n; i++) {
		if (fscanf(fp, FLOATRW, &d->c[i]) != 1)
			return -1;
	}
	A = malloc(sizeof(AMatrix));
	A->p = malloc(sizeof(scs_int) * (d->n + 1));
	for (i = 0; i < d->n + 1; i++) {
		if (fscanf(fp, INTRW, &A->p[i]) != 1)
			return -1;
	}
	Anz = A->p[d->n];
	A->i = malloc(sizeof(scs_int) * Anz);
	for (i = 0; i < Anz; i++) {
		if (fscanf(fp, INTRW, &A->i[i]) != 1)
			return -1;
	}
	A->x = malloc(sizeof(scs_float) * Anz);
	for (i = 0; i < Anz; i++) {
		if (fscanf(fp, FLOATRW, &A->x[i]) != 1)
			return -1;
	}
	d->A = A;
	return 0;
}

scs_int open_file(scs_int argc, char ** argv, scs_int idx, char * default_file, FILE ** fb) {
	if (argc < idx + 1) {
		printf("Not enough arguments supplied, using %s as default\n", default_file);
	} else {
		*fb = fopen(argv[idx], "r");
		if (*fb != NULL)
			return 0;
		else {
			printf("Couldn't open file %s, using %s as default\n", argv[idx], default_file);
			fclose(*fb);
		}
	}
	*fb = fopen(default_file, "r");
	if (*fb == NULL) {
		printf("Couldn't open %s\n", default_file);
		return -1;
	}
	return 0;
}
