#include "scs.h"
#include "run_scs.h"

#ifndef DEMO_PATH
#define DEMO_PATH "../data_sparse"
#endif 

#define NUM_TRIALS 1 
#define RHOX 1e-3

int main(int argc, char **argv)
{
	FILE * fp;
	if(open_file(argc, argv, 1, DEMO_PATH, &fp)==-1) return -1;
	Cone * k = scs_calloc(1,sizeof(Cone));
	Data * d = scs_calloc(1,sizeof(Data));
	if (read_in_data(fp,d,k) == -1){
        printf("Error reading in data, aborting.\n");
        return -1;
    }
	fclose(fp);
	Sol * sol = scs_malloc(sizeof(Sol));
	Info * info = scs_malloc(sizeof(Info));
	int i;
	for (i=0;i<NUM_TRIALS;i++)
	{
		scs(d,k,sol,info);
	}
	//printSol(d,sol,info);
    freeData(d,k);
	freeSol(sol);
	scs_free(info);
	return 0;
}

int read_in_data(FILE * fp,Data * d, Cone * k){
	/* MATRIX IN DATA FILE MUST BE IN COLUMN COMPRESSED FORMAT */
	d->RHO_X = RHOX;
    if(fscanf(fp, "%i", &(d->n)) != 1) return -1;
	if(fscanf(fp, "%i", &(d->m))!= 1) return -1;
    if(fscanf(fp, "%i", &(k->f))!= 1) return -1;
	if(fscanf(fp, "%i", &(k->l))!= 1) return -1;
    if(fscanf(fp, "%i", &(k->qsize))!= 1) return -1;
	if(fscanf(fp, "%i", &(k->ssize))!= 1) return -1;

    // allow arbitrary additional cones, simply add to below:
    int len = 32;
    char s[len];
    if( fgets (s, len, fp) == NULL ) return -1;
    char* token = strtok(s, " ");
    if(token) k->ep = atoi(token); token = strtok(NULL, " ");
    if(token) k->ed = atoi(token); token = strtok(NULL, " ");
    
	//if(fscanf(fp, "%i", &(k->ep))!= 1) return -1;
	//if(fscanf(fp, "%i", &(k->ed))!= 1) return -1;
	
	if(fscanf(fp, "%i", &(d->MAX_ITERS))!= 1) return -1;
	if(fscanf(fp, "%i", &(d->VERBOSE))!= 1) return -1;
	if(fscanf(fp, "%i", &(d->NORMALIZE))!= 1) return -1;
	if(fscanf(fp, "%lf", &(d->ALPH))!= 1) return -1;
	if(fscanf(fp, "%lf", &(d->UNDET_TOL))!= 1) return -1;
	if(fscanf(fp, "%lf", &(d->EPS_ABS))!= 1) return -1;
	if(fscanf(fp, "%i", &(d->Anz))!= 1) return -1;

	k->q = malloc(sizeof(int)*k->qsize);
	for(int i = 0; i < k->qsize; i++)
	{ 
		if(fscanf(fp, "%i", &k->q[i])!= 1) return -1;
	}
    k->s = malloc(sizeof(int)*k->ssize);
    for(int i = 0; i < k->ssize; i++)
    {   
        if(fscanf(fp, "%i", &k->s[i])!= 1) return -1;
    }   
	d->b = malloc(sizeof(double)*d->m);
	for(int i = 0; i < d->m; i++)
	{ 
		if(fscanf(fp, "%lf", &d->b[i])!= 1) return -1;
	}
	d->c = malloc(sizeof(double)*d->n);
	for(int i = 0; i < d->n; i++)
	{ 
		if(fscanf(fp, "%lf", &d->c[i])!= 1) return -1;
	}
	d->Ai = malloc(sizeof(int)*(d->Anz));
	for(int i = 0; i < d->Anz; i++)
	{
		if(fscanf(fp, "%i", &d->Ai[i])!= 1) return -1;
	}
	d->Ap = malloc(sizeof(int)*(d->n+1));
	for(int i = 0; i < d->n+1; i++) 
	{
		if(fscanf(fp, "%i", &d->Ap[i])!= 1) return -1;
	}
	d->Ax = malloc(sizeof(double)*(d->Anz));
	for(int i = 0; i < (d->Anz); i++)
	{
		if(fscanf(fp, "%lf", &d->Ax[i])!= 1) return -1;
	}
	//		fscanf(fp, "%zu", &NNZ);
	//		int *Kr = malloc(sizeof(int)*NNZ);
	//		for(int i = 0; i < NNZ; i++)
	//		{
	//		fscanf(fp, "%i", &Kr[i]);
	//		}
	//		int *Kp=malloc(sizeof(int)*(w->l+1));
	//		for(int i = 0; i < w->l+1; i++)
	//		{
	//		fscanf(fp, "%i", &Kp[i]);
	//		}
	//		double *Kx=malloc(sizeof(double)*NNZ);
	//		for(int i = 0; i < NNZ; i++)
	//		{
	//		fscanf(fp, "%lf", &Kx[i]);
	//		}
    return 0;
}

void freeData(Data * d, Cone * k){
	if(d) {
		if(d->b) scs_free(d->b);
		if(d->c) scs_free(d->c);
		if(d->Ax) scs_free(d->Ax);
		if(d->Ai) scs_free(d->Ai);
		if(d->Ap) scs_free(d->Ap);
		scs_free(d);
	}
	if(k) {
		if(k->q) scs_free(k->q);
		if(k->s) scs_free(k->s);
		scs_free(k);
	}
	d = NULL; k = NULL;
}

void freeSol(Sol *sol){
	if(sol) {
		if(sol->x) scs_free(sol->x);
		if(sol->y) scs_free(sol->y);
		if(sol->s) scs_free(sol->s);
    scs_free(sol);
	}
	sol = NULL;
}


int open_file(int argc, char ** argv, int idx, char * default_file, FILE ** fb) 
{
	if (argc<idx+1){
		printf("Not enough arguments supplied, using %s as default\n", default_file);
	}
	else{
		*fb = fopen(argv[idx], "r");
		if (*fb != NULL) return 0;
		else{
			printf("Couldn't open file %s, using %s as default\n", argv[idx],default_file);
			fclose(*fb);
		}
	}
	*fb = fopen(default_file, "r");
	if (*fb == NULL){
		printf("Couldn't open %s\n",default_file);
		return -1;
	}
	return 0;
}
