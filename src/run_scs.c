#include "scs.h"

int main(int argc, char **argv);
idxint read_in_data(FILE * fp,Data * d, Cone * k);
idxint open_file(idxint argc, char ** argv, idxint idx, char * default_file, FILE ** fb);

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
	idxint i;
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

#ifdef DLONG
    #define INTRW "%ld"
#else
    #define INTRW "%i"
#endif

#ifndef DFLOAT
    #define FLOATRW "%lf"
#else
    #define FLOATRW "%f"
#endif


idxint read_in_data(FILE * fp,Data * d, Cone * k){
	/* MATRIX IN DATA FILE MUST BE IN COLUMN COMPRESSED FORMAT */
	d->RHO_X = RHOX;
    if(fscanf(fp, INTRW, &(d->n)) != 1) return -1;
	if(fscanf(fp, INTRW, &(d->m))!= 1) return -1;
    if(fscanf(fp, INTRW, &(k->f))!= 1) return -1;
	if(fscanf(fp, INTRW, &(k->l))!= 1) return -1;
    if(fscanf(fp, INTRW, &(k->qsize))!= 1) return -1;
	if(fscanf(fp, INTRW, &(k->ssize))!= 1) return -1;

    // allow arbitrary additional cones, simply add to below:
    int len = 32;
    char s[len];
    if( fgets (s, len, fp) == NULL ) return -1;
    char* token = strtok(s, " ");
    if(token) k->ep = atoi(token); token = strtok(NULL, " ");
    if(token) k->ed = atoi(token); token = strtok(NULL, " ");
    
	//if(fscanf(fp, INTRW, &(k->ep))!= 1) return -1;
	//if(fscanf(fp, INTRW, &(k->ed))!= 1) return -1;
	
	if(fscanf(fp, INTRW, &(d->MAX_ITERS))!= 1) return -1;
	if(fscanf(fp, INTRW, &(d->VERBOSE))!= 1) return -1;
	if(fscanf(fp, INTRW, &(d->NORMALIZE))!= 1) return -1;
	if(fscanf(fp, FLOATRW, &(d->ALPHA))!= 1) return -1;
	if(fscanf(fp, FLOATRW, &(d->UNDET_TOL))!= 1) return -1;
	if(fscanf(fp, FLOATRW, &(d->EPS))!= 1) return -1;
	k->q = malloc(sizeof(idxint)*k->qsize);
	for(idxint i = 0; i < k->qsize; i++)
	{ 
		if(fscanf(fp, INTRW, &k->q[i])!= 1) return -1;
	}
    k->s = malloc(sizeof(idxint)*k->ssize);
    for(idxint i = 0; i < k->ssize; i++)
    {   
        if(fscanf(fp, INTRW, &k->s[i])!= 1) return -1;
    }   
	d->b = malloc(sizeof(pfloat)*d->m);
	for(idxint i = 0; i < d->m; i++)
	{ 
		if(fscanf(fp, FLOATRW, &d->b[i])!= 1) return -1;
	}
	d->c = malloc(sizeof(pfloat)*d->n);
	for(idxint i = 0; i < d->n; i++)
	{ 
		if(fscanf(fp, FLOATRW, &d->c[i])!= 1) return -1;
	}
    d->Ap = malloc(sizeof(idxint)*(d->n+1));
	for(idxint i = 0; i < d->n+1; i++) 
	{
		if(fscanf(fp, INTRW, &d->Ap[i])!= 1) return -1;
	}
    idxint Anz = d->Ap[d->n];
	d->Ai = malloc(sizeof(idxint)*Anz);
	for(idxint i = 0; i < Anz; i++)
	{
		if(fscanf(fp, INTRW, &d->Ai[i])!= 1) return -1;
	}
	d->Ax = malloc(sizeof(pfloat)*Anz);
	for(idxint i = 0; i < Anz; i++)
	{
		if(fscanf(fp, FLOATRW, &d->Ax[i])!= 1) return -1;
	}
	//		fscanf(fp, "%zu", &NNZ);
	//		idxint *Kr = malloc(sizeof(idxint)*NNZ);
	//		for(idxint i = 0; i < NNZ; i++)
	//		{
	//		fscanf(fp, INTRW, &Kr[i]);
	//		}
	//		idxint *Kp=malloc(sizeof(idxint)*(w->l+1));
	//		for(idxint i = 0; i < w->l+1; i++)
	//		{
	//		fscanf(fp, INTRW, &Kp[i]);
	//		}
	//		pfloat *Kx=malloc(sizeof(pfloat)*NNZ);
	//		for(idxint i = 0; i < NNZ; i++)
	//		{
	//		fscanf(fp, FLOATRW, &Kx[i]);
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


idxint open_file(idxint argc, char ** argv, idxint idx, char * default_file, FILE ** fb) 
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
