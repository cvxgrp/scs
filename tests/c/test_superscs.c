#include "test_superscs.h"
#include "linsys/amatrix.h"

bool test_superscs(char** str) {

    const scs_int n = 3;
    const scs_int m = 4;
    const scs_int nnz = 5;
    Sol* sol;
    Data * data;
    AMatrix * A;
    Info info = {0};

    data = malloc(sizeof (Data));

    data->c = malloc(n * sizeof (scs_float));
    data->c[0] = 1.0;
    data->c[1] = -2.0;
    data->c[2] = -3.0;

    data->b = malloc(m * sizeof (scs_float));
    data->b[0] = 0.2;
    data->b[1] = 0.1;
    data->b[2] = -0.1;
    data->b[3] = 0.1;

    data->m = m;
    data->n = n;


    A = malloc(sizeof (AMatrix));
    A->m = m;
    A->n = n;
    A->p = malloc((n + 1) * sizeof (scs_int));
    A->i = malloc(nnz * sizeof (scs_int));
    A->x = malloc(nnz * sizeof (scs_float));
    
    A->p[0] = 0;
    A->p[1] = 2;
    A->p[2] = 4;
    A->p[3] = 5;

    A->i[0] = 0;
    A->i[1] = 3;
    A->i[2] = 1;
    A->i[3] = 3;
    A->i[4] = 2;

    A->x[0] = 0.3;
    A->x[1] = -0.5;
    A->x[2] = 0.7;
    A->x[3] = 0.9;
    A->x[4] = 0.2;
       
    printAMatrix(A);
    
    data->A = A;
    data->stgs = scs_malloc(sizeof(Settings));
    data->stgs->max_iters = 1000;
    
    data->stgs->alpha = 1.5;
    data->stgs->beta = 0.5;
    data->stgs->c1 = C1;
    data->stgs->c_bl = 0.999;
    data->stgs->eps = .00100;
    data->stgs->k0 = 1;
    data->stgs->k1 = 1;
    data->stgs->k2 = 1;
    data->stgs->ls = 10;
    data->stgs->normalize = 0;
    data->stgs->warm_start = 0;    
    data->stgs->rho_x = 1;
    data->stgs->scale = 1;    
    data->stgs->verbose = 1;                    
    data->stgs->sigma = 1e-2;
    data->stgs->thetabar = 0.1;
    data->stgs->sse = 1e-3;
    data->stgs->memory = 10;

    Cone * cone = malloc(sizeof (Cone));
    cone->qsize = 1;
    cone->q = malloc(4*sizeof(scs_int));
    cone->q[0] = 4;
    initCone(cone);

    sol = scs_calloc(1, sizeof (Sol));
    data->stgs->do_super_scs = 1;
    scs(data, cone, sol, &info);

    SUCCEED(str);
}