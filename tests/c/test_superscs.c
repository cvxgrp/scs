#include "test_superscs.h"
#include "linsys/amatrix.h"

bool test_superscs(char** str) {
    
    scs_int status;
    const scs_int n = 3;
    const scs_int m = 4;
    const scs_int nnz = 5;
    Sol* sol;
    Data * data;
    AMatrix * A;
    Info * info;
    Cone * cone;

    data = initData();

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


    data->A = A;
    data->stgs->eps = 1e-9;
    data->stgs->rho_x = 1.0;   
    data->stgs->sse = 0.7;
    data->stgs->direction = (direction_type)restarted_broyden;
    data->stgs->k0 = 1;
    data->stgs->k1 = 1;
    data->stgs->k2 = 1;
    data->stgs->ls = 1000;
    data->stgs->verbose = 2;
    data->stgs->scale=1;
    data->stgs->normalize = 1;
    
    cone = malloc(sizeof (Cone));
    cone->ssize = 0;
    cone->ed = 0;
    cone->ep = 0;
    cone->f = 0;
    cone->l = 0;
    cone->psize = 0;
    cone->ssize = 0;
    cone->qsize = 1;
    cone->q = malloc(4 * sizeof (scs_int));
    cone->q[0] = 4;

    cone->p = SCS_NULL;
    cone->s = SCS_NULL;

    sol = initSol();
    info = initInfo();

    data->stgs->do_super_scs = 1;
    status = scs(data, cone, sol, info);
    ASSERT_EQAUL_INT_OR_FAIL(status, SCS_SOLVED, str, "Problem not solved");

    ASSERT_EQAUL_FLOAT_OR_FAIL(sol->x[0], -16.874896969005714, 1e-6, str, "x_star[0] wrong");
    ASSERT_EQAUL_FLOAT_OR_FAIL(sol->x[1], -5.634341514927034, 1e-6, str, "x_star[1] wrong");
    ASSERT_EQAUL_FLOAT_OR_FAIL(sol->x[2], 3.589737499286473, 1e-6, str, "x_star[2] wrong");
            
  
    ASSERT_EQAUL_FLOAT_OR_FAIL(sol->y[0], 96.506238327408525, 1e-6, str, "y_star[0] wrong");
    ASSERT_EQAUL_FLOAT_OR_FAIL(sol->y[1], -74.161955281143605, 1e-6, str, "y_star[0] wrong");
    ASSERT_EQAUL_FLOAT_OR_FAIL(sol->y[2], 15.000000000002315, 1e-6, str, "y_star[0] wrong");
    ASSERT_EQAUL_FLOAT_OR_FAIL(sol->y[3], 59.903742996445253, 1e-6, str, "y_star[0] wrong");
    
    ASSERT_EQAUL_FLOAT_OR_FAIL(info->pobj, -16.375426437011065, 1e-7, str, "pobj wrong");
    
    ASSERT_EQAUL_FLOAT_OR_FAIL(info->pobj, info->dobj, 1e-4, str, "P not equal to D");
    ASSERT_TRUE_OR_FAIL(info->relGap<1e-10, str, "relative gap too high");
    ASSERT_EQAUL_INT_OR_FAIL(strcmp(info->status, "Solved"), 0, str, "problem not 'Solved'");
    ASSERT_EQAUL_INT_OR_FAIL(info->statusVal, SCS_SOLVED, str, "problem status not SCS_SOLVED");    
        
    freeData(data, cone);
    freeSol(sol);
    scs_free(info);

    SUCCEED(str);
}