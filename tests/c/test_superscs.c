#include "test_superscs.h"
#include "linsys/amatrix.h"
#include "linsys/common.h"
#include "linsys/direct/external/amd_internal.h"

static void prepare_data(Data ** data) {
    const scs_int n = 3;
    const scs_int m = 4;
    const scs_int nnz = 5;

    AMatrix * A;

    *data = initData();
    (*data)->c = malloc(n * sizeof (scs_float));
    (*data)->c[0] = 1.0;
    (*data)->c[1] = -2.0;
    (*data)->c[2] = -3.0;
    (*data)->b = malloc(m * sizeof (scs_float));
    (*data)->b[0] = 0.2;
    (*data)->b[1] = 0.1;
    (*data)->b[2] = -0.1;
    (*data)->b[3] = 0.1;

    (*data)->m = m;
    (*data)->n = n;

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

    (*data)->A = A;
}

static void prepare_cone(Cone ** cone) {
    *cone = malloc(sizeof (Cone));
    (*cone)->ssize = 0;
    (*cone)->ed = 0;
    (*cone)->ep = 0;
    (*cone)->f = 0;
    (*cone)->l = 0;
    (*cone)->psize = 0;
    (*cone)->ssize = 0;
    (*cone)->qsize = 1;
    (*cone)->q = malloc(4 * sizeof (scs_int));
    (*cone)->q[0] = 4;
    (*cone)->p = SCS_NULL;
    (*cone)->s = SCS_NULL;
}

bool test_superscs_solve(char** str) {

    scs_int status;
    Sol* sol;
    Data * data;
    Info * info;
    Cone * cone;

    prepare_data(&data);
    prepare_cone(&cone);

    data->stgs->eps = 1e-9;
    data->stgs->rho_x = 1.0;
    data->stgs->direction = (direction_type) restarted_broyden;
    data->stgs->verbose = 0;

    sol = initSol();
    info = initInfo();

    data->stgs->do_super_scs = 1;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_SOLVED, str, "Problem not solved");

    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[0], -16.874896969005714, 1e-6, str, "x_star[0] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[1], -5.634341514927034, 1e-6, str, "x_star[1] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[2], 3.589737499286473, 1e-6, str, "x_star[2] wrong");


    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->y[0], 96.506238327408525, 1e-6, str, "y_star[0] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->y[1], -74.161955281143605, 1e-6, str, "y_star[1] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->y[2], 15.000000000002315, 1e-6, str, "y_star[2] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->y[3], 59.903742996445253, 1e-6, str, "y_star[3] wrong");

    ASSERT_EQUAL_FLOAT_OR_FAIL(info->pobj, -16.375426437011065, 1e-7, str, "pobj wrong");

    ASSERT_EQUAL_FLOAT_OR_FAIL(info->pobj, info->dobj, 1e-4, str, "P not equal to D");
    ASSERT_TRUE_OR_FAIL(info->relGap < 1e-10, str, "relative gap too high");
    ASSERT_EQUAL_INT_OR_FAIL(strcmp(info->status, "Solved"), 0, str, "problem not 'Solved'");
    ASSERT_EQUAL_INT_OR_FAIL(info->statusVal, SCS_SOLVED, str, "problem status not SCS_SOLVED");

    freeData(data, cone);
    freeSol(sol);
    freeInfo(info);

    SUCCEED(str);
}

bool test_superscs_000(char** str) {
    scs_int status;
    Sol* sol;
    Data * data;
    Info * info;
    Cone * cone;

    prepare_data(&data);
    prepare_cone(&cone);

    data->stgs->sse = 0.5;
    data->stgs->eps = 1e-4;
    data->stgs->rho_x = 1.0;
    data->stgs->direction = (direction_type) restarted_broyden;
    data->stgs->verbose = 0;
    data->stgs->k0 = 0;
    data->stgs->k1 = 0;
    data->stgs->k2 = 0;
    data->stgs->ls = 0;
    data->stgs->max_iters = 1;
    data->stgs->do_super_scs = 1;

    sol = initSol();
    info = initInfo();

    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(info->iter, data->stgs->max_iters, str, "no iterations");
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_SOLVED_INACCURATE, str, "wrong status");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[0], -0.109053491087962, 1e-10, str, "x[0] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[1], -0.003683620781908, 1e-10, str, "x[1] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[2], 2.645438455229390, 1e-10, str, "x[2] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->y[0], 19.912548935347708, 1e-10, str, "y[0] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->y[1], 2.294092141293139, 1e-10, str, "y[1] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->y[2], 19.262127935028715, 1e-10, str, "y[2] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->y[3], 4.496351161159519, 1e-10, str, "y[3] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->s[0], 0.200927752555239, 1e-10, str, "s[0] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->s[1], -0.023148557203865, 1e-10, str, "s[1] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->s[2], -0.194364673652928, 1e-10, str, "s[2] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->s[3], -0.045370471477269, 1e-10, str, "s[3] wrong");

    data->stgs->max_iters = 2;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(info->iter, data->stgs->max_iters, str, "no iterations");
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_SOLVED_INACCURATE, str, "wrong status");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[0], -1.261463537904218, 1e-10, str, "x[0] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->y[2], 92.364496490679642, 1e-10, str, "y[2] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->s[3], -0.209261231033243, 1e-10, str, "s[3] wrong");

    data->stgs->max_iters = 10;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(info->iter, data->stgs->max_iters, str, "no iterations");
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_SOLVED_INACCURATE, str, "wrong status");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[0], -2.232691713491094, 1e-10, str, "x[0] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[1], -0.590885686812609, 1e-10, str, "x[1] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->y[1], -21.635616132075267, 1e-10, str, "y[1] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->y[3], 20.962816440385115, 1e-10, str, "y[3] wrong");

    data->stgs->max_iters = 1000;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(info->iter, 121, str, "Erroneous no. iter.");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[0], -16.871782513122774, 1e-10, str, "x[0] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[1], -5.633253618312680, 1e-10, str, "x[1] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[2], 3.589570393256817, 1e-10, str, "x[2] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->y[0], 96.497976811110945, 1e-10, str, "y[0] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->y[1], -74.155009946640732, 1e-10, str, "y[1] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->y[2], 15.000840046600869, 1e-10, str, "y[2] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->y[3], 59.898821577286611, 1e-10, str, "y[3] wrong");
    ASSERT_EQUAL_INT_OR_FAIL(info->statusVal, SCS_SOLVED, str, "problem status not SCS_SOLVED");
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_SOLVED, str, "wrong status");

    freeData(data, cone);
    freeSol(sol);
    freeInfo(info);

    SUCCEED(str);
}

bool test_superscs_001_fpr(char** str) {
    scs_int status;
    Sol* sol;
    Data * data;
    Info * info;
    Cone * cone;
    scs_int i;

    prepare_data(&data);
    prepare_cone(&cone);

    data->stgs->sse = 0.5;
    data->stgs->eps = 1e-4;
    data->stgs->rho_x = 1.0;
    data->stgs->direction = (direction_type) fixed_point_residual;
    data->stgs->verbose = 0;
    data->stgs->k0 = 0;
    data->stgs->k1 = 0;
    data->stgs->k2 = 1;
    data->stgs->ls = 10;
    data->stgs->max_iters = 1;
    data->stgs->do_super_scs = 1;

    sol = initSol();
    info = initInfo();

    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_UNBOUNDED_INACCURATE, str, "wrong status");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[0], 0.274057420504456, 1e-10, str, "x[0] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[1], -0.058098186140208, 1e-10, str, "x[1] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[2], 0.463417930928291, 1e-10, str, "x[2] wrong");
    for (i = 0; i < 4; ++i) {
        ASSERT_TRUE_OR_FAIL(isnan(sol->y[i]), str, "y should be nan");
    }
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->s[0], -0.191928792495329, 1e-10, str, "s[0] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->s[1], -0.047508022860835, 1e-10, str, "s[1] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->s[2], 0.182152982731530, 1e-10, str, "s[2] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->s[3], 0.037415772537480, 1e-10, str, "s[3] wrong");

    data->stgs->max_iters = 2;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(info->iter, data->stgs->max_iters, str, "no iterations");
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_SOLVED_INACCURATE, str, "wrong status");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[0], -1.052622260714879, 1e-10, str, "x[0] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->y[0], 23.429483256003490, 1e-10, str, "y[0] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->s[0], 0.402310583070969, 1e-10, str, "s[0] wrong");

    data->stgs->max_iters = 3;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(info->iter, data->stgs->max_iters, str, "no iterations");
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_SOLVED_INACCURATE, str, "wrong status");
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_SOLVED_INACCURATE, str, "wrong status");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[1], -0.3708797830183899, 1e-10, str, "x[1] wrong");

    data->stgs->max_iters = 10;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(info->iter, data->stgs->max_iters, str, "no iterations");
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_SOLVED_INACCURATE, str, "wrong status");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->y[0], 34.954357266943035, 1e-10, str, "y[0] wrong");

    data->stgs->max_iters = 80;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(info->statusVal, SCS_SOLVED, str, "problem status not SCS_SOLVED");
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_SOLVED, str, "wrong status");

    freeData(data, cone);
    freeSol(sol);
    scs_free(info);

    SUCCEED(str);
}

bool test_superscs_001_rbroyden(char** str) {
    scs_int status;
    Sol* sol;
    Data * data;
    Info * info;
    Cone * cone;
    scs_int i;

    prepare_data(&data);
    prepare_cone(&cone);

    data->stgs->sse = 0.5;
    data->stgs->eps = 1e-4;
    data->stgs->rho_x = 1.;
    data->stgs->direction = (direction_type) restarted_broyden;
    data->stgs->verbose = 0;
    data->stgs->k0 = 0;
    data->stgs->k1 = 0;
    data->stgs->k2 = 1;
    data->stgs->ls = 10;
    data->stgs->max_iters = 1;
    data->stgs->do_super_scs = 1;
    data->stgs->memory = 10;
    data->stgs->sigma = 1e-2;
    data->stgs->c1 = 1.0 - 1e-4;
    data->stgs->c_bl = 0.999;
    data->stgs->beta = 0.5;
    data->stgs->normalize = 1;
    data->stgs->scale = 1;
    data->stgs->alpha = 1.5;

    sol = initSol();
    info = initInfo();


    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_UNBOUNDED_INACCURATE, str, "wrong status");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[0], 0.274057420504456, 1e-10, str, "x[0] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[1], -0.058098186140208, 1e-10, str, "x[1] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[2], 0.463417930928291, 1e-10, str, "x[2] wrong");
    for (i = 0; i < 4; ++i) {
        ASSERT_TRUE_OR_FAIL(isnan(sol->y[i]), str, "y should be nan");
    }
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->s[0], -0.191928792495329, 1e-10, str, "s[0] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->s[1], -0.047508022860835, 1e-10, str, "s[1] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->s[2], 0.182152982731530, 1e-10, str, "s[2] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->s[3], 0.037415772537480, 1e-10, str, "s[3] wrong");

    data->stgs->max_iters = 2;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(info->iter, data->stgs->max_iters, str, "no iterations");
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_SOLVED_INACCURATE, str, "wrong status");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[0], -0.465406066728364, 1e-10, str, "x[0] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[1], -0.166978364590537, 1e-10, str, "x[1] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[2], 1.116606860418411, 1e-10, str, "x[2] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->y[0], 7.224785302606174, 1e-10, str, "y[0] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->s[0], 0.326281791938008, 1e-10, str, "s[0] wrong");


    data->stgs->max_iters = 11;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(info->iter, data->stgs->max_iters, str, "no iterations");
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_UNBOUNDED_INACCURATE, str, "wrong status");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[0], -1.046552668150064, 1e-10, str, "x[0] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[1], -0.353299417556677, 1e-10, str, "x[1] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[2], 0.220015388987763, 1e-10, str, "x[2] wrong");
    for (i = 0; i < 4; ++i) {
        ASSERT_TRUE_OR_FAIL(isnan(sol->y[i]), str, "y should be nan");
    }
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->s[0], 0.380569177488686, 1e-10, str, "s[0] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->s[1], 0.283386052682735, 1e-10, str, "s[1] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->s[2], -0.094498466741551, 1e-10, str, "s[2] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->s[3], -0.235786521630921, 1e-10, str, "s[3] wrong");

    data->stgs->max_iters = 40;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(info->iter, data->stgs->max_iters, str, "no iterations");
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_SOLVED_INACCURATE, str, "wrong status");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[0], -18.660744885301725, 1e-10, str, "x[0] wrong");

    /*
     * Here I'm modifying the maximum number of iterations to make sure that  
     * those tricks with stgs->previous_max_iter indeed work.
     */

    data->stgs->max_iters = 1000;
    data->stgs->eps = 1e-4;
    data->stgs->rho_x = 0.5;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(info->statusVal, SCS_SOLVED, str, "problem status not SCS_SOLVED");
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_SOLVED, str, "wrong status");
    ASSERT_TRUE_OR_FAIL(info->progress_dcost == SCS_NULL, str, "progress not NULL");
    ASSERT_TRUE_OR_FAIL(info->progress_pcost == SCS_NULL, str, "progress not NULL");
    ASSERT_TRUE_OR_FAIL(info->progress_relgap == SCS_NULL, str, "progress not NULL");
    ASSERT_TRUE_OR_FAIL(info->progress_respri == SCS_NULL, str, "progress not NULL");
    ASSERT_TRUE_OR_FAIL(info->progress_resdual == SCS_NULL, str, "progress not NULL");

    data->stgs->max_iters = 2000;
    data->stgs->do_record_progress = 1;
    data->stgs->rho_x = .1;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(info->statusVal, SCS_SOLVED, str, "problem status not SCS_SOLVED");
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_SOLVED, str, "wrong status");
    ASSERT_TRUE_OR_FAIL(info->progress_dcost != SCS_NULL, str, "progress NULL");
    ASSERT_TRUE_OR_FAIL(info->progress_pcost != SCS_NULL, str, "progress NULL");
    ASSERT_TRUE_OR_FAIL(info->progress_relgap != SCS_NULL, str, "progress NULL");
    ASSERT_TRUE_OR_FAIL(info->progress_respri != SCS_NULL, str, "progress NULL");
    ASSERT_TRUE_OR_FAIL(info->progress_resdual != SCS_NULL, str, "progress NULL");
    ASSERT_EQUAL_INT_OR_FAIL(data->stgs->max_iters, 2000, str, "Wrong previous no. iter");
    ASSERT_EQUAL_INT_OR_FAIL(data->stgs->previous_max_iters, 2000, str, "Wrong previous no. iter");

    data->stgs->max_iters = 3000;
    data->stgs->rho_x = .01;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(info->statusVal, SCS_SOLVED, str, "problem status not SCS_SOLVED");
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_SOLVED, str, "wrong status");

    data->stgs->max_iters = 2000;
    data->stgs->rho_x = .001;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(info->statusVal, SCS_SOLVED, str, "problem status not SCS_SOLVED");
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_SOLVED, str, "wrong status");

    data->stgs->max_iters = 3100;
    data->stgs->rho_x = .0001;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(info->statusVal, SCS_SOLVED, str, "problem status not SCS_SOLVED");
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_SOLVED, str, "wrong status");

    data->stgs->max_iters = 3200;
    data->stgs->rho_x = 10;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(info->statusVal, SCS_SOLVED, str, "problem status not SCS_SOLVED");
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_SOLVED, str, "wrong status");

    data->stgs->max_iters = 3300;
    data->stgs->rho_x = 0.001;
    data->stgs->do_super_scs = 0;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(info->statusVal, SCS_SOLVED, str, "problem status not SCS_SOLVED");
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_SOLVED, str, "wrong status");

    freeData(data, cone);
    freeSol(sol);
    freeInfo(info);

    SUCCEED(str);
}

bool test_superscs_100_rbroyden(char** str) {

    scs_int status;
    Sol* sol;
    Data * data;
    Info * info;
    Cone * cone;

    prepare_data(&data);
    prepare_cone(&cone);

    data->stgs->eps = 1e-4;
    data->stgs->rho_x = 1.0;
    data->stgs->direction = (direction_type) restarted_broyden;
    data->stgs->verbose = 0;
    data->stgs->k0 = 1;
    data->stgs->k1 = 0;
    data->stgs->k2 = 0;
    data->stgs->ls = 10;

    data->stgs->verbose = 0;
    data->stgs->do_super_scs = 1;

    sol = initSol();
    info = initInfo();

    data->stgs->max_iters = 2;
    status = scs(data, cone, sol, info);

    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_SOLVED_INACCURATE, str, "wrong status");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[0], -0.349018320302040, 1e-10, str, "x[0] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[1], 0.015102755569314, 1e-10, str, "x[1] wrong");
    ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[2], 1.778110351429428, 1e-10, str, "x[2] wrong");


    freeData(data, cone);
    freeSol(sol);
    freeInfo(info);

    SUCCEED(str);
}

bool test_superscs_011_progress(char** str) {
    scs_int status;
    Sol* sol;
    Data * data;
    Info * info = SCS_NULL;
    Cone * cone;
    scs_int mode_exp[20] = {1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2};
    scs_int ls_exp[20] = {1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4};
    scs_float gap_exp[6] = {0.521136588669913, 0.452292411994033, 0.390774095584746, 0.084086837089436, 0.010215190495717, 0.737805654686660};
    scs_float pres_exp[6] = {5.651585413034721, 0.881043969782380, 0.271916642815964, 0.052153221225419, 0.045003598675170, 0.480719322491409};
    scs_float dres_exp[6] = {9.739130818329103, 5.430982744844291, 0.738951762380322, 0.262551525410492, 0.246763246541060, 0.634137521605225};
    scs_int i;

    prepare_data(&data);
    prepare_cone(&cone);

    data->stgs->eps = 1e-8;
    data->stgs->do_super_scs = 1;
    data->stgs->verbose = 0;
    data->stgs->memory = 10;
    data->stgs->do_record_progress = 1;

    info = initInfo();
    sol = initSol();

    status = scs(data, cone, sol, info);

    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_SOLVED, str, "wrong status");
    for (i = 0; i < 20; ++i) {
        ASSERT_EQUAL_INT_OR_FAIL(info->progress_mode[i], mode_exp[i], str, "mode");
        ASSERT_EQUAL_INT_OR_FAIL(info->progress_ls[i], ls_exp[i], str, "ls");
    }
    for (i = 0; i < 6; ++i) {
        ASSERT_EQUAL_FLOAT_OR_FAIL(info->progress_relgap[i], gap_exp[i], 1e-10, str, "gap");
        ASSERT_EQUAL_FLOAT_OR_FAIL(info->progress_respri[i], pres_exp[i], 1e-10, str, "pres");
        ASSERT_EQUAL_FLOAT_OR_FAIL(info->progress_resdual[i], dres_exp[i], 1e-10, str, "dres");
    }

    freeData(data, cone);
    freeSol(sol);
    freeInfo(info);

    SUCCEED(str);
}

bool test_residuals(char** str) {
    scs_int status;
    Sol* sol;
    Data * data;
    Info * info;
    Cone * cone;

    scs_float relgap_expected[12] = {
        0,
        0.641360567339623,
        0.258326003751872,
        0.427755914304124,
        0.070601296495286,
        0.136391692925419,
        0.110228818969576,
        0.116212468002787,
        0.100073649960616,
        0.037913746742520,
        0.031013566758557,
        0.031786667245133,
    };

    prepare_data(&data);
    prepare_cone(&cone);

    data->stgs->eps = 1e-8;
    data->stgs->k0 = 0;
    data->stgs->k1 = 1;
    data->stgs->k2 = 1;
    data->stgs->ls = 10;
    data->stgs->rho_x = 1.0;
    data->stgs->direction = 100;
    data->stgs->sse = 0.999;
    data->stgs->sigma = 1e-2;
    data->stgs->c_bl = 0.999;
    data->stgs->c1 = 1.0 - 1e-4;
    data->stgs->beta = 0.5;
    data->stgs->normalize = 1;
    data->stgs->scale = 1;
    data->stgs->alpha = 1.5;
    data->stgs->do_super_scs = 1;
    data->stgs->verbose = 1;
    data->stgs->do_override_streams = 1;
    data->stgs->output_stream = stderr;
    data->stgs->do_record_progress = 1;
    data->stgs->max_iters = 120;

    sol = initSol();
    info = initInfo();

    status = scs(data, cone, sol, info);
    ASSERT_TRUE_OR_FAIL(isnan(info->progress_relgap[0]), str, "rel gap [0] not NAN");
    ASSERT_EQUAL_ARRAY_OR_FAIL(info->progress_relgap + 1, relgap_expected + 1, 11, 1e-13, str, "relative gap");
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_SOLVED, str, "wrong status");

    freeData(data, cone);
    freeSol(sol);
    freeInfo(info);

    SUCCEED(str);
}

bool test_rho_x(char** str) {
    scs_int status;
    Sol* sol;
    Data * data;
    Info * info;
    Cone * cone;
    scs_int i, j, k, l;
    scs_float rho;

    prepare_data(&data);
    prepare_cone(&cone);

    data->stgs->eps = 1e-8;
    data->stgs->do_super_scs = 1;
    data->stgs->alpha = 1.5;
    data->stgs->scale = 1.0;
    data->stgs->verbose = 0;
    data->stgs->normalize = 1;
    data->stgs->direction = restarted_broyden;
    data->stgs->beta = 0.5;
    data->stgs->c1 = 0.9999;
    data->stgs->c_bl = 0.999;

    data->stgs->sigma = 1e-2;
    data->stgs->thetabar = 0.1;
    data->stgs->memory = 10;
    data->stgs->sse = 0.999;
    data->stgs->do_record_progress = 1;
    data->stgs->max_iters = 2000;


    sol = initSol();
    info = initInfo();

    /* Make sure it converges for different values of rho_x */
    data->stgs->rho_x = 0.2;
    for (i = 0; i <= 1; ++i) {
        for (j = 0; j <= 1; ++j) {
            for (k = 0; k <= 1; ++k) {
                for (l = 0; l <= 5; ++l) {
                    for (rho = 0.001; rho < 0.5; rho *= 1.8) {
                        /* Test for all possible choices of k0, k1, k2 */
                        data->stgs->k0 = i;
                        data->stgs->k1 = j;
                        data->stgs->k2 = k;
                        /* Test for different values of ls */
                        data->stgs->ls = l;
                        /* Test for different values of rho_x */
                        data->stgs->rho_x = rho;
                        status = scs(data, cone, sol, info);
                        ASSERT_EQUAL_INT_OR_FAIL(status, SCS_SOLVED, str, "wrong status");
                        /* verify that the residuals are lower than eps */
                        ASSERT_TRUE_OR_FAIL(info->resPri < data->stgs->eps, str, "primal residual too high");
                        ASSERT_TRUE_OR_FAIL(info->resDual < data->stgs->eps, str, "dual residual too high");
                        ASSERT_TRUE_OR_FAIL(info->relGap < data->stgs->eps, str, "duality gap too high");
                        ASSERT_TRUE_OR_FAIL(info->iter < data->stgs->max_iters, str, "too many iterations");
                        /* verify the solution */
                        ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[0], -16.874896967418358, 1e-5, str, "x[0] wrong");
                        ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[1], -5.634341514381794, 1e-5, str, "x[1] wrong");
                        ASSERT_EQUAL_FLOAT_OR_FAIL(sol->x[2], 3.589737499067709, 1e-5, str, "x[2] wrong");
                        ASSERT_EQUAL_FLOAT_OR_FAIL(sol->y[0], 96.506238321412667, 1e-5, str, "y[0] wrong");
                        ASSERT_EQUAL_FLOAT_OR_FAIL(sol->y[1], -74.161955276422589, 1e-5, str, "y[1] wrong");
                        ASSERT_EQUAL_FLOAT_OR_FAIL(sol->y[2], 14.999999999635596, 1e-5, str, "y[2] wrong");
                        ASSERT_EQUAL_FLOAT_OR_FAIL(sol->y[3], 59.903742992722336, 1e-5, str, "y[3] wrong");
                        ASSERT_EQUAL_FLOAT_OR_FAIL(sol->s[0], 5.262469090219521, 1e-5, str, "s[0] wrong");
                        ASSERT_EQUAL_FLOAT_OR_FAIL(sol->s[1], 4.044039060072064, 1e-5, str, "s[1] wrong");
                        ASSERT_EQUAL_FLOAT_OR_FAIL(sol->s[2], -0.817947499813188, 1e-5, str, "s[2] wrong");
                        ASSERT_EQUAL_FLOAT_OR_FAIL(sol->s[3], -3.266541120769288, 1e-5, str, "s[3] wrong");
                    }
                }
            }
        }
    }

    freeData(data, cone);
    freeSol(sol);
    freeInfo(info);

    SUCCEED(str);
}

bool test_validation(char** str) {
    scs_int status;
    Sol* sol;
    Data * data;
    Info * info;
    Cone * cone;

    prepare_data(&data);
    prepare_cone(&cone);
    sol = initSol();
    info = initInfo();

    data->stgs->do_override_streams = 1;
    data->stgs->k0 = 2;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_FAILED, str, "k0");
    setDefaultSettings(data);

    data->stgs->do_override_streams = 1;
    data->stgs->k1 = -1;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_FAILED, str, "k1");
    setDefaultSettings(data);

    data->stgs->do_override_streams = 1;
    data->stgs->k2 = 5;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_FAILED, str, "k2");
    setDefaultSettings(data);

    data->stgs->do_override_streams = 1;
    data->stgs->alpha = -0.1;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_FAILED, str, "alpha < 0");
    setDefaultSettings(data);

    data->stgs->do_override_streams = 1;
    data->stgs->alpha = 2.1;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_FAILED, str, "alpha > 2");
    setDefaultSettings(data);

    data->stgs->do_override_streams = 1;
    data->stgs->beta = 1.01;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_FAILED, str, "beta>1");
    setDefaultSettings(data);

    data->stgs->do_override_streams = 1;
    data->stgs->beta = 1;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_FAILED, str, "beta=1");
    setDefaultSettings(data);

    data->stgs->do_override_streams = 1;
    data->stgs->ls = 40;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_FAILED, str, "ls=40");
    setDefaultSettings(data);

    data->stgs->do_override_streams = 1;
    data->stgs->ls = -1;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_FAILED, str, "ls=-1");
    setDefaultSettings(data);

    data->stgs->do_override_streams = 1;
    data->stgs->sigma = -0.0001;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_FAILED, str, "sigma < 0");
    setDefaultSettings(data);

    data->stgs->do_override_streams = 1;
    data->stgs->c_bl = -1e-4;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_FAILED, str, "c_bl");
    setDefaultSettings(data);

    data->stgs->do_override_streams = 1;
    data->stgs->c1 = -1e-4;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_FAILED, str, "c1");
    setDefaultSettings(data);

    freeData(data, cone);
    freeSol(sol);
    freeInfo(info);

    SUCCEED(str);
}

bool test_no_normalization(char** str) {
    scs_int status;
    Sol* sol;
    Data * data;
    Info * info;
    Cone * cone;
    Settings *s;

    prepare_data(&data);
    prepare_cone(&cone);
    sol = initSol();
    info = initInfo();

    s = data->stgs;

    s->normalize = 0;
    s->k0 = 1;
    s->do_super_scs = 1;
    s->eps = 1e-10;
    s->max_iters = 120;
    s->do_override_streams = 1;
    s->verbose = 1;
    s->output_stream = stderr;
    status = scs(data, cone, sol, info);
    ASSERT_EQUAL_INT_OR_FAIL(status, SCS_SOLVED, str, "wrong status");
    ASSERT_TRUE_OR_FAIL(info->resPri < data->stgs->eps, str, "primal residual too high");
    ASSERT_TRUE_OR_FAIL(info->resDual < data->stgs->eps, str, "dual residual too high");
    ASSERT_TRUE_OR_FAIL(info->relGap < data->stgs->eps, str, "duality gap too high");
    ASSERT_TRUE_OR_FAIL(info->iter < data->stgs->max_iters, str, "too many iterations");

    freeData(data, cone);
    freeSol(sol);
    freeInfo(info);

    SUCCEED(str);
}
