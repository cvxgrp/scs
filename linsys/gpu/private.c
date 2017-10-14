#include "private.h"

#ifdef __cplusplus
extern "C" {
#endif

#define CG_BEST_TOL 1e-9
#define CG_MIN_TOL 1e-1

#define CUDA_CHECK_ERR                                                         \
    do {                                                                       \
        cuda_error_t err = cuda_get_last_error();                                  \
        if (err != cuda_success) {                                              \
            printf("%s:%d:%s\n ERROR_CUDA: %s\n", __FILE__, __LINE__,          \
                   __func__, cuda_get_error_string(err));                         \
        }                                                                      \
    } while (0)

#ifndef EXTRAVERBOSE
#ifndef FLOAT
#define CUBLAS(x) cublas_d##x
#define CUSPARSE(x) cusparse_d##x
#else
#define CUBLAS(x) cublas_s##x
#define CUSPARSE(x) cusparse_s##x
#endif
#else
#ifndef FLOAT
#define CUBLAS(x)                                                              \
    CUDA_CHECK_ERR;                                                            \
    cublas_d##x
#define CUSPARSE(x)                                                            \
    CUDA_CHECK_ERR;                                                            \
    cusparse_d##x
#else
#define CUBLAS(x)                                                              \
    CUDA_CHECK_ERR;                                                            \
    cublas_s##x
#define CUSPARSE(x)                                                            \
    CUDA_CHECK_ERR;                                                            \
    cusparse_s##x
#endif
#endif

/*
 CUDA matrix routines only for CSR, not CSC matrices:
    CSC             CSR             GPU     Mult
    A (m x n)       A' (n x m)      Ag      accum_by_a_trans_gpu
    A'(n x m)       A  (m x n)      Agt     accum_by_a_gpu
*/

void accum_by_atrans_gpu(const Priv *p, const scs_float *x, scs_float *y) {
    /* y += A'*x
       x and y MUST be on GPU already
    */
    const scs_float onef = 1.0;
    AMatrix *Ag = p->Ag;
    CUSPARSE(csrmv)(p->cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE, Ag->n,
                    Ag->m, p->Annz, &onef, p->descr, Ag->x, Ag->p, Ag->i, x,
                    &onef, y);
}

void accum_by_a_gpu(const Priv *p, const scs_float *x, scs_float *y) {
    /* y += A*x
       x and y MUST be on GPU already
     */
    const scs_float onef = 1.0;
    AMatrix *Agt = p->Agt;
    CUSPARSE(csrmv)(p->cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE, Agt->n,
                    Agt->m, p->Annz, &onef, p->descr, Agt->x, Agt->p, Agt->i, x,
                    &onef, y);
}

/* do not use within pcg, reuses memory */
void accum_by_atrans(const AMatrix *A, Priv *p, const scs_float *x,
                   scs_float *y) {
    scs_float *v_m = p->tmp_m;
    scs_float *v_n = p->r;
    cuda_memcpy(v_m, x, A->m * sizeof(scs_float), cuda_memcpy_host_to_device);
    cuda_memcpy(v_n, y, A->n * sizeof(scs_float), cuda_memcpy_host_to_device);
    accum_by_atrans_gpu(p, v_m, v_n);
    cuda_memcpy(y, v_n, A->n * sizeof(scs_float), cuda_memcpy_device_to_host);
}

/* do not use within pcg, reuses memory */
void accum_by_a(const AMatrix *A, Priv *p, const scs_float *x, scs_float *y) {
    scs_float *v_m = p->tmp_m;
    scs_float *v_n = p->r;
    cuda_memcpy(v_n, x, A->n * sizeof(scs_float), cuda_memcpy_host_to_device);
    cuda_memcpy(v_m, y, A->m * sizeof(scs_float), cuda_memcpy_host_to_device);
    accum_by_a_gpu(p, v_n, v_m);
    cuda_memcpy(y, v_m, A->m * sizeof(scs_float), cuda_memcpy_device_to_host);
}

char *get_lin_sys_method(const AMatrix *A, const Settings *s) {
    char *str = (char *)scs_malloc(sizeof(char) * 128);
    sprintf(str, "sparse-indirect GPU, nnz in A = %li, CG tol ~ 1/iter^(%2.2f)",
            (long)A->p[A->n], s->cg_rate);
    return str;
}

char *get_lin_sys_summary(Priv *p, const Info *info) {
    char *str = (char *)scs_malloc(sizeof(char) * 128);
    sprintf(str,
            "\tLin-sys: avg # CG iterations: %2.2f, avg solve time: %1.2es\n",
            (scs_float)p->tot_cg_its / (info->iter + 1),
            p->total_solve_time / (info->iter + 1) / 1e3);
    p->tot_cg_its = 0;
    p->total_solve_time = 0;
    return str;
}

void cuda_free_a_matrix(AMatrix *A) {
    if (A->x)
        cuda_free(A->x);
    if (A->i)
        cuda_free(A->i);
    if (A->p)
        cuda_free(A->p);
}

void free_priv(Priv *p) {
    if (p) {
        if (p->p)
            cuda_free(p->p);
        if (p->r)
            cuda_free(p->r);
        if (p->Gp)
            cuda_free(p->Gp);
        if (p->bg)
            cuda_free(p->bg);
        if (p->tmp_m)
            cuda_free(p->tmp_m);
        if (p->z)
            cuda_free(p->z);
        if (p->M)
            cuda_free(p->M);
        if (p->Ag) {
            cuda_free_a_matrix(p->Ag);
            scs_free(p->Ag);
        }
        if (p->Agt) {
            cuda_free_a_matrix(p->Agt);
            scs_free(p->Agt);
        }
        cusparse_destroy(p->cusparse_handle);
        cublas_destroy(p->cublas_handle);
        cuda_device_reset();
        scs_free(p);
    }
}

/*y = (RHO_X * I + A'A)x */
static void mat_vec(const AMatrix *A, const Settings *s, Priv *p,
                   const scs_float *x, scs_float *y) {
    /* x and y MUST already be loaded to GPU */
    scs_float *tmp_m = p->tmp_m; /* temp memory */
    cuda_memset(tmp_m, 0, A->m * sizeof(scs_float));
    accum_by_a_gpu(p, x, tmp_m);
    cuda_memset(y, 0, A->n * sizeof(scs_float));
    accum_by_atrans_gpu(p, tmp_m, y);
    CUBLAS(axpy)(p->cublas_handle, A->n, &(s->rho_x), x, 1, y, 1);
}

/* M = inv ( diag ( RHO_X * I + A'A ) ) */
void get_preconditioner(const AMatrix *A, const Settings *stgs, Priv *p) {
    scs_int i;
    scs_float *M = (scs_float *)scs_malloc(A->n * sizeof(scs_float));

#if EXTRAVERBOSE > 0
    scs_printf("getting pre-conditioner\n");
#endif

    for (i = 0; i < A->n; ++i) {
        M[i] = 1 / (stgs->rho_x +
                    calc_norm_sq(&(A->x[A->p[i]]), A->p[i + 1] - A->p[i]));
        /* M[i] = 1; */
    }
    cuda_memcpy(p->M, M, A->n * sizeof(scs_float), cuda_memcpy_host_to_device);
    scs_free(M);

#if EXTRAVERBOSE > 0
    scs_printf("finished getting pre-conditioner\n");
#endif
}

Priv *init_priv(const AMatrix *A, const Settings *stgs) {
    cuda_error_t err;
    Priv *p = (Priv *)scs_calloc(1, sizeof(Priv));
    p->Annz = A->p[A->n];
    p->cublas_handle = 0;
    p->cusparse_handle = 0;
    p->descr = 0;

    p->total_solve_time = 0;
    p->tot_cg_its = 0;

    /* Get handle to the CUBLAS context */
    cublas_create(&p->cublas_handle);

    /* Get handle to the CUSPARSE context */
    cusparse_create(&p->cusparse_handle);

    /* Matrix description */
    cusparse_create_mat_descr(&p->descr);
    cusparse_set_mat_type(p->descr, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparse_set_mat_index_base(p->descr, CUSPARSE_INDEX_BASE_ZERO);

    AMatrix *Ag = (AMatrix *)scs_malloc(sizeof(AMatrix));
    Ag->n = A->n;
    Ag->m = A->m;
    p->Ag = Ag;

    AMatrix *Agt = (AMatrix *)scs_malloc(sizeof(AMatrix));
    Agt->n = A->m;
    Agt->m = A->n;
    p->Agt = Agt;

    cuda_malloc((void **)&Ag->i, (A->p[A->n]) * sizeof(scs_int));
    cuda_malloc((void **)&Ag->p, (A->n + 1) * sizeof(scs_int));
    cuda_malloc((void **)&Ag->x, (A->p[A->n]) * sizeof(scs_float));

    cuda_malloc((void **)&p->p, A->n * sizeof(scs_float));
    cuda_malloc((void **)&p->r, A->n * sizeof(scs_float));
    cuda_malloc((void **)&p->Gp, A->n * sizeof(scs_float));
    cuda_malloc((void **)&p->bg, (A->n + A->m) * sizeof(scs_float));
    cuda_malloc((void **)&p->tmp_m,
               A->m * sizeof(scs_float)); /* intermediate result */
    cuda_malloc((void **)&p->z, A->n * sizeof(scs_float));
    cuda_malloc((void **)&p->M, A->n * sizeof(scs_float));

    cuda_memcpy(Ag->i, A->i, (A->p[A->n]) * sizeof(scs_int),
               cuda_memcpy_host_to_device);
    cuda_memcpy(Ag->p, A->p, (A->n + 1) * sizeof(scs_int),
               cuda_memcpy_host_to_device);
    cuda_memcpy(Ag->x, A->x, (A->p[A->n]) * sizeof(scs_float),
               cuda_memcpy_host_to_device);

    cuda_malloc((void **)&Agt->i, (A->p[A->n]) * sizeof(scs_int));
    cuda_malloc((void **)&Agt->p, (A->m + 1) * sizeof(scs_int));
    cuda_malloc((void **)&Agt->x, (A->p[A->n]) * sizeof(scs_float));

    get_preconditioner(A, stgs, p);

    /* transpose Ag into Agt for faster multiplies */
    /* TODO: memory intensive, could perform transpose in CPU and copy to GPU */
    CUSPARSE(csr2csc)(p->cusparse_handle, A->n, A->m, A->p[A->n], Ag->x, Ag->p,
                      Ag->i, Agt->x, Agt->i, Agt->p, CUSPARSE_ACTION_NUMERIC,
                      CUSPARSE_INDEX_BASE_ZERO);

    err = cuda_get_last_error();
    if (err != cuda_success) {
        printf("%s:%d:%s\nERROR_CUDA: %s\n", __FILE__, __LINE__, __func__,
               cuda_get_error_string(err));
        free_priv(p);
        return SCS_NULL;
    }
    return p;
}

static void apply_pre_conditioner(cublas_handle_t cublas_handle, scs_float *M,
                                scs_float *z, scs_float *r, scs_int n) {
    cuda_memcpy(z, r, n * sizeof(scs_float), cuda_memcpy_device_to_device);
    CUBLAS(tbmv)(cublas_handle, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N,
                 CUBLAS_DIAG_NON_UNIT, n, 0, M, 1, z, 1);
}

/* solves (I+A'A)x = b, s warm start, solution stored in bg (on GPU) */
static scs_int pcg(const AMatrix *A, const Settings *stgs, Priv *pr,
                   const scs_float *s, scs_float *bg, scs_int max_its,
                   scs_float tol) {
    scs_int i, n = A->n;
    scs_float alpha, nrm_r, p_gp, neg_alpha, beta, ipzr, ipzr_old;
    scs_float onef = 1.0, neg_onef = -1.0;
    scs_float *p = pr->p;   /* cg direction */
    scs_float *Gp = pr->Gp; /* updated CG direction */
    scs_float *r = pr->r;   /* cg residual */
    scs_float *z = pr->z;   /* preconditioned */
    scs_float *M = pr->M;   /* preconditioner */
    cublas_handle_t cublas_handle = pr->cublas_handle;

    if (s == SCS_NULL) {
        cuda_memcpy(r, bg, n * sizeof(scs_float), cuda_memcpy_device_to_device);
        cuda_memset(bg, 0, n * sizeof(scs_float));
    } else {
        /* p contains bg temporarily */
        cuda_memcpy(p, bg, n * sizeof(scs_float), cuda_memcpy_device_to_device);
        /* bg contains s */
        cuda_memcpy(bg, s, n * sizeof(scs_float), cuda_memcpy_host_to_device);
        mat_vec(A, stgs, pr, bg, r);
        CUBLAS(axpy)(cublas_handle, n, &neg_onef, p, 1, r, 1);
        CUBLAS(scal)(cublas_handle, n, &neg_onef, r, 1);
    }

    /* for some reason nrm2 is VERY slow */
    /* CUBLAS(nrm2)(cublas_handle, n, r, 1, &nrm_r); */
    CUBLAS(dot)(cublas_handle, n, r, 1, r, 1, &nrm_r);
    nrm_r = SQRTF(nrm_r);
    /* check to see if we need to run CG at all */
    if (nrm_r < MIN(tol, 1e-18)) {
        return 0;
    }

    apply_pre_conditioner(cublas_handle, M, z, r, n);
    CUBLAS(dot)(cublas_handle, n, r, 1, z, 1, &ipzr);
    /* put z in p, replacing temp mem */
    cuda_memcpy(p, z, n * sizeof(scs_float), cuda_memcpy_device_to_device);

    for (i = 0; i < max_its; ++i) {
        mat_vec(A, stgs, pr, p, Gp);

        CUBLAS(dot)(cublas_handle, n, p, 1, Gp, 1, &p_gp);

        alpha = ipzr / p_gp;
        neg_alpha = -alpha;

        CUBLAS(axpy)(cublas_handle, n, &alpha, p, 1, bg, 1);
        CUBLAS(axpy)(cublas_handle, n, &neg_alpha, Gp, 1, r, 1);

        /* for some reason nrm2 is VERY slow */
        /* CUBLAS(nrm2)(cublas_handle, n, r, 1, &nrm_r); */
        CUBLAS(dot)(cublas_handle, n, r, 1, r, 1, &nrm_r);
        nrm_r = SQRTF(nrm_r);
        if (nrm_r < tol) {
            i++;
            break;
        }
        ipzr_old = ipzr;
        apply_pre_conditioner(cublas_handle, M, z, r, n);
        CUBLAS(dot)(cublas_handle, n, r, 1, z, 1, &ipzr);

        beta = ipzr / ipzr_old;
        CUBLAS(scal)(cublas_handle, n, &beta, p, 1);
        CUBLAS(axpy)(cublas_handle, n, &onef, z, 1, p, 1);
    }
#if EXTRAVERBOSE > 0
    scs_printf("tol: %.4e, resid: %.4e, iters: %li\n", tol, nrm_r, (long)i + 1);
#endif
    return i;
}

#ifdef TEST_GPU_MAT_MUL
void accum_by_atrans_host(const AMatrix *A, Priv *p, const scs_float *x,
                       scs_float *y) {
    _accum_by_atrans(A->n, A->x, A->i, A->p, x, y);
}

void accum_by_a_host(const AMatrix *A, Priv *p, const scs_float *x, scs_float *y) {
    _accum_by_a(A->n, A->x, A->i, A->p, x, y);
}

void test_gpu_mat_mul(const AMatrix *A, Priv *p, scs_float *b) {
    /* test to see if matrix multiplication codes agree */
    scs_float t[A->n + A->m], u[A->n + A->m], *bg;
    cuda_malloc((void **)&bg, (A->n + A->m) * sizeof(scs_float));

    cuda_memcpy(bg, b, (A->n + A->m) * sizeof(scs_float),
               cuda_memcpy_host_to_device);
    memcpy(t, b, (A->n + A->m) * sizeof(scs_float));

    accum_by_atrans_gpu(p, &(bg[A->n]), bg);
    accum_by_atrans_host(A, p, &(t[A->n]), t);
    cuda_memcpy(u, bg, (A->n + A->m) * sizeof(scs_float),
               cuda_memcpy_device_to_host);
    printf("A trans multiplication err %2.e\n", calc_norm_diff(u, t, A->n));

    accum_by_a_gpu(p, bg, &(bg[A->n]));
    accum_by_a_host(A, p, t, &(t[A->n]));
    cuda_memcpy(u, bg, (A->n + A->m) * sizeof(scs_float),
               cuda_memcpy_device_to_host);
    printf("A multiplcation err %2.e\n",
           calc_norm_diff(&(u[A->n]), &(t[A->n]), A->m));
    cuda_free(bg);
}
#endif

scs_int solve_lin_sys(const AMatrix *A, const Settings *stgs, Priv *p,
                    scs_float *b, const scs_float *s, scs_int iter) {
    scs_int cg_its;
    timer linsys_timer;
    scs_float *bg = p->bg;
    scs_float neg_onef = -1.0;
    scs_float cg_tol =
        calc_norm(b, A->n) *
        (iter < 0 ? CG_BEST_TOL
                  : CG_MIN_TOL / POWF((scs_float)iter + 1, stgs->cg_rate));

    tic(&linsys_timer);
/* solves Mx = b, for x but stores result in b */
/* s contains warm-start (if available) */

#ifdef TEST_GPU_MAT_MUL
    test_gpu_mat_mul(A, p, b);
#endif

    /* all on GPU */
    cuda_memcpy(bg, b, (A->n + A->m) * sizeof(scs_float),
               cuda_memcpy_host_to_device);
    accum_by_atrans_gpu(p, &(bg[A->n]), bg);
    /* solves (I+A'A)x = b, s warm start, solution stored in b */
    cg_its = pcg(A, stgs, p, s, bg, A->n, MAX(cg_tol, CG_BEST_TOL));
    CUBLAS(scal)(p->cublas_handle, A->m, &neg_onef, &(bg[A->n]), 1);
    accum_by_a_gpu(p, bg, &(bg[A->n]));
    cuda_memcpy(b, bg, (A->n + A->m) * sizeof(scs_float),
               cuda_memcpy_device_to_host);

    if (iter >= 0) {
        p->tot_cg_its += cg_its;
    }

    p->total_solve_time += tocq(&linsys_timer);
#if EXTRAVERBOSE > 0
    scs_printf("linsys solve time: %1.2es\n", tocq(&linsys_timer) / 1e3);
#endif
    return 0;
}

#ifdef __cplusplus
}
#endif
