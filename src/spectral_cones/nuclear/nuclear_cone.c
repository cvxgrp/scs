#include "cones.h"
// #include "scs.h"
#include "linalg.h"
#include "scs_blas.h"
#include "scs_types.h"
#include "util.h" // just for timer

/*
 * Spectral matrix cone projections, from "Projection onto Spectral Matrix
 * Cones" by Daniel Cederberg and Stephen Boyd, 2024.
 *
 * If you have any questions on the code, please reach out to the code author
 * Daniel Cederberg.
 * 
 * This file implements code for projecting onto the nuclear norm cone.
 * 
 * Last modified: 25 August 2024.
 */

void BLAS(gesvd)(const char *jobu, const char *jobvt, const blas_int *m,
                 const blas_int *n, scs_float *a, const blas_int *lda,
                 scs_float *s, scs_float *u, const blas_int *ldu, scs_float *vt,
                 const blas_int *ldvt, scs_float *work, const blas_int *lwork,
                 blas_int *info);

// forward declaration from ell1_cone.c
scs_int ell1_cone_proj_sorted(scs_float t0, const scs_float *x0, scs_float *proj,
                              scs_int n);

// X is of size m x n, stored column major. It is assumed that m >= n.
scs_int SCS(proj_nuclear_cone)(scs_float *tX, size_t m, size_t n, ScsConeWork *c)
{
    assert(m >= n);
    scs_float *X = tX + 1;
    blas_int bm = m;
    blas_int bn = n;

    // -------------------------------------------------------------------------
    //                            Compute SVD
    // -------------------------------------------------------------------------
    scs_float *s = c->s_nuc;
    scs_float *u = c->u_nuc;
    scs_float *vt = c->vt_nuc;

    scs_float *work = c->work_nuc;
    int lwork = c->lwork_nuc;
    int info = 0;

    dgesvd_("S", "A", &bm, &bn, X, &bm, s, u, &bm, vt, &bn, work,
            &lwork, &info);
    if (info != 0)
    {
        printf("WARN: LAPACK gesvd error, info = %i\n", (int)info);
        if (info < 0)
        {
            return info;
        }
    }
    // -------------------------------------------------------------------------
    //                  Project onto spectral *vector* cone
    // -------------------------------------------------------------------------
    scs_int status = ell1_cone_proj_sorted(tX[0], s, tX, n);
    
    if (status < 0)
    {
        return status;
    }

    // -------------------------------------------------------------------------
    //  Recover projection onto spectral *matrix* cone
    // -------------------------------------------------------------------------
    int one = 1;
    for (scs_int i = 0; i < n; ++i)
    {
        BLAS(scal)
        (&bm, &tX[i + 1], &u[i * m], &one);
    }

    char trans = 'N';
    scs_float alpha = 1.0;
    scs_float beta = 0.0;
    BLAS(gemm)
    (&trans, &trans, &bm, &bn, &bn, &alpha, u, &bm, vt, &bn,
     &beta, tX + 1, &bm);

    return 0;
}
