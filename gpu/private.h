#ifndef PRIV_H_GUARD
#define PRIV_H_GUARD

#include "glbopts.h"
#include "scs.h"
#include "cs.h"
#include <math.h>
#include "common.h"
#include "linAlg.h"

#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>

#ifdef XXXXXXX
#define CUDA_CHECK_ERR() \
    do { \
        cudaError_t err = cudaGetLastError(); \
        if (err != cudaSuccess) { \
            printf("%s:%d:%s\n ERROR_CUDA: %s\n", __FILE__, __LINE__, __func__, \
                    cudaGetErrorString(err)); \
            exit(EXIT_FAILURE); \
        } \
    } while (0)
#else
#define CUDA_CHECK_ERR()
#endif

struct PRIVATE_DATA {
    scs_float * p; /* cg iterate, n  */
	scs_float * r; /* cg residual, n */
	scs_float * Gp; /* G * p, n */
    scs_float * bg; /* b, n */
	scs_float * tmp_m; /* m, used in matVec */
	scs_float * tmp_n; /* n */
    scs_int Annz;
    AMatrix * Ag;   /* A matrix on GPU */
	AMatrix * Agt; /* A trans matrix on GPU */

    /* CUDA */
    // XXX: needed?
    cublasHandle_t cublasHandle;
    cublasStatus_t cublasStatus;
    cusparseHandle_t cusparseHandle; 
    cusparseStatus_t cusparseStatus;
    cusparseMatDescr_t descr;
};

#endif
