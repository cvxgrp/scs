#ifndef CUDSS_KERNELS_H
#define CUDSS_KERNELS_H

#include <cuda_runtime.h>
#include "glbopts.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Scatter n+m new diagonal values into the KKT matrix on the GPU.
 * d_kkt_val: device KKT values array
 * d_idxs: device array of scatter indices (length len)
 * d_new_vals: device array of new diagonal values (length len)
 * len: number of values to scatter (n + m)
 * stream: CUDA stream for async execution
 * Returns cudaSuccess on successful launch, error code otherwise.
 */
cudaError_t scs_scatter_diag(scs_float *d_kkt_val, const scs_int *d_idxs,
                             const scs_float *d_new_vals, scs_int len,
                             cudaStream_t stream);

#ifdef __cplusplus
}
#endif

#endif
