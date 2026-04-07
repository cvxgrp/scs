#include "kernels.h"

static __global__ void scatter_diag_kernel(scs_float *kkt_val,
                                           const scs_int *idxs,
                                           const scs_float *new_vals,
                                           scs_int len) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < len) {
    kkt_val[idxs[i]] = new_vals[i];
  }
}

extern "C" cudaError_t scs_scatter_diag(scs_float *d_kkt_val,
                                        const scs_int *d_idxs,
                                        const scs_float *d_new_vals,
                                        scs_int len, cudaStream_t stream) {
  if (len <= 0) {
    return cudaSuccess;
  }
  int threads = 256;
  int blocks = (len + threads - 1) / threads;
  scatter_diag_kernel<<<blocks, threads, 0, stream>>>(d_kkt_val, d_idxs,
                                                       d_new_vals, len);
  return cudaPeekAtLastError();
}
