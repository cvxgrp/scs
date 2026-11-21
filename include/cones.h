#ifndef CONES_H_GUARD
#define CONES_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#include "glbopts.h"
#include "scs_blas.h"
#include "scs_work.h"

#ifdef USE_SPECTRAL_CONES
#include "util_spectral_cones.h" /* for Newton_stats */

/* macro for time measurements of SpectralSCS */
#ifdef SPECTRAL_TIMING_FLAG
#define SPECTRAL_TIMING(action) action
#else
#define SPECTRAL_TIMING(action)
#endif
#endif

/* Workspace for the eigendecomposition of one matrix cone at a time (cones.c):
 * one per thread of the PSD projection. */
typedef struct SCS_PSD_WORK ScsPsdWork;
/* LAPACK work-array sizes of one PSD block (cones.c) */
typedef struct SCS_PSD_LAPACK_SIZES ScsPsdLapackSizes;
#ifdef USE_LAPACK
struct SCS_PSD_WORK {
  /* Xs: the block, then its eigenvectors (syevd / heevd overwrite it with
   * them); e: the eigenvalues; Z: the result of the log-determinant and
   * sum-of-largest-eigenvalues projections (in psd_tw[0] when there are such
   * cones, NULL otherwise); the PSD projections put theirs in work / cwork */
  scs_float *Xs, *Z, *e, *work;
  scs_float *gamma; /* PSD diagonal-congruence factors */
  scs_complex_float *cXs, *cwork;
  blas_int *iwork;
  /* lengths of work, cwork and iwork: the largest sizes of the blocks that
   * use this workspace, and in psd_tw[0] at least lwork_syev for the
   * spectral cones (lcwork is 0 without a complex block) */
  blas_int lwork, lcwork, liwork;
};

struct SCS_PSD_LAPACK_SIZES {
  blas_int lwork;  /* work (real block) or rwork (complex block) */
  blas_int lcwork; /* complex work array, 0 for a real block */
  blas_int liwork; /* integer work array */
};
#endif

#ifdef SCS_FORCE_PSD_FALLBACK
/* test-only build flag (cones.c): number of PSD projections that fell back
 * to syev / heev */
extern scs_int SCS(psd_fallbacks);
#endif

/* private data to help cone projection step */
struct SCS_CONE_WORK {
  /*
   * cone_boundaries will contain array of indices of rows of A corresponding to
   * cone boundaries, boundaries[0] is starting index for cones of size larger
   * than 1
   */
  ScsCone *k; /* original cone information */
  scs_int *cone_boundaries;
  scs_int cone_boundaries_len;
  /* Per boundary entry: matrix dimension n if the block is a real PSD
   * cone (svec length n(n+1)/2), else 0. PSD blocks admit a richer
   * cone-invariant metric than a single scalar: any svec-diagonal
   * weight of the rank-one form w_ij = delta_i * delta_j (a diagonal
   * congruence X -> DXD). enforce_cone_boundaries fits this form for
   * PSD blocks instead of collapsing them to one value. */
  scs_int *cone_boundaries_psd_n;
  scs_int scaled_cones; /* boolean, whether the cones have been scaled */
  scs_float *s;         /* used for Moreau decomposition in projection */
  scs_int m;            /* total length of cone */
  /* box cone quantities */
  scs_float box_t_warm_start;
  scs_float *r_box_inv; /* precomputed 1/r_box[j] for the box projection */
  /* kinks of the box projection's derivative, 2 (bsize - 1) entries (see
   * proj_box_cone in cones.c); NULL without a box cone of size 2 or more */
  void *box_kinks;

  /* Eigendecomposition workspaces (with LAPACK), one per thread: psd_tw[t]
   * serves the PSD blocks of list t and is sized for them alone; psd_tw[0]
   * also serves the logdet and sum-of-largest-eigenvalues cones. */
  ScsPsdWork *psd_tw;
  /* PSD blocks: block b < ssize is real block b, else complex block b - ssize,
   * at offset psd_off[b] from the first PSD entry (psd_len entries in all).
   * List t of psd_nthreads is psd_blk[psd_tptr[t]] ..
   * psd_blk[psd_tptr[t + 1] - 1], the blocks of one thread (a smaller team
   * shares the lists out). psd_nthreads is 1 for a serial projection, one
   * list with all blocks in order, and 0 without PSD blocks. */
  scs_int psd_nthreads, psd_len;
  scs_int *psd_blk, *psd_tptr, *psd_off;
  /* Number of PSD blocks that a serial projection runs with the OpenMP
   * thread count and MKL's thread-local count at one (PSD_SERIAL_HOLD_N in
   * cones.c): with USE_PARALLEL_PSD, the blocks of 2 <= n <=
   * PSD_SERIAL_HOLD_N, also when psd_nthreads > 1 (the projection is serial
   * if the OpenMP thread count is lowered to one); 0 otherwise. */
  scs_int psd_nhold;
  /* With LAPACK, the work-array sizes of PSD block b are psd_lsz[b] (indexed
   * like psd_off), whichever workspace projects it (NULL without LAPACK). */
  ScsPsdLapackSizes *psd_lsz;

#ifdef USE_SPECTRAL_CONES
  /* if the projection onto the logarithmic cone should be warmstarted*/
  bool *log_cone_warmstarts;

  /* Needed for ell1 norm cone projection */
  Value_index *work_ell1;
  scs_float *work_ell1_proj;

  /* used for timing spectral vector cone and spectral matrix cone projections */
  SPECTRAL_TIMING(scs_float tot_time_mat_cone_proj;)
  SPECTRAL_TIMING(scs_float tot_time_vec_cone_proj;)

  /* workspace for singular value decompositions: */
  scs_float *s_nuc, *u_nuc, *vt_nuc, *work_nuc;
  blas_int lwork_nuc;

  /* the work-array size that the logdet and sum-of-largest projections pass
   * to syev (in psd_tw[0]->work): syev's workspace query for the largest of
   * these cones, which does not depend on the PSD blocks or their threads */
  blas_int lwork_syev;

  /* workspace that is used internally in the logdet projection (for example,
     the gradient and Hessian of the objective function in the projection
     problem are stored using this memory) */
  scs_float *work_logdet;

  /* workspace to store the projection onto the logarithm cone */
  scs_float *saved_log_projs;

  /* Stats for spectral projections, assuming there is only one spectral cone */
  Newton_stats newton_stats;

  /* workspace for projection onto sum-largest-evals cone */
  scs_float *work_sum_of_largest;
#endif
};

void SCS(free_cone)(ScsCone *k);
scs_int SCS(deep_copy_cone)(ScsCone *dest, const ScsCone *src);
ScsConeWork *SCS(init_cone)(ScsCone *k, scs_int m);
char *SCS(get_cone_header)(const ScsCone *k);
scs_int SCS(validate_cones)(const ScsData *d, const ScsCone *k);
scs_int SCS(proj_dual_cone)(scs_float *x, ScsConeWork *c,
                            const ScsScaling *scal, scs_float *r_y);
void SCS(finish_cone)(ScsConeWork *c);
void SCS(set_r_y)(const ScsConeWork *c, scs_float scale, scs_float *r_y);
void SCS(enforce_cone_boundaries)(const ScsConeWork *c, scs_float *vec,
                                  scs_float (*f)(const scs_float *, scs_int),
                                  scs_int psd_rank1);
void SCS(reset_cone_cache)(ScsConeWork *c);

#ifdef __cplusplus
}
#endif
#endif
