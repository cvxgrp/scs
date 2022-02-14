/* This file contains the outward facing SCS API. */
/* It includes all the input/output data structs and the API functions. */

#ifndef SCS_H_GUARD
#define SCS_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

/* Contains definitions of primitive types `scs_int` and `scs_float`. */
#include "scs_types.h"

#define SCS_NULL 0 /* NULL type */

/* The following abstract structs are implemented later. */

/** Struct containing acceleration workspace. Implemented by acceleration. */
typedef struct ACCEL_WORK AaWork;
/** Struct containing cone projection workspace. Implemented by cones. */
typedef struct SCS_CONE_WORK ScsConeWork;
/** Struct containing linear system workspace. Implemented by linear solver. */
typedef struct SCS_LIN_SYS_WORK ScsLinSysWork;
/** Struct containing SCS workspace. Implemented in `scs_work.h`. */
typedef struct SCS_WORK ScsWork;

/* SCS returns one of the following integer exit flags:            */
#define SCS_INFEASIBLE_INACCURATE (-7)
#define SCS_UNBOUNDED_INACCURATE (-6)
#define SCS_SIGINT (-5)
#define SCS_FAILED (-4)
#define SCS_INDETERMINATE (-3)
#define SCS_INFEASIBLE (-2) /* primal infeasible, dual unbounded   */
#define SCS_UNBOUNDED (-1)  /* primal unbounded, dual infeasible   */
#define SCS_UNFINISHED (0)  /* never returned, used as placeholder */
#define SCS_SOLVED (1)
#define SCS_SOLVED_INACCURATE (2)

/** This defines the data matrices which should be supplied in compressed
 *  sparse column format with zero based indexing.
 */
typedef struct {
  /** Matrix values, size: number of non-zeros. */
  scs_float *x;
  /** Matrix row indices, size: number of non-zeros. */
  scs_int *i;
  /** Matrix column pointers, size: `n+1`. */
  scs_int *p;
  /** Number of rows. */
  scs_int m;
  /** Number of columns. */
  scs_int n;
} ScsMatrix;

/** Struct containing all settings. */
typedef struct {
  /** Whether to heuristically rescale the data before solve. */
  scs_int normalize;
  /** Initial dual scaling factor (may be updated if adaptive_scale is on). */
  scs_float scale;
  /** Whether to adaptively update `scale`. */
  scs_int adaptive_scale;
  /** Primal constraint scaling factor. */
  scs_float rho_x;
  /** Maximum iterations to take. */
  scs_int max_iters;
  /** Absolute convergence tolerance. */
  scs_float eps_abs;
  /** Relative convergence tolerance. */
  scs_float eps_rel;
  /** Infeasible convergence tolerance. */
  scs_float eps_infeas;
  /** Douglas-Rachford relaxation parameter. */
  scs_float alpha;
  /** Time limit in secs (can be fractional). */
  scs_float time_limit_secs;
  /** Whether to log progress to stdout. */
  scs_int verbose;
  /** Whether to use warm start (put initial guess in ScsSolution struct). */
  scs_int warm_start;
  /** Memory for acceleration. */
  scs_int acceleration_lookback;
  /** Interval to apply acceleration. */
  scs_int acceleration_interval;
  /** String, if set will dump raw prob data to this file. */
  const char *write_data_filename;
  /** String, if set will log data to this csv file (makes SCS very slow). */
  const char *log_csv_filename;
} ScsSettings;

/** Struct containing problem data. */
typedef struct {
  /** A has `m` rows. */
  scs_int m;
  /** A has `n` cols, P has `n` cols and `n` rows. */
  scs_int n;
  /** A is supplied in CSC format (size `m` x `n`). */
  ScsMatrix *A;
  /** P is supplied in CSC format (size `n` x `n`), must be symmetric positive
   * semidefinite. Only pass in the upper triangular entries. If `P = 0` then
   * set `P = SCS_NULL`. */
  ScsMatrix *P;
  /** Dense array for b (size `m`). */
  scs_float *b;
  /** Dense array for c (size `n`). */
  scs_float *c;
} ScsData;

/** Cone data. Rows of data matrix `A` must be specified in this exact order. */
typedef struct {
  /** Number of linear equality constraints (primal zero, dual free). */
  scs_int z;
  /** Number of positive orthant cones. */
  scs_int l;
  /** Upper box values, `len(bu) = len(bl) = max(bsize-1, 0)`. */
  scs_float *bu;
  /** Lower box values, `len(bu) = len(bl) = max(bsize-1, 0)`. */
  scs_float *bl;
  /** Total length of box cone (includes scale `t`). */
  scs_int bsize;
  /** Array of second-order cone constraints, `len(q) = qsize`. */
  scs_int *q;
  /** Length of second-order cone array `q`. */
  scs_int qsize;
  /** Array of semidefinite cone constraints, `len(s) = ssize`. */
  scs_int *s;
  /** Length of semidefinite constraints array `s`. */
  scs_int ssize;
  /** Number of primal exponential cone triples. */
  scs_int ep;
  /** Number of dual exponential cone triples. */
  scs_int ed;
  /** Array of power cone params, must be in `[-1, 1]`, negative values are
   * interpreted as specifying the dual cone, `len(p) = psize ` */
  scs_float *p;
  /** Number of (primal and dual) power cone triples. */
  scs_int psize;
} ScsCone;

/** Contains primal-dual solution arrays or a certificate of infeasibility.
 *  Check the exit flag to determine whether this contains a solution or a
 *  certificate. If when passed into SCS the members `x`, `y`, `s` are
 *  NULL then SCS will allocate memory for them which should be managed
 *  by the user to prevent memory leaks.
 */
typedef struct {
  /** Primal variable. */
  scs_float *x;
  /** Dual variable. */
  scs_float *y;
  /** Slack variable. */
  scs_float *s;
} ScsSolution;

/** Contains information about the solve run at termination. */
typedef struct {
  /** Number of iterations taken. */
  scs_int iter;
  /** Status string, e.g. 'solved'. */
  char status[128];
  /** Linear system solver used. */
  char lin_sys_solver[128];
  /** Status as scs_int, defined in glbopts.h. */
  scs_int status_val;
  /** Number of updates to scale. */
  scs_int scale_updates;
  /** Primal objective. */
  scs_float pobj;
  /** Dual objective. */
  scs_float dobj;
  /** Primal equality residual. */
  scs_float res_pri;
  /** Dual equality residual. */
  scs_float res_dual;
  /** Duality gap. */
  scs_float gap;
  /** Infeasibility cert residual. */
  scs_float res_infeas;
  /** Unbounded cert residual. */
  scs_float res_unbdd_a;
  /** Unbounded cert residual. */
  scs_float res_unbdd_p;
  /** Time taken for setup phase (milliseconds). */
  scs_float setup_time;
  /** Time taken for solve phase (milliseconds). */
  scs_float solve_time;
  /** Final scale parameter. */
  scs_float scale;
  /** Complementary slackness. */
  scs_float comp_slack;
  /** Number of rejected AA steps. */
  scs_int rejected_accel_steps;
  /** Number of accepted AA steps. */
  scs_int accepted_accel_steps;
  /** Total time (milliseconds) spent in the linear system solver. */
  scs_float lin_sys_time;
  /** Total time (milliseconds) spent in the cone projection. */
  scs_float cone_time;
  /** Total time (milliseconds) spent in the acceleration routine. */
  scs_float accel_time;
} ScsInfo;

/*
 * Main library API.
 */

/**
 * Initialize SCS and allocate memory.
 *
 * All the inputs must be already allocated in memory before calling. After
 * this function returns then the memory associated with `d`, `k`, and `stgs`
 * can be freed as SCS maintains deep copies of these internally.
 *
 * It performs:
 * - data and settings validation
 * - problem data scaling
 * - automatic parameters tuning (if enabled)
 * - setup linear system solver:
 *      - direct solver: KKT matrix factorization is performed here
 *      - indirect solver: KKT matrix preconditioning is performed here.
 *
 *
 * @param  d      Problem data.
 * @param  k      Cone data.
 * @param  stgs   SCS solve settings.
 * @return        Solver workspace.
 */
ScsWork *scs_init(const ScsData *d, const ScsCone *k, const ScsSettings *stgs);

/**
 * Update the `b` vector, `c` vector, or both, before another solve call.
 *
 * After a solve we can reuse the SCS workspace in another solve if the only
 * problem data that has changed are the `b` and `c` vectors.
 *
 * @param  w            SCS workspace from scs_init (modified in-place).
 * @param  b            New `b` vector (can be `SCS_NULL` if unchanged).
 * @param  c            New `c` vector (can be `SCS_NULL` if unchanged).
 *
 * @return              0 if update successful.
 */
scs_int scs_update(ScsWork *w, scs_float *b, scs_float *c);

/**
 * Solve quadratic cone program initialized by scs_init.
 *
 * @param  w            Workspace allocated by scs_init.
 * @param  sol          Solution will be stored here. If members `x`, `y`, `s`
 *                      are NULL then SCS will allocate memory for them which
 *                      must be freed by the caller.
 * @param  info         Information about the solve will be stored here.
 * @param  warm_start   Whether to use the entries of `sol` as warm-start for
 *                      the solve.
 *
 * @return       Flag containing problem status (see \a glbopts.h).
 */
scs_int scs_solve(ScsWork *w, ScsSolution *sol, ScsInfo *info,
                  scs_int warm_start);

/**
 * Clean up allocated SCS workspace.
 *
 * @param  w  Workspace allocated by init, will be deallocated.
 */
void scs_finish(ScsWork *w);

/**
 * Solve quadratic cone program defined by data in d and cone k.
 *
 * All the inputs must already be allocated in memory before calling.
 *
 * @param  d      Problem data.
 * @param  k      Cone data.
 * @param  stgs   SCS solver settings.
 * @param  sol    Solution will be stored here. If members `x`, `y`, `s` are
 *                NULL then SCS will allocate memory for them.
 * @param  info   Information about the solve will be stored here.
 * @return        Flag containing problem status (see \a glbopts.h).
 */
scs_int scs(const ScsData *d, const ScsCone *k, const ScsSettings *stgs,
            ScsSolution *sol, ScsInfo *info);

/**
 * Helper function to set all settings to default values (see \a glbopts.h).
 *
 * @param  stgs   Settings struct that will be populated.
 */
void scs_set_default_settings(ScsSettings *stgs);

/**
 * Helper function simply returns the current version of SCS as a string.
 *
 * @return       SCS version as a string.
 */
const char *scs_version(void);

#ifdef __cplusplus
}
#endif
#endif
