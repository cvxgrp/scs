#ifndef CONSTANTS_H_GUARD
#define CONSTANTS_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

/* SCS VERSION NUMBER ----------------------------------------------    */
#define SCS_VERSION                                                            \
    ("1.2.6") /* string literals automatically null-terminated */

/* SCS returns one of the following integers:                           */
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

/* DEFAULT SOLVER PARAMETERS AND SETTINGS --------------------------    */
#define MAX_ITERS (2500)
#define EPS (1E-3)
#define ALPHA (1.5)
#define RHO_X (1E-3)
#define SCALE (1.0)
#define CG_RATE (2.0)
#define VERBOSE (1)
#define NORMALIZE (1)
#define WARM_START (0)

#ifdef __cplusplus
}
#endif
#endif
