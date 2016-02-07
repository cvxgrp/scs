SCS
====

[![Build Status](https://travis-ci.org/cvxgrp/scs.svg?branch=master)](https://travis-ci.org/cvxgrp/scs)
[![Build status](https://ci.appveyor.com/api/projects/status/4542u6kom5293qpm)](https://ci.appveyor.com/project/bodono/scs)

SCS (`splitting conic solver`) is a numerical optimization package for solving
large-scale convex cone problems, based on our paper [Conic Optimization via Operator Splitting
and Homogeneous Self-Dual Embedding](http://www.stanford.edu/~boyd/papers/scs.html). It is written in C
and can be used in other C, C++, Python, Matlab, R, Julia, Java, and Scala
programs via included interfaces (Julia interface available
[here](https://github.com/JuliaOpt/SCS.jl)). It can also be called as a solver from
convex optimization toolboxes [CVX](http://cvxr.com/cvx/) (3.0 or later),
[CVXPY](https://github.com/cvxgrp/cvxpy),
[Convex.jl](https://github.com/JuliaOpt/Convex.jl), and
[Yalmip](https://github.com/johanlofberg/YALMIP).

----
SCS numerically solves convex cone programs using the alternating direction
method of multipliers ([ADMM](http://web.stanford.edu/~boyd/papers/admm_distr_stats.html)).
It returns solutions to both the primal and dual problems if the problem
is feasible, or a certificate of infeasibility otherwise. It solves
the following primal cone problem:

```
minimize        c'x
subject to      Ax + s = b
                s in K
```
over variables `x` and `s`, where `A`, `b` and `c` are user-supplied data and `K` is a user-defined convex cone.
The dual problem is given by
```
maximize        -b'y
subject to      -A'y == c
                y in K^*
```
over variable `y`, where `K^*` denotes the dual cone to `K`.

The cone `K` can be any Cartesian product of the following primitive cones:
+ zero cone `{x | x = 0 }` (dual to the free cone `{x | x in R}`)
+ positive orthant `{x | x >= 0}`
+ second-order cone `{(t,x) | ||x||_2 <= t}`
+ positive semidefinite cone `{ X | min(eig(X)) >= 0, X = X^T }`
+ exponential cone `{(x,y,z) | y e^(x/y) <= z, y>0 }`
+ dual exponential cone `{(u,v,w) | −u e^(v/u) <= e w, u<0}`
+ power cone `{(x,y,z) | x^a * y^(1-a) >= |z|, x>=0, y>=0}`
+ dual power cone `{(u,v,w) | (u/a)^a * (v/(1-a))^(1-a) >= |w|, u>=0, v>=0}`

The rows of the data matrix `A` correspond to the cones in `K`.
**The rows of `A` must be in the order of the cones given above, i.e., first come the
rows that correspond to the zero/free cones, then those that correspond to the
positive orthants, then SOCs, etc.** For a `k` dimensional semidefinite cone
when interpreting the rows of the data matrix `A`
SCS assumes that the `k x k` matrix variable has been vectorized by scaling the
off-diagonal entries by `sqrt(2)` and stacking the **lower triangular elements column-wise**
to create a vector of length `k(k+1)/2`. See the section on semidefinite programming
below.

At termination SCS returns solution `(x*, s*, y*)` if the problem is feasible,
or a certificate of infeasibility otherwise. See [here](http://web.stanford.edu/~boyd/cvxbook/)
for (much) more details about cone programming and certificates of infeasibility.

**Semidefinite Programming**

SCS assumes that the matrix variables and the input data corresponding to
semidefinite cones have been vectorized by **scaling the off-diagonal entries by
`sqrt(2)`** and stacking the lower triangular elements **column-wise**. For a `k x k`
matrix variable (or data matrix) this operation would create a vector of length
`k(k+1)/2`. Scaling by `sqrt(2)` is required to preserve the inner-product.

**To recover the matrix solution this operation must be inverted on the
components of the vector returned by SCS corresponding to semidefinite cones**.
That is, the off-diagonal entries must be scaled by `1/sqrt(2)` and the upper
triangular entries are filled in by copying the values of lower triangular
entries.

More explcitly, we want to express
`Tr(C X)` as `vec(C)'*vec(X)`, where the `vec` operation takes the `k x k` matrix
```
X = [ X11 X12 ... X1k
      X21 X22 ... X2k
      ...
      Xk1 Xk2 ... Xkk ]
```
and produces a vector consisting of
```
vec(X) = (X11, sqrt(2)*X21, ..., sqrt(2)*Xk1, X22, sqrt(2)*X32, ..., Xkk).
```


**Linear equation solvers**

Each iteration of SCS requires the solution of a set of linear equations.  This
package includes two implementations for solving linear equations: a direct
solver which uses a cached LDL factorization and an indirect solver based on
conjugate gradients. The indirect solver can be run on either the cpu or
gpu.

The direct solver uses the LDL and AMD packages numerical linear
algebra packages by Timothy Davis and others, the necessary files are included.
See [here](http://www.cise.ufl.edu/research/sparse/) for more information about
these packages.

### Using SCS in C
Typing `make` at the command line will compile the code and create
SCS libraries in the `out` folder. It will also produce six demo binaries in the
`out` folder named `demo_direct`, `demo_indirect`, `demo_gpu`, `demo_SOCP_direct`,
`demo_SOCP_indirect`, and `demo_SOCP_gpu`.

If `make` completes successfully, it will produce three static library files,
`libscsdir.a`, `libscsindir.a`, and `libscsgpu.a` under the `out` folder and
three dynamic library files `libscsdir.ext`, `libscsindir.ext`, and
`libscsgpu.ext` (where `.ext` extension is platform dependent) in the same
folder. To include the libraries in your own source code, compile with the
linker option with `-L(PATH_TO_SCS)\lib` and `-lscsdir` or `-lscsindir` (as
needed).

These libraries (and `scs.h`) expose only four API functions:

* `Work * scs_init(const Data * d, const Cone * k, Info * info);`
    
    This initializes the Work struct containing the workspace that scs will
    use, and performs the necessary preprocessing (e.g. matrix factorization).
    All inputs `d`, `k`, and `info` must be memory allocated before calling.

* `scs_int scs_solve(Work * w, const Data * d, const Cone * k, Sol * sol, Info * info);`
    
    This solves the problem as defined by Data `d` and Cone `k` using the workspace
    in `w`. The solution is returned in `sol` and information about the solve
    is returned in `info` (outputs must have memory allocated before calling).
    None of the inputs can be NULL. You can call `scs_solve` many times for one
    call to `scs_init`, so long as the matrix A
    does not change (b and c can change).

* `void scs_finish(Work * w);`
    
    Called after all solves completed to free allocated memory and other cleanup.

* `scs_int scs(const Data * d, const Cone * k, Sol * sol, Info * info);`

    Convenience method that simply calls all the above routines in order, for cases where
    the workspace does not need to be reused. All inputs must have memory allocated
    before this call.

The relevant data structures are:
```C

    typedef struct SCS_PROBLEM_DATA Data;
    typedef struct SCS_SETTINGS Settings;
    typedef struct SCS_SOL_VARS Sol;
    typedef struct SCS_INFO Info;
    typedef struct SCS_SCALING Scaling;
    typedef struct SCS_WORK Work;
    typedef struct SCS_CONE Cone;

    /* defined in linSys.h, can be overriden by user */
    typedef struct A_DATA_MATRIX AMatrix;

    /* this struct defines the data matrix A */
    struct A_DATA_MATRIX {
        /* A is supplied in column compressed format */
        scs_float * x;      /* A values, size: NNZ A */
        scs_int * i;        /* A row index, size: NNZ A */
        scs_int * p;        /* A column pointer, size: n+1 */
        scs_int m, n;       /* m rows, n cols */
    };
    
    /* struct containing problem data */
    struct SCS_PROBLEM_DATA {
        /* these cannot change for multiple runs for the same call to scs_init */
        scs_int m, n;       /* A has m rows, n cols */
        AMatrix * A;        /* A is supplied in data format specified by linsys solver */

        /* these can change for multiple runs for the same call to scs_init */
        scs_float * b, *c;  /* dense arrays for b (size m), c (size n) */

        Settings * stgs;    /* contains solver settings specified by user */
    };

    /* struct containing solver settings */
    struct SCS_SETTINGS {
        /* settings parameters: default suggested input */

        /* these *cannot* change for multiple runs with the same call to scs_init */
        scs_int normalize;  /* boolean, heuristic data rescaling: 1 */
        scs_float scale;    /* if normalized, rescales by this factor: 5 */
        scs_float rho_x;    /* x equality constraint scaling: 1e-3 */

        /* these can change for multiple runs with the same call to scs_init */
        scs_int max_iters;  /* maximum iterations to take: 2500 */
        scs_float eps;      /* convergence tolerance: 1e-3 */
        scs_float alpha;    /* relaxation parameter: 1.8 */
        scs_float cg_rate;  /* for indirect, tolerance goes down like (1/iter)^cg_rate: 2 */
        scs_int verbose;    /* boolean, write out progress: 1 */
        scs_int warm_start; /* boolean, warm start (put initial guess in Sol struct): 0 */
    };   

    /* contains primal-dual solution arrays */
    struct SCS_SOL_VARS {
        scs_float * x, *y, *s;
    };

    /* contains terminating information */
    struct SCS_INFO {
        scs_int iter;       /* number of iterations taken */
        char status[32];    /* status string, e.g. 'Solved' */
        scs_int statusVal;  /* status as scs_int, defined below */
        scs_float pobj;     /* primal objective */
        scs_float dobj;     /* dual objective */
        scs_float resPri;   /* primal equality residual */
        scs_float resDual;  /* dual equality residual */
        scs_float resInfeas;/* infeasibility cert residual */
        scs_float resUnbdd; /* unbounded cert residual */
        scs_float relGap;   /* relative duality gap */
        scs_float setupTime;/* time taken for setup phase (milliseconds) */
        scs_float solveTime;/* time taken for solve phase (milliseconds) */
    };


    /* contains normalization variables */
    struct SCS_SCALING {
        scs_float *D, *E; /* for normalization */
        scs_float meanNormRowA, meanNormColA;
    };

    /* NB: rows of data matrix A must be specified in this exact order */
    struct SCS_CONE {
        scs_int f;          /* number of linear equality constraints */
        scs_int l;          /* length of LP cone */
        scs_int *q;         /* array of second-order cone constraints */
        scs_int qsize;      /* length of SOC array */
        scs_int *s;         /* array of SD constraints */
        scs_int ssize;      /* length of SD array */
        scs_int ep;         /* number of primal exponential cone triples */
        scs_int ed;         /* number of dual exponential cone triples */
        scs_int psize;      /* number of (primal and dual) power cone triples */
        scs_float * p;      /* array of power cone params, must be \in [-1, 1],
                               negative values are interpreted as specifying the dual cone */
    };

    /* SCS returns one of the following integers: (zero never returned)     */
    #define SCS_INFEASIBLE_INACCURATE   (-7)
    #define SCS_UNBOUNDED_INACCURATE    (-6)
    #define SCS_SIGINT                  (-5)
    #define SCS_FAILED                  (-4)
    #define SCS_INDETERMINATE           (-3)
    #define SCS_INFEASIBLE              (-2) /* primal infeasible, dual unbounded   */
    #define SCS_UNBOUNDED               (-1) /* primal unbounded, dual infeasible   */
    #define SCS_UNFINISHED              (0)  /* never returned, used as placeholder */
    #define SCS_SOLVED                  (1)
    #define SCS_SOLVED_INACCURATE       (2)
```

The types `scs_float` and `scs_int` can be specified by the user, they default
to `double` and `int` respectively.

The data matrix `A` is specified in column-compressed format and the vectors
`b` and `c` are specified as dense arrays. The solutions `x` (primal), `s`
(slack), and `y` (dual) are returned as dense arrays. Cones are specified as
the struct above, the rows of `A` must correspond to the cones in the
exact order as specified by the cone struct (i.e. put linear cones before
second-order cones etc.).

**Warm-start**

You can warm-start SCS (supply a guess of the solution) by setting warm_start in
Data to `1` and supplying the warm-starts in the Sol struct (`x`,`y`, and `s`). All
inputs must be warm-started if any one is. These
are used to initialize the iterates in `scs_solve`.
 
**Re-using matrix factorization**

If using the direct version you can factorize the matrix once and solve many
times. Simply call scs_init once, and use `scs_solve` many times with the same
workspace, changing the input data `b` and `c` (and optionally warm-starts) for
each iteration. See run_scs.c for an example.

**Using your own linear system solver**

To use your own linear system solver simply implement all the methods and the
two structs in `include/linSys.h` and plug it in.

**Solving SDPs**

In order to solve SDPs you must have BLAS and LAPACK installed.
Edit `scs.mk` to set `USE_LAPACK = 1` and point to the location of these
libraries. Without these you can still solve problems using the other cones.

### Using SCS in Matlab
Running `make_scs` in Matlab under the `matlab` folder will produce two mex
files, one for each of the direct and indirect solvers.

Remember to include the `matlab` directory in your Matlab path if you wish to
use the mex file in your Matlab code. The calling sequence is (for the direct version):

	[x,y,s,info] = scs_direct(data,cones,settings)

where data is a struct containing `A`, `b`, and `c`, settings is a struct containing 
solver options (see matlab file, can be empty),
and cones is a struct that contains one or more of:
+ `f`  (num primal zero / dual free cones, i.e. primal equality constraints)
+ `l`  (num linear cones)
+ `q`  (array of SOCs sizes)
+ `s`  (array of SDCs sizes)
+ `ep` (num primal exponential cones)
+ `ed` (num dual exponential cones)
+ `p`  (array of primal/dual power params).

Type `help scs_direct` at the Matlab prompt to see its documentation.

### Using SCS in Python

To install SCS as a python package type:
```shell
cd <scs-directory>/python
python setup.py install
```
You may need `sudo` privileges for a global installation. Running SCS requires numpy and
scipy to be installed. You can install the gpu interface using
```shell
python setup.py install --scs --gpu
```

After installing the SCS interface, you import SCS using
```python
import scs
```
This module provides a single function `scs` with the following call signature:
```python
sol = scs(data, cone, [use_indirect=false, gpu=false, verbose=true, normalize=true, max_iters=2500, scale=5, eps=1e-3, cg_rate=2, alpha=1.8, rho_x=1e-3])
```
Arguments in the square brackets are optional, and default to the values on the right of their respective equals signs.
The argument `data` is a python dictionary with three elements `A`, `b`, and
`c` where `b` and `c` are NUMPY arrays (i.e., matrices with a single column)
and `A` is a SCIPY **sparse matrix in CSC format**; if they are not of the proper
format, SCS will attempt to convert them.

The argument `cone` is a dictionary with fields `f`, `l`, `q`, `s`, `ep`,
`ed`, and `p` (all of which are optional) corresponding to the supported cone types.

The returned object is a dictionary containing the fields `sol['x']`, `sol['y']`, `sol['s']`, and `sol['info']`.
The first three are NUMPY arrays containing the relevant solution. The last field contains a dictionary with solver information.

### Using SCS in Java / Scala

SCS can be called from Java and Scala via the Java Native Interface (JNI).
Typing `make` in the `java` directory will compile the native libraries into
`java/bin` and build `scs.jar`. You will need to add these files to your java
projects to use SCS. To test a random cone program type `make testproblem`.

To use with eclipse follow the above steps, then import `scs.jar` into eclipse and
link the dynamic libraries (`libscsjdir.*` `libscsjindir.*`):

```
Project > Properties > Build Path > Open scs.jar > Native Library > Edit > Select the folder with libraries
```

To solve a problem create a new instance of `ConeProgram` with constructor

```Java
public ConeProgram(Data d, Cone k, Settings p, IConeSolver solver);
```

where the `IConeSolver` interface can be one of `DirectSolver` or `IndirectSolver`.
Then call

```Java
ConeProgram.solve();
```

on your instance of `ConeProgram`, which will return an instance of
`Solution`, containing the solution and information about the run.

### Using SCS in R

To install SCS as an R packge cd into the `r` subdirectory and type
```shell
./install.sh
```
and to test type `R --no-save < demo/randomLp.R` (this is also called by the install script).

To call SCS use the function `scs`, defined as
```R
scs <- function(A, b, c, cone, params)
```
where `A` is a `Matrix` object, `b` and `c` are vectors, and `cone` and `params`
are lists specifying the cone dimensions and input params respectively.

### Using SCS in Julia
See usage instructions [here](https://github.com/JuliaOpt/SCS.jl).

