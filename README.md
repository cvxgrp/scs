SCS
====

[![Build Status](https://travis-ci.org/cvxgrp/scs.svg?branch=master)](https://travis-ci.org/cvxgrp/scs)
[![Build status](https://ci.appveyor.com/api/projects/status/4542u6kom5293qpm)](https://ci.appveyor.com/project/bodono/scs)

SCS = `splitting cone solver`

SCS is a C package for solving large-scale convex cone problems,
based on ["Operator Splitting for Conic Optimization via Homogeneous Self-Dual Embedding"](http://www.stanford.edu/~boyd/papers/scs.html) by 
Brendan O’Donoghue, Eric Chu, Neal Parikh, and Stephen Boyd

----
This code provides a solver for convex cone problems. It is an implementation
of the algorithm described in [this
paper](http://www.stanford.edu/~boyd/papers/scs.html). It provides both a
direct and an indirect solver in the form of a static library for inclusion in
other projects.

It simultaneously solves the primal cone program

```
minimize        c'*x 
subject to      A*x + s = b 
                s in K
```

and its dual

```
maximize        -b'*y 
subject to      -A'*y == c 
                y in K^*
```   

where `K` is a product cone of zero cones, linear cones `{ x | x >= 0 }`, 
second-order cones `{ (t,x) | ||x||_2 <= t }`, semi-definite cones `{ X | X psd }`,
and exponential cones `{(x,y,z) | y e^(x/y) <= z, y>0 }`.
`K^*` denotes the dual cone to `K`.

This package uses the LDL and AMD packages numerical linear
algebra packages by Tomothy Davis and others, the necessary files are included.
See [here](http://www.cise.ufl.edu/research/sparse/) for more information about
these packages.

----
This project was originally called `coneOS`; see previous version history
[here](https://github.com/cvxgrp/coneos).

This code is in alpha, if you're having trouble running this, or getting
erroneous results, please contact us.


Installing 
---------- 
Typing `make` at the command line
will produce two libaries, `libscsdir.a` and `libscsindir.a` found under the
`lib` folder. As a byproduct, it will also produce two demo binaries under the
`bin` folder called `demo_direct` and `demo_indirect`.

One caveat: if you have a 32-bit version of Matlab and use the build process
below (for Matlab), then if you try to make the libraries (on a 64-bit
machine), you must `make purge` before `make` again.

File an issue with us if the build process fails for you.

### Compiling a Matlab mex file
Running `make_scs` in Matlab under the
`matlab` folder will produce two usable mex files

If `make_scs` fails and complains about an incompatible architecture, edit the
`make_scs` file according to the comments.

Remember to include the `matlab` directory in your Matlab path if you wish to
use the mex file in your Matlab code. The calling sequence is (for the direct version):

	[x,y,s,info] = scs_direct(data,cones,params)

where data contains `A`, `b`, `c`  
params contains various options (see matlab file, can be empty)  
cones contains one or more of:  
+ `f` (num primal zero / dual free cones, i.e. primal equality constraints)
+ `l` (num linear cones)
+ `q` (array of SOCs sizes)
+ `s` (array of SDCs sizes)
+ `ep` (num primal exponential cones) 
+ `ed` (num dual exponential cones).

Type `help scs_direct` at the Matlab prompt to see its documentation.

### Installing a CVX solver
For users familiar with [CVX](http://cvxr.com), we supply a CVX shim which can be easily installed by invoking the following in the Matlab command line under the `matlab` directory.

    >> cvx_install_scs
    
You can select the scs solver with CVX as follows 
    
    >> cvx_solver 'scs'
    >> cvx_begin
    >> ... 
    >> cvx_end

### Using scs in Python

To create the Python interface, the following lines of code should work:
```
cd <scs-directory>/python
python setup.py install
```
You may need `sudo` privileges for a global installation. 

After installing the scs interface, you must import the module with
```
import scs
```
This module provides a single function `scs` with the following calling sequences:
```
sol = scs(data, cone, opts = None, USE_INDIRECT = False)
```
The argument `data` is a python dictionary with three elements `A`, `b`, and `c` where `b` and `c` are NUMPY arrays (i.e., matrices with a single column)  and `A` is a SCIPY *sparse* matrix in CSC format; if they are not of the proper format, scs will attempt to convert them. 

The argument `cone` is a dictionary with fields `f`, `l`, `q`, `s`, `ep` and `ed` (all of which are optional) corresponding to the supported cone types.

The argument `opts` is optional and is a dictionary with fields `MAX_ITERS`, `EPS`, `ALPHA`, `UNDET_TOL`, `VERBOSE`, and `NORMALIZE`. If `opts` is missing, then the algorithm uses default settings.

Finally set `USE_INDIRECT = True` to use the indirect linear equation solver.

The returned object is a dictionary containing the fields `sol['x']`, `sol['y']`, `sol['s']`, and `sol['info']`. 
The first four are NUMPY arrays containing the relevant solution. The last field contains a dictionary with the same fields as the `info` struct in the MATLAB interface.

Usage in C 
---------- 
If `make` completes successfully, it will produce two static library files,
`libscsdir.a` and `libscsindir.a` under the `lib` folder. To include the
libraries in your own source code, compile with the linker option with
`-L(PATH_TO_scs)\lib` and `-lscsdir` or `-lscsindir` (as needed).

These libraries (and `scs.h`) expose only four API functions:

* `Work * scs_init(Data * d, Cone * k, Info * info);`
    
    This initializes the Work struct containing the workspace that scs will use, and performs the necessary preprocessing (e.g. matrix factorization).

* `idxint scs_solve(Work * w, Data * d, Cone * k, Sol * sol, Info * info);`
    
    This solves the problem as defined by Data and Cone using workspace in w. The solution is returned in sol and information about the solve is retu.rned in info. None of the inputs can be NULL. You can call scs_solve many times for one call to scs_init, so long as the matrix A does not change (b and c can change).

* `void scs_finish(Data * d, Work * w);`
    
    Called after all solves completed, to free data and cleanup.
* `idxint scs(Data * d, Cone * k, Sol * sol, Info * info);`

    Simply calls all the above routines, so can't use reuse the workspace (for, e.g., factorization caching).

The five relevant data structures are:
```
    typedef struct PROBLEM_DATA Data;
    typedef struct SOL_VARS Sol;
    typedef struct INFO Info;
    typedef struct CONE Cone;

    /* defined in linSys.h, can be overriden by user */
    typedef struct A_DATA_MATRIX AMatrix;


    /* this struct defines the data matrix A */
    struct A_DATA_MATRIX {
        /* A is supplied in column compressed format */
        pfloat * x; /* A values, size: NNZ A */
        idxint * i; /* A row index, size: NNZ A */
        idxint * p; /* A column pointer, size: n+1 */
    };

    /* struct that containing standard problem data */
    struct PROBLEM_DATA {
    	/* problem dimensions */
    	idxint m, n;        /* A has m rows, n cols*/
   
        AMatrix * A;        /* A is supplied in data format specified by linsys solver */
    	pfloat * b, *c;     /* dense arrays for b (size m), c (size n) */
    
    	/* other input parameters: default suggested input */
    	idxint MAX_ITERS;   /* maximum iterations to take: 2500 */
    	pfloat EPS;         /* convergence tolerance: 1e-3 */
    	pfloat ALPHA;       /* relaxation parameter: 1.8 */
    	pfloat RHO_X;       /* x equality constraint scaling: 1e-3 */
    	pfloat SCALE;       /* if normalized, rescales by this factor: 1 */
        pfloat CG_RATE;     /* for indirect, tolerance goes down like (1/iter)^CG_RATE: 2 */
        idxint VERBOSE;     /* boolean, write out progress: 1 */
    	idxint NORMALIZE;   /* boolean, heuristic data rescaling: 1 */
    	idxint WARM_START;  /* boolean, warm start with guess in Sol struct: 0 */
    };
    
    /* contains primal-dual solution arrays */
    struct SOL_VARS {
    	pfloat * x, *y, *s;
    };
    
    /* contains terminating information */
    struct INFO {
    	idxint iter;        /* number of iterations taken */
    	char status[32];    /* status string, e.g. Solved */
    	idxint statusVal;   /* status as idxint, defined below */
    	pfloat pobj;        /* primal objective */
    	pfloat dobj;        /* dual objective */
    	pfloat resPri;      /* primal equality residual */
    	pfloat resDual;     /* dual equality residual */
    	pfloat relGap;      /* relative duality gap */
    	pfloat setupTime;   /* time taken for setup phase */
    	pfloat solveTime;   /* time taken for solve phase */
    };
   
    struct CONE {
        idxint f;           /* number of primal linear equality constraints */
        idxint l;           /* length of LP cone */
        idxint *q;   	    /* array of second-order cone constraints */
        idxint qsize;       /* length of SOC array */
    	idxint *s;			/* array of SD constraints */
    	idxint ssize;		/* length of SD array */
        idxint ep;          /* number of primal exponential cone triples */
        idxint ed;          /* number of dual exponential cone triples */
    };
```        
The data matrix `A` is specified in column-compressed format, and the vectors
`b` and `c` are specified as dense arrays. The solutions `x` (primal), `s`
(slack), and `y` (dual) are returned as dense arrays. Cones are specified as
the struct above, the rows of `A` must correspond to the cones in the
exact order as specified by the cone struct (i.e. put linear cones before
soc cones etc.).

### Warm Start
You can warm-start (supply a guess of the solution) by setting WARM_START in Data to 1 and supplying the warm-starts in the Sol struct (x,y and s). These are used to initialize the iterates in scs_solve.
 
### Re-using matrix factorization
To factorize the matrix once and solve many times, simply call scs_init once, and use scs_solve many times with the same workspace, changing the input data (and optionally warm-starts) for each iteration. See run_scs.c for an example.

### Using your own linear system solver
Simply implement all the methods and the two structs in `include/linSys.h` and plug it in.

Solving SDPs
---------- 
In order to solve SDPs you must have BLAS and LAPACK installed.
Point `scs.mk` to the location of these libraries. Without
these you can still solve SOCPs, LPs, and ECPs.

Scalability
----------- 
Note that this code is merely meant as an
implementation of the ideas in our paper. The actual code does not use more
than a single CPU. Nevertheless, for problems that fit in memory on a single
computer, this code will (attempt to) solve them.

To scale this solver, one must either provide a distributed solver for linear
systems or a distributed matrix-vector multiplication.
