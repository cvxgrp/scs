SCS
====

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

The argument `opts` is optional and is a dictionary with fields `MAX_ITERS`, `EPS`, `ALPHA`, `VERBOSE`, and `NORMALIZE`. If `opts` is missing, then the algorithm uses default settings.

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

* `scs_int scs_solve(Work * w, Data * d, Cone * k, Sol * sol, Info * info);`
    
    This solves the problem as defined by Data and Cone using workspace in w. The solution is returned in sol and information about the solve is retu.rned in info. None of the inputs can be NULL. You can call scs_solve many times for one call to scs_init, so long as the matrix A does not change (b and c can change).

* `void scs_finish(Data * d, Work * w);`
    
    Called after all solves completed, to free data and cleanup.
* `scs_int scs(Data * d, Cone * k, Sol * sol, Info * info);`

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
        scs_float * x; /* A values, size: NNZ A */
        scs_int * i; /* A row index, size: NNZ A */
        scs_int * p; /* A column pointer, size: n+1 */
    };

    /* struct that containing standard problem data */
    struct PROBLEM_DATA {
    	/* problem dimensions */
    	scs_int m, n;        /* A has m rows, n cols*/
   
        AMatrix * A;        /* A is supplied in data format specified by linsys solver */
    	scs_float * b, *c;     /* dense arrays for b (size m), c (size n) */
    
    	/* other input parameters: default suggested input */
    	scs_int MAX_ITERS;   /* maximum iterations to take: 2500 */
    	scs_float EPS;         /* convergence tolerance: 1e-3 */
    	scs_float ALPHA;       /* relaxation parameter: 1.8 */
    	scs_float RHO_X;       /* x equality constraint scaling: 1e-3 */
    	scs_float SCALE;       /* if normalized, rescales by this factor: 1 */
        scs_float CG_RATE;     /* for indirect, tolerance goes down like (1/iter)^CG_RATE: 2 */
        scs_int VERBOSE;     /* boolean, write out progress: 1 */
    	scs_int NORMALIZE;   /* boolean, heuristic data rescaling: 1 */
    	scs_int WARM_START;  /* boolean, warm start with guess in Sol struct: 0 */
    };
    
    /* contains primal-dual solution arrays */
    struct SOL_VARS {
    	scs_float * x, *y, *s;
    };
    
    /* contains terminating information */
    struct INFO {
    	scs_int iter;        /* number of iterations taken */
    	char status[32];    /* status string, e.g. Solved */
    	scs_int statusVal;   /* status as scs_int, defined below */
    	scs_float pobj;        /* primal objective */
    	scs_float dobj;        /* dual objective */
    	scs_float resPri;      /* primal equality residual */
    	scs_float resDual;     /* dual equality residual */
    	scs_float relGap;      /* relative duality gap */
    	scs_float setupTime;   /* time taken for setup phase */
    	scs_float solveTime;   /* time taken for solve phase */
    };
   
    struct CONE {
        scs_int f;           /* number of primal linear equality constraints */
        scs_int l;           /* length of LP cone */
        scs_int *q;   	    /* array of second-order cone constraints */
        scs_int qsize;       /* length of SOC array */
    	scs_int *s;			/* array of SD constraints */
    	scs_int ssize;		/* length of SD array */
        scs_int ep;          /* number of primal exponential cone triples */
        scs_int ed;          /* number of dual exponential cone triples */
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
