scs 
====

scs = `splitting cone solver`

scs is a C package for solving large-scale convex cone problems,
based on ["Operator Splitting for Conic Optimization via Homogeneous Self-Dual Embedding"](http://www.stanford.edu/~boyd/papers/scs.html) by 
Brendan O’Donoghue, Eric Chu, Neal Parikh, and Stephen Boyd

----
This code provides a solver for convex cone problems. It is an implementation
of the algorithm described in [this
paper](http://www.stanford.edu/~boyd/papers/scs.html). It provides both a
direct and an indirect solver in the form of a static library for inclusion in
other projects.

It simultaneously solves the primal cone program

	minimize     c'*x subject to   A*x + s == b s in K 
                 
and its dual

	maximize     -b'*y subject to   -A'*y == c y in K^* 

where `K` is a product cone of free cones, linear cones `{ x | x >= 0 }`, 
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
use the mex file in your Matlab code. The calling sequence is

	[x,y,info] = scs_direct(data,cones,params)

where data contains `A`, `b`, `c`
params contains various options (see matlab file, can be empty)
cones contains one or more of:
+ `f` (num free/zero cones)
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


Usage in C 
---------- 
If `make` completes successfully, it will produce two static library files,
`libscsdir.a` and `libscsindir.a` under the `lib` folder. To include the
libraries in your own source code, compile with the linker option with
`-L(PATH_TO_scs)\lib` and `-lscsdir` or `-lscsindir` (as needed).

These libraries (and `scs.h`) expose only three API functions:

* int scs(Data * d, Cone * k, Sol * sol, Info * info); 
    
	This solves the problem specified in the `Data` and `Cone` structures,
    and returns the solution in the Sol struct and various information about the run in
    the Info struct.

* void freeData(Data \* d, Cone \* k)
    
	This frees the `Data` and `Cone` structures.
    
* void freeSol(Sol \* sol)

	This frees the `Sol` structure.
    
The four relevant data structures are:

    
    /* struct that containing standard problem data */
    typedef struct PROBLEM_DATA {
      int n, m; /* problem dimensions */
      /* problem data, A, b, c: */
      double * Ax;
      int * Ai, * Ap;
      int Anz;
      double * b, * c;
      int MAX_ITERS;
      double EPS_ABS, ALPH, UNDET_TOL, RHO_X;
      int VERBOSE, NORMALIZE;  // boolean
    } Data;
    
    /* contains primal-dual solution vectors */
    typedef struct SOL_VARS {
      double * x, * y, *s;
    } Sol;
    
    /* contains terminating information */
    typedef struct INFO {
    	int iter;
    	char status[16];
    	int stint; // status as int
        double pobj;
    	double dobj;
    	double resPri;
    	double resDual;
    	double relGap;
    	double time;
    } Info;
    
    typedef struct Cone_t {
        int f;          /* number of linear equality constraints */
        int l;          /* length of LP cone */
        int *q;         /* array of second-order cone constraints */
        int qsize;      /* length of SOC array */
        int *s;         /* array of SD constraints */
        int ssize;      /* length of SD array */
        int ep;         /* number of triples in exponential cone */
        int ed;         /* number of triples in dual exponential cone */
    } Cone;
        
    
The data matrix `A` is specified in column-compressed format, and the vectors
`b` and `c` are specified as dense arrays. The solutions `x` (primal), `s`
(slack), and `y` (dual) are returned as dense arrays. Cones are specified as
the struct above, the rows of `A` must correspond to the cones in the
exact order as specified by the cone struct (i.e. put linear cones before
soc cones etc.).

Solving SDPs
---------- 
In order to solve SDPs you must have BLAS and LAPACK installed.
Point scs.mk to the location of these libraries. Without
these you can still solve SOCPs, LPs, and EXPs.

Scalability
----------- 
Note that this code is merely meant as an
implementation of the ideas in our paper. The actual code does not use more
than a single CPU. Nevertheless, for problems that fit in memory on a single
computer, this code will (attempt to) solve them.

To scale this solver, one must either provide a distributed solver for linear
systems or a distributed matrix-vector multiplication.
