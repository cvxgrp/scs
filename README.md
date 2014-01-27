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
The argument `data` is a python dictionary with three elements `A`, `b`, and `c` where `c` and `b` are NUMPY arrays (i.e., matrices with a single column).  and `A` is a SCIPY *sparse* matrix in CSC format; if they are not of the proper format, scs will attempt to convert them. 

The argument `cone` is a dictionary with fields `f`, `l`, `q`, `s`, `ep` and `ed` (any of which are optional) corresponding to the supported cone types.

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

These libraries (and `scs.h`) expose only three API functions:

* idxint scs(Data * d, Cone * k, Sol * sol, Info * info); 
    
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
      idxint n, m; /* problem dimensions */
      /* problem data, A, b, c: */
      pfloat * Ax;
      idxint * Ai, * Ap;
      idxint Anz;
      pfloat * b, * c;
      idxint MAX_ITERS;
      pfloat EPS_ABS, ALPH, UNDET_TOL, RHO_X;
      idxint VERBOSE, NORMALIZE;  // boolean
    } Data;
    
    /* contains primal-dual solution vectors */
    typedef struct SOL_VARS {
      pfloat * x, * y, *s;
    } Sol;
    
    /* contains terminating information */
    typedef struct INFO {
    	idxint iter;
    	char status[16];
    	idxint stint; // status as int
        pfloat pobj;
    	pfloat dobj;
    	pfloat resPri;
    	pfloat resDual;
    	pfloat relGap;
    	pfloat time;
    } Info;
    
    typedef struct Cone_t {
        idxint f;          /* number of linear equality constraints */
        idxint l;          /* length of LP cone */
        idxint *q;         /* array of second-order cone constraints */
        idxint qsize;      /* length of SOC array */
        idxint *s;         /* array of SD constraints */
        idxint ssize;      /* length of SD array */
        idxint ep;         /* number of triples in exponential cone */
        idxint ed;         /* number of triples in dual exponential cone */
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
Poidxint `scs.mk` to the location of these libraries. Without
these you can still solve SOCPs, LPs, and EXPs.

Scalability
----------- 
Note that this code is merely meant as an
implementation of the ideas in our paper. The actual code does not use more
than a single CPU. Nevertheless, for problems that fit in memory on a single
computer, this code will (attempt to) solve them.

To scale this solver, one must either provide a distributed solver for linear
systems or a distributed matrix-vector multiplication.
