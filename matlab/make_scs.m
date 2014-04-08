% WARNING: OPENMP WITH MATLAB CAN CAUSE ERRORS AND CRASH, USE WITH CAUTION:
% openmp parallelizes the matrix multiply for the indirect solver (using CG):
flags.COMPILE_WITH_OPENMP = false;

flags.BLASLIB = '-lmwblas -lmwlapack';
flags.LCFLAG = '-DLAPACK_LIB_FOUND -DBLAS64';
flags.INCS = '';
flags.LOCS = '';

compile_direct(flags);
compile_indirect(flags);

%%
clear data cones
disp('Example run:');
m = 9;
n = 3;
data.A = sparse(randn(m,n));
data.b = randn(m,1);
data.c = randn(n,1);
cones.l = m;
[x,y,info] = scs_direct(data,cones,[]);
[x,y,info] = scs_indirect(data,cones,[]);
disp('SUCCESSFULLY INSTALLED SCS')

