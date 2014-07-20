% WARNING: OPENMP WITH MATLAB CAN CAUSE ERRORS AND CRASH, USE WITH CAUTION:
% openmp parallelizes the matrix multiply for the indirect solver (using CG):
flags.COMPILE_WITH_OPENMP = false;

flags.BLASLIB = '-lmwblas -lmwlapack';
flags.LCFLAG = '-DMATLAB_MEX_FILE -DLAPACK_LIB_FOUND -DDLONG';
flags.INCS = '';
flags.LOCS = '';

common_scs = '../src/linAlg.c ../src/cones.c ../src/cs.c ../src/util.c ../src/scs.c ../linsys/common.c scs_mex.c';
if (~isempty (strfind (computer, '64')))
    flags.arr = '-largeArrayDims';
else
    flags.arr = '';
end

if ( isunix && ~ismac )
    flags.link = '-lm -lrt';
elseif  ( ismac )
    flags.link = '-lm';
else
    flags.link = '';
    flags.LCFLAG = sprintf('-DNOBLASUNDERSCORE %s', flags.LCFLAG);
end

compile_direct(flags, common_scs);
compile_indirect(flags, common_scs);

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

