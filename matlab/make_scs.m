% WARNING: OPENMP WITH MATLAB CAN CAUSE ERRORS AND CRASH, USE WITH CAUTION:
% openmp parallelizes the matrix multiply for the indirect solver (using CG):
flags.COMPILE_WITH_OPENMP = false;

% try to use blas libs; if not available then scs cannot solve SDPs:
try
    
    % EDIT THESE TO POINT TO YOUR BLAS + LAPACK LIBS:
    flags.BLASLIB = '-lopenblas -llapack -llapacke';
    flags.INCS = '-I/opt/local/include -I/usr/local/include';
    flags.LOCS = '-L/opt/local/lib -L/usr/local/lib -L/usr/lib';
    flags.LCFLAG = '-DLAPACK_LIB_FOUND';
    
    compile_direct(flags);
    compile_indirect(flags);
    
catch err
    
    flags.BLASLIB = '';
    flags.INCS = '';
    flags.LOCS = '';
    flags.LCFLAG = '';
    
    compile_direct(flags);
    compile_indirect(flags);    
    
    disp('Compiled without lapack support - unable to solve SDPs (can solve LPs, QPs, SOCPs, EXPs)')
    disp('To solve SDPs you must install cblas + lapacke and point the flags to the right locations')
end

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
