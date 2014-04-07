function compile_indirect(flags)

common_scs = '../src/linAlg.c ../src/cones.c ../src/cs.c ../src/util.c ../src/scs.c scs_mex.c';
d = '-fPIC -DDLONG';
if (~isempty (strfind (computer, '64')))
    arr = '-largeArrayDims';
else
    arr = '';
end
if ( isunix && ~ismac ) 
    link = '-lm -lrt';
else
    link = '-lm';
end


% compile indirect
if (flags.COMPILE_WITH_OPENMP)
    cmd = sprintf('mex -O %s COMPFLAGS="/openmp \\$COMPFLAGS" CFLAGS="\\$CFLAGS -O3 -fopenmp -pthread -DMATLAB_MEX_FILE %s %s" ../linsys/indirect/private.c %s -I../include %s %s %s %s -output scs_indirect',  arr, d, flags.LCFLAG, common_scs, flags.INCS, link, flags.LOCS, flags.BLASLIB);
else
    cmd = sprintf('mex -O %s CFLAGS="-O3 -DMATLAB_MEX_FILE %s %s" ../linsys/indirect/private.c %s -I../include %s %s %s %s -output scs_indirect',  arr, d, flags.LCFLAG, common_scs, flags.INCS, link, flags.LOCS, flags.BLASLIB);
end
eval(cmd);
