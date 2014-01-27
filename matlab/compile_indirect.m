function compile_indirect(flags)

common_scs = '../src/linAlg.c ../src/cones.c ../src/cs.c ../src/util.c ../src/scs.c scs_mex.c';
if (~isempty (strfind (computer, '64')))
    d = '-fPIC -DDLONG';
    arr = '-largeArrayDims';
else
    d = '-fPIC -m32';
    arr = '';
end

% compile indirect
if (flags.COMPILE_WITH_OPENMP)
    cmd = sprintf('mex -O %s COMPFLAGS="/openmp \\$COMPFLAGS" CFLAGS="\\$CFLAGS -std=c99 -O3 -fopenmp -pthread -DMATLAB_MEX_FILE %s %s" ../linsys/indirect/private.c %s -I../include %s -output scs_indirect "LINKLIBS="\\$LINKLIBS -lm  %s %s"',  arr, d, flags.LCFLAG, common_scs, flags.INCS, flags.LOCS, flags.BLASLIB);
else
    cmd = sprintf('mex -O %s CFLAGS="-std=c99 -O3 -pthread -DMATLAB_MEX_FILE %s %s" ../linsys/indirect/private.c %s -I../include %s -output scs_indirect LINKLIBS="\\$LINKLIBS -lm  %s %s"',  arr, d, flags.LCFLAG, common_scs, flags.INCS, flags.LOCS, flags.BLASLIB);
end
eval(cmd);
