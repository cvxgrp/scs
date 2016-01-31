function compile_gpu(flags, common_scs)

flags.link = sprintf('-lcudart -lcublas -lcusparse %s', flags.link);
flags.INCS = sprintf('-I/usr/local/cuda/include %s', flags.INCS);
if ( ismac )
flags.link = sprintf('-L/usr/local/cuda/lib %s', flags.link);
else
% TODO probably not right for windows
flags.link = sprintf('-L/usr/local/cuda/lib64 %s', flags.link);
end

% compile gpu
if (flags.COMPILE_WITH_OPENMP)
    cmd = sprintf('mex -O %s %s %s -DOPENMP COMPFLAGS="/openmp \\$COMPFLAGS" CFLAGS="\\$CFLAGS -fopenmp" ../linsys/gpu/private.c %s -I.. -I../include %s %s %s -output scs_gpu',  flags.arr, flags.LCFLAG, common_scs, flags.INCS, flags.link, flags.LOCS, flags.BLASLIB);
else
    cmd = sprintf('mex -O %s %s %s ../linsys/gpu/private.c %s -I.. -I../include %s %s %s -output scs_gpu',  flags.arr, flags.LCFLAG, common_scs, flags.INCS, flags.link, flags.LOCS, flags.BLASLIB);
end
eval(cmd);
