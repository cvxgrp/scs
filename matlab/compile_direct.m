function compile_direct(flags, common_scs)

if (flags.COMPILE_WITH_OPENMP)
    cmd = sprintf('mex -O %s %s %s %s -DOPENMP COMPFLAGS="/openmp \\$COMPFLAGS" CFLAGS="\\$CFLAGS -fopenmp" -I.. -I../include %s', flags.arr, flags.LCFLAG, flags.INCS, flags.INT);
else
    cmd = sprintf ('mex -O %s %s %s %s -I.. -I../include %s', flags.arr, flags.LCFLAG, flags.INCS, flags.INT);
end
cmd = sprintf ('%s %s ../linsys/direct/private.c %s %s %s -output scs_direct', cmd, common_scs, flags.link, flags.LOCS, flags.BLASLIB) ;
eval(cmd);
