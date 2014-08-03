function compile_direct(flags, common_scs)

cmd = sprintf ('mex -O %s %s %s -I.. -I../include %s', flags.arr, flags.LCFLAG, flags.INCS) ;
cmd = sprintf ('%s %s ../linsys/direct/private.c %s %s %s -output scs_direct', cmd, common_scs, flags.link, flags.LOCS, flags.BLASLIB) ;
eval(cmd);
