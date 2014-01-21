function compileDirect(flags)

common_scs = '../src/linAlg.c ../src/cones.c ../src/cs.c ../src/util.c ../src/scs.c scs_mex.c';
if (~isempty (strfind (computer, '64')))
    d = '-fPIC';
    arr = '-largeArrayDims';
else
    d = '-fPIC -m32';
    arr = '';
end

cmd = sprintf ('mex -O %s CFLAGS="-std=c99 -O3 -DMATLAB_MEX_FILE %s %s" -I"../src/include" %s', arr, d, flags.LCFLAG, flags.INCS) ;
amd_files = {'amd_order', 'amd_dump', 'amd_postorder', 'amd_post_tree', ...
    'amd_aat', 'amd_2', 'amd_1', 'amd_defaults', 'amd_control', ...
    'amd_info', 'amd_valid', 'amd_global', 'amd_preprocess' } ;
for i = 1 : length (amd_files)
    cmd = sprintf ('%s ../direct/dirsrc/%s.c', cmd, amd_files {i}) ;
end
cmd = sprintf ('%s ../direct/dirsrc/ldl.c %s ../direct/private.c LINKLIBS="\\$LINKLIBS -lm %s %s" -output scs_direct', cmd, common_scs, flags.LOCS, flags.BLASLIB) ;
eval(cmd);
