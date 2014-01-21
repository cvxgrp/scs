clear all
make_scs

[ ver, isoctave, fs, ps ] = cvx_version;
scs_install_path = strcat(fs,'/scs');
mkdir(scs_install_path);

%scs_bin_loc = sprintf('%s/',scs_install_path);
copyfile('scs_direct.m*',scs_install_path);
copyfile('scs_indirect.m*',scs_install_path);

copyfile('cvx_scs.m', strcat(fs,'/shims/cvx_scs.m'));
copyfile('scs.m', scs_install_path);
cvx_setup
