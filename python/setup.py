from __future__ import print_function
from glob import glob
from numpy import get_include
from numpy.distutils.system_info import get_info, BlasNotFoundError
from platform import system
from setuptools import setup, Extension
import argparse
import os
import sys

SCS_ARG_MARK = '--scs'

parser = argparse.ArgumentParser(description='Compilation args for SCS.')
parser.add_argument(SCS_ARG_MARK, dest='scs', action='store_true',
                    default=False, help='Put this first to ensure following arguments are parsed correctly')
parser.add_argument('--gpu', dest='gpu', action='store_true',
                    default=False, help='Also compile the GPU CUDA version of SCS')
parser.add_argument('--openmp', dest='openmp', action='store_true',
                    default=False, help='Compile with OpenMP multi-threading')
parser.add_argument('--float', dest='float32', action='store_true',
                    default=False, help='Use 32 bit (single precision) floats, default is 64 bit')
parser.add_argument('--extraverbose', dest='extraverbose', action='store_true',
                    default=False, help='Extra verbose SCS (for debugging)')
parser.add_argument('--int', dest='int32', action='store_true',
                    default=False, help='Use 32 bit ints, default is 64 bit (GPU code always uses 32 bit ints)')
parser.add_argument('--blas64', dest='blas64', action='store_true',
                    default=False, help='Use 64 bit ints for the blas/lapack libs')
parser.add_argument('--root_dir', dest='root_dir', action='store',
                    default='../', help='Location of SCS src directory')
args, unknown = parser.parse_known_args()

env_lib_dirs = os.environ.get('BLAS_LAPACK_LIB_PATHS', [])
env_libs = os.environ.get('BLAS_LAPACK_LIBS', [])

root_dir = args.root_dir
print(args)

# necessary to remove SCS args before passing to setup:
if SCS_ARG_MARK in sys.argv:
    sys.argv = sys.argv[0:sys.argv.index(SCS_ARG_MARK)]

def run_install():
    if env_lib_dirs or env_libs:
        print("using environment variables for blas/lapack libraries")
        env_vars = {}
        if env_lib_dirs:
            env_vars['library_dirs'] = env_lib_dirs.split(':')
        if env_libs:
            env_vars['libraries'] = env_libs.split(':')
        install_scs(blas_info=env_vars, lapack_info={})
        return

    # environment variables not set, using defaults instead
    try:
        print("using blas_opt / lapack_opt")
        install_scs(blas_info=get_info('blas_opt'), lapack_info=get_info('lapack_opt'))
        return
    except:
        pass # fall back to blas / lapack (not opt)

    print("blas_opt / lapack_opt install failed, trying blas / lapack")
    try:
        install_scs(blas_info=get_info('blas'), lapack_info=get_info('lapack'))
    except:
        install_scs(blas_info={}, lapack_info={})
        print("###############################################################################################")
        print("# failed to find blas/lapack libs, SCS cannot solve SDPs but can solve LPs, SOCPs, ECPs, PCPs #")
        print("# install blas/lapack and run this install script again to allow SCS to solve SDPs            #")
        print("#                                                                                             #")
        print("# scs will use environment variables BLAS_LAPACK_LIB_PATHS and BLAS_LAPACK_LIBS if set        #")
        print("# use this to link against blas/lapack libs that scs can't find on it's own, usage ex:        #")
        print("#        >> export BLAS_LAPACK_LIB_PATHS=/usr/lib/:/other/dir                                 #")
        print("#        >> export BLAS_LAPACK_LIBS=blas:lapack                                               #")
        print("#        >> python setup.py install                                                           #")
        print("###############################################################################################")


def install_scs(**kwargs):
    blas_info = kwargs['blas_info']
    lapack_info = kwargs['lapack_info']

    extra_compile_args = ["-O3"]
    library_dirs = []
    extra_link_args = []
    libraries = []
    extra_define_macros = []
    sources = ['scsmodule.c', ] + glob(root_dir + 'src/*.c') + glob(root_dir + 'linsys/*.c')
    include_dirs = [root_dir, root_dir + 'include', get_include(), root_dir + 'linsys']
    define_macros = [('PYTHON', None), ('CTRLC', 1), ('COPYAMATRIX', None)]

    if system() == 'Linux':
        libraries += ['rt']
    if args.float32:
        define_macros += [('FLOAT', 1)] # single precision floating point
    if args.extraverbose:
        define_macros += [('EXTRAVERBOSE', 999)] # for debugging
    if args.openmp:
        define_macros += [('OPENMP', 1)] # openMP multi-threading
        library_dirs += ['/usr/local/gfortran/lib'] # TODO not for all systems
        extra_compile_args += ['-fopenmp']
        extra_link_args += ['-lgomp']
    if args.blas64:
        define_macros += [('BLAS64', 1)] # 64 bit blas
    if blas_info or lapack_info:
        define_macros += [('LAPACK_LIB_FOUND', None)] + blas_info.pop('define_macros', []) + lapack_info.pop('define_macros', [])
        include_dirs += blas_info.pop('include_dirs', []) + lapack_info.pop('include_dirs', [])
        library_dirs += blas_info.pop('library_dirs', []) + lapack_info.pop('library_dirs', [])
        libraries += blas_info.pop('libraries', []) + lapack_info.pop('libraries', [])
        extra_link_args += blas_info.pop('extra_link_args', []) + lapack_info.pop('extra_link_args', [])
        extra_compile_args += blas_info.pop('extra_compile_args', []) + lapack_info.pop('extra_compile_args', [])
    if not args.int32:
        extra_define_macros += [('DLONG', 1)] # longs for integer type

    _scs_direct = Extension(
                        name='_scs_direct',
                        sources=sources + glob(root_dir + 'linsys/direct/*.c') + glob(root_dir + 'linsys/direct/external/*.c'),
                        define_macros=define_macros + extra_define_macros,
                        include_dirs=include_dirs + [root_dir + 'linsys/direct/', root_dir + 'linsys/direct/external/'],
                        library_dirs=library_dirs,
                        libraries=libraries,
                        extra_link_args=extra_link_args,
                        extra_compile_args=extra_compile_args
                        )

    _scs_indirect = Extension(
                        name='_scs_indirect',
                        sources=sources + glob(root_dir + 'linsys/indirect/*.c'),
                        define_macros=define_macros + extra_define_macros + [('INDIRECT', None)],
                        include_dirs=include_dirs + [root_dir + 'linsys/indirect/'],
                        library_dirs=library_dirs,
                        libraries=libraries,
                        extra_link_args=extra_link_args,
                        extra_compile_args=extra_compile_args
                        )

    ext_modules = [_scs_direct, _scs_indirect]

    if args.gpu:
        _scs_gpu = Extension(
                        name='_scs_gpu',
                        sources=sources + glob(root_dir + 'linsys/gpu/*.c'),
                        define_macros=define_macros + [('GPU', None)],
                        include_dirs=include_dirs + [root_dir + 'linsys/gpu/', '/usr/local/cuda/include'],
                        library_dirs=library_dirs + ['/usr/local/cuda/lib', '/usr/local/cuda/lib64'], # TODO probably not right for windows
                        libraries=libraries + ['cudart', 'cublas', 'cusparse'],
                        extra_link_args=extra_link_args,
                        extra_compile_args=extra_compile_args
                        )
        ext_modules += [_scs_gpu]

    setup(name='scs',
            version='1.2.2',
            author = 'Brendan O\'Donoghue',
            author_email = 'bodonoghue85@gmail.com',
            url = 'http://github.com/cvxgrp/scs',
            description='scs: splitting conic solver',
            py_modules=['scs'],
            ext_modules=ext_modules,
            install_requires=["numpy >= 1.7","scipy >= 0.13.2"],
            license = "MIT",
            zip_safe=False,
            long_description="Solves convex cone programs via operator splitting. Can solve: linear programs (LPs), second-order cone programs (SOCPs), semidefinite programs (SDPs), exponential cone programs (ECPs), and power cone programs (PCPs), or problems with any combination of those cones. See http://github.com/cvxgrp/scs for more details."
            )

run_install()
