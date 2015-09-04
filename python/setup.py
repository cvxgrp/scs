from __future__ import print_function
import os
from setuptools import setup, Extension
from glob import glob
from platform import system
from numpy import get_include
from numpy.distutils.system_info import get_info, BlasNotFoundError 

def install_scs(USE_64_BIT_BLAS, blas_info, lapack_info, USE_OPENMP, rootDir): 
    libraries = []
    if system() == 'Linux':
        libraries += ['rt']
   
    sources = ['scsmodule.c', ] + glob(rootDir + 'src/*.c') + glob(rootDir + 'linsys/*.c')
    include_dirs = [rootDir, rootDir + 'include', get_include(), rootDir + 'linsys']
    
    define_macros = [('PYTHON', None), ('DLONG', None), ('CTRLC', 1), ('COPYAMATRIX', None)]
    # define_macros += [('EXTRAVERBOSE', 999)] # for debugging
    extra_compile_args = ["-O3"]
    library_dirs = []
    extra_link_args = []
    
    if USE_OPENMP:
        define_macros += [('OPENMP', None)]
        extra_compile_args += ['-fopenmp']
        extra_link_args += ['-lgomp']
    
    if USE_64_BIT_BLAS:
        define_macros += [('BLAS64', None)]
    
    if blas_info or lapack_info:
        define_macros += [('LAPACK_LIB_FOUND', None)] + blas_info.pop('define_macros', []) + lapack_info.pop('define_macros', [])
        include_dirs += blas_info.pop('include_dirs', []) + lapack_info.pop('include_dirs', [])
        library_dirs += blas_info.pop('library_dirs', []) + lapack_info.pop('library_dirs', [])
        libraries += blas_info.pop('libraries', []) + lapack_info.pop('libraries', [])
        extra_link_args += blas_info.pop('extra_link_args', []) + lapack_info.pop('extra_link_args', [])
        extra_compile_args += blas_info.pop('extra_compile_args', []) + lapack_info.pop('extra_compile_args', [])
    
    _scs_direct = Extension(
                        name='_scs_direct',
                        sources=sources + glob(rootDir + 'linsys/direct/*.c') + glob(rootDir + 'linsys/direct/external/*.c'),
                        define_macros=define_macros,
                        include_dirs=include_dirs + [rootDir + 'linsys/direct/', rootDir + 'linsys/direct/external/'],
                        library_dirs=library_dirs,
                        libraries=libraries,
                        extra_link_args=extra_link_args,
                        extra_compile_args=extra_compile_args
                        )
    
    _scs_indirect = Extension(
                        name='_scs_indirect',
                        sources=sources + glob(rootDir + 'linsys/indirect/*.c'),
                        define_macros=define_macros + [('INDIRECT', None)],
                        include_dirs=include_dirs + [rootDir + 'linsys/indirect/'],
                        library_dirs=library_dirs,
                        libraries=libraries,
                        extra_link_args=extra_link_args,
                        extra_compile_args=extra_compile_args
                        )
    setup(name='scs',
            version='1.1.6',
            author = 'Brendan O\'Donoghue',
            author_email = 'bodonoghue85@gmail.com',
            url = 'http://github.com/cvxgrp/scs',
            description='scs: splitting conic solver',
            py_modules=['scs'],
            ext_modules=[_scs_direct, _scs_indirect],
            install_requires=["numpy >= 1.7","scipy >= 0.13.2"],
            license = "MIT",
            long_description="Solves convex cone programs via operator splitting. Can solve: linear programs (LPs), second-order cone programs (SOCPs), semidefinite programs (SDPs), exponential cone programs (ECPs), and power cone programs (PCPs), or problems with any combination of those cones. See http://github.com/cvxgrp/scs for more details."
            )

#####################################

# location of SCS root directory, containing 'src/' etc.
rootDir = '../'

# use 'export OMP_NUM_THREADS=16' to control num of threads (in that case, 16)
USE_OPENMP = False

# set to true if linking against blas/lapack libraries that use longs instead of ints for indices:
USE_64_BIT_BLAS = False

BLAS_LAPACK_LIB_PATHS='BLAS_LAPACK_LIB_PATHS'
BLAS_LAPACK_LIBS='BLAS_LAPACK_LIBS'

env_lib_dirs = os.environ.get(BLAS_LAPACK_LIB_PATHS, [])
env_libs = os.environ.get(BLAS_LAPACK_LIBS, [])

if env_lib_dirs or env_libs:
    print("using environment variables for blas/lapack libraries")
    env_vars = {}
    if env_lib_dirs:
        env_vars['library_dirs'] = env_lib_dirs.split(':')
    if env_libs:
        env_vars['libraries'] = env_libs.split(':') 
    install_scs(USE_64_BIT_BLAS=USE_64_BIT_BLAS, blas_info=env_vars, lapack_info={},  USE_OPENMP=USE_OPENMP, rootDir=rootDir)
else:
    # environment variables not set, using defaults instead
    try:
        print("using blas_opt / lapack_opt")
        install_scs(USE_64_BIT_BLAS=USE_64_BIT_BLAS, blas_info=get_info('blas_opt'), lapack_info=get_info('lapack_opt'), USE_OPENMP=USE_OPENMP, rootDir=rootDir)
    except SystemExit as e: # catch permission denied error
        print("SystemExit")
        print(e)
    except:
        # print("error:", sys.exc_info()[0])
        print("blas_opt / lapack_opt install failed, trying blas / lapack")
        try:
            install_scs(USE_64_BIT_BLAS=USE_64_BIT_BLAS, blas_info=get_info('blas'), lapack_info=get_info('lapack'), USE_OPENMP=USE_OPENMP, rootDir=rootDir)
        except:
            install_scs(USE_64_BIT_BLAS=USE_64_BIT_BLAS, blas_info={}, lapack_info={}, USE_OPENMP=USE_OPENMP, rootDir=rootDir)
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

