from __future__ import print_function
from distutils.core import setup, Extension
from distutils.errors import *
#import os
from glob import glob
from platform import system
from numpy import get_include
from numpy.distutils.system_info import get_info, BlasNotFoundError 

def install_scs(USE_LAPACK, USE_64_BIT_BLAS, BLAS_STR, LAPACK_STR, USE_OPENMP, rootDir): 
    libraries = []
    if system() == 'Linux':
        libraries += ['rt']
    
    sources = ['scsmodule.c', ] + glob(rootDir + 'src/*.c') + glob(rootDir + 'linsys/*.c')
    include_dirs = [rootDir, rootDir + 'include', get_include(), rootDir + 'linsys']
    
    define_macros = [('PYTHON', None), ('DLONG', None)]
    extra_compile_args = ["-O3"]
    library_dirs = []
    extra_link_args = []
    
    if USE_OPENMP:
        define_macros += [('OPENMP', None)]
        extra_compile_args += ['-fopenmp']
        extra_link_args += ['-lgomp']
    
    if USE_64_BIT_BLAS:
        define_macros += [('BLAS64', None)]
    
    blas_info = get_info(BLAS_STR)
    lapack_info = get_info(LAPACK_STR)
    if blas_info and lapack_info and USE_LAPACK:
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
            version='1.0.5',
            author = 'Brendan O\'Donoghue',
            author_email = 'bodonoghue85@gmail.com',
            url = 'http://github.com/cvxgrp/scs',
            description='scs: splittling cone solver',
            py_modules=['scs'],
            ext_modules=[_scs_direct, _scs_indirect],
            requires=["numpy (>= 1.7)","scipy (>= 0.13.2)"],
            license = "GPLv3",
            long_description="Solves convex cone programs via operator splitting. Can solve: linear programs (LPs) second-order cone programs (SOCPs), semidefinite programs (SDPs), and exponential cone programs (ECPs). See http://github.com/cvxgrp/scs for more details."
            )

#####################################

# location of SCS root directory, containing 'src/' etc.
rootDir = '../'

# use 'export OMP_NUM_THREADS=16' to control num of threads (in that case use 16)
USE_OPENMP = False

# set to true if linking against blas/lapack libraries that use longs instead of ints for indices:
USE_64_BIT_BLAS = False

try:
    install_scs(USE_LAPACK=True, USE_64_BIT_BLAS=USE_64_BIT_BLAS, BLAS_STR='blas_opt', LAPACK_STR='lapack_opt', USE_OPENMP=USE_OPENMP, rootDir=rootDir)
except:
    try:
        install_scs(USE_LAPACK=True, USE_64_BIT_BLAS=USE_64_BIT_BLAS, BLAS_STR='blas', LAPACK_STR='lapack', USE_OPENMP=USE_OPENMP, rootDir=rootDir)
    except:
        install_scs(USE_LAPACK=False, USE_64_BIT_BLAS=USE_64_BIT_BLAS, BLAS_STR='', LAPACK_STR='', USE_OPENMP=USE_OPENMP, rootDir=rootDir)
        print("#############################################################################################")
        print("# failed to find blas/lapack libs, SCS cannot solve SDPs but can solve LPs, SOCPs, and ECPs #")
        print("# install blas/lapack and run this install script again to allow SCS to solve SDPs          #")
        print("#############################################################################################")

