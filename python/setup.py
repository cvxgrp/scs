from distutils.core import setup, Extension
from glob import glob
from platform import system
from numpy import get_include
from numpy.distutils.system_info import get_info, BlasNotFoundError 

### if you're having errors linking blas/lapack, set this to false:
USE_LAPACK = True

lib = ['m']
if system() == 'Linux':
    lib += ['rt']

sources = ['scsmodule.c', ] + glob('../src/*.c')
define_macros = [('PYTHON', None), ('DLONG', None)] # add ('BLAS64', None) for 64 bit blas libs
include_dirs = ['../include', get_include()]
libraries = lib
extra_link_args = []
extra_compile_args = []

blas_info = get_info('blas_opt')
lapack_info = get_info('lapack_opt')
if blas_info and lapack_info and USE_LAPACK:
    include_dirs += blas_info.pop('include_dirs', []) + lapack_info.pop('include_dirs', [])
    define_macros += [('LAPACK_LIB_FOUND', None)] + blas_info.pop('define_macros', []) + lapack_info.pop('define_macros', [])
    libraries += blas_info.pop('libraries', []) + lapack_info.pop('libraries', [])
    extra_link_args += blas_info.pop('extra_link_args', []) + lapack_info.pop('extra_link_args', [])
    extra_compile_args += blas_info.pop('extra_compile_args', []) + lapack_info.pop('extra_compile_args', [])

_scs_direct = Extension(
                    name='_scs_direct',
                    sources=sources + glob('../linsys/direct/*.c') + glob('../linsys/direct/external/*.c'),
                    define_macros=define_macros,
                    include_dirs=include_dirs + ['../linsys/direct/', '../linsys/direct/external/'],
                    libraries=libraries,
                    extra_link_args=extra_link_args,
                    extra_compile_args=extra_compile_args
                    )

_scs_indirect = Extension(
                    name='_scs_indirect',
                    sources=sources + glob('../linsys/indirect/*.c'),
                    define_macros=define_macros + [('INDIRECT', None)],
                    include_dirs=include_dirs + ['../linsys/indirect/'],
                    libraries=libraries,
                    extra_link_args=extra_link_args,
                    extra_compile_args=extra_compile_args
                     )

# _scs_direct = Extension('_scs_direct', libraries=lib,
#                    # define LDL and AMD to use long ints
#                    # also define that we are building a python module
#                    define_macros=[
#                        ('PYTHON', None),
#                        ('DLONG', None)],
#                    include_dirs=['../include', get_include(),
#                        '../linsys/direct/',
#                        '../linsys/direct/external/'],
#                    sources=['scsmodule.c',
#                    ] + glob('../linsys/direct/*.c')
#                      + glob('../linsys/direct/external/*.c')
#                      + glob('../src/*.c'),
#                    extra_compile_args=[])

# _scs_indirect = Extension('_scs_indirect', libraries=lib,
#                    # define LDL and AMD to use long ints
#                    # also define that we are building a python module
#                    define_macros=[
#                        ('INDIRECT', None),
#                        ('PYTHON', None),
#                        ('DLONG', None)],
#                    include_dirs=['../include', get_include(),
#                        '../linsys/indirect/'],
#                    sources=['scsmodule.c',
#                    ] + glob('../linsys/indirect/*.c')
#                      + glob('../src/*.c'),
#                    extra_compile_args=[])

setup(name='scs',
        version='1.0',
        author = 'Brendan O\'Donoghue',
        author_email = 'bodonoghue85@gmail.com',
        url = 'http://github.com/cvxgrp/scs',
        description='This is the Python package for scs: splittling cone solver. See Github page for more information.',
        py_modules=['scs'],
        ext_modules=[_scs_direct, _scs_indirect],
        requires=["numpy (>= 1.7)"])
