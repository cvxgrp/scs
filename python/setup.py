"""Setup SCS."""
from __future__ import print_function
from distutils import msvccompiler
import setuptools
import glob
import platform
import shutil
import tempfile
import argparse
import os
import subprocess
import sys
from setuptools.command.build_ext import build_ext
import numpy as np
from numpy.distutils import system_info

env_lib_dirs = os.environ.get('BLAS_LAPACK_LIB_PATHS', [])
env_libs = os.environ.get('BLAS_LAPACK_LIBS', [])

if env_lib_dirs or env_libs:
  print('using environment variables for blas/lapack libraries')
  env_vars = {}
  if env_lib_dirs:
    env_vars['library_dirs'] = [env_lib_dirs]
  if env_libs:
    env_vars['libraries'] = env_libs.split(':')
  install_scs(blas_info=env_vars, lapack_info={})
  return

# environment variables not set, using defaults instead
blas_info = system_info.get_info('blas_opt')
# ugly hack due to scipy bug
if 'libraries' in blas_info:
  if 'mkl_intel_lp64' in blas_info['libraries']:
    blas_info = system_info.get_info('blas_mkl')
if not blas_info:
  blas_info = system_info.get_info('blas')
print(blas_info)

lapack_info = system_info.get_info('lapack_opt')
# ugly hack due to scipy bug
if 'libraries' in lapack_info:
  if 'mkl_intel_lp64' in lapack_info['libraries']:
    lapack_info = system_info.get_info('lapack')
if not lapack_info:
  lapack_info = system_info.get_info('lapack')
print(lapack_info)

blas_info = kwargs['blas_info']
lapack_info = kwargs['lapack_info']

extra_compile_args = ['-O3']
library_dirs = []
extra_link_args = []
libraries = []
extra_define_macros = []
sources = [
    'scsmodule.c',
] + glob.glob('../src/*.c') + ('../src/linsys/*.c')
include_dirs = [
    '../include',
    '../linsys',
    np.get_include(),
]
define_macros = [('PYTHON', None), ('CTRLC', 1), ('COPYAMATRIX', None),
                 ('DLONG', 1)]

if platform.system() == 'Linux':
  libraries += ['rt']
define_macros += [('LAPACK_LIB_FOUND', None)] + blas_info.pop(
    'define_macros', []) + lapack_info.pop('define_macros', [])
include_dirs += blas_info.pop('include_dirs', []) + lapack_info.pop(
    'include_dirs', [])
library_dirs += blas_info.pop('library_dirs', []) + lapack_info.pop(
    'library_dirs', [])
libraries += blas_info.pop('libraries', []) + lapack_info.pop('libraries', [])
extra_link_args += blas_info.pop('extra_link_args', []) + lapack_info.pop(
    'extra_link_args', [])
extra_compile_args += blas_info.pop('extra_compile_args', []) + lapack_info.pop(
    'extra_compile_args', [])

_scs_direct = setuptools.Extension(
    name='_scs_direct',
    sources=sources + glob.glob(os.path.join(root_dir, 'linsys/direct/*.c')) +
    glob.glob(os.path.join(root_dir, 'linsys/direct/external/*.c')),
    define_macros=define_macros + extra_define_macros,
    include_dirs=include_dirs + [
        os.path.join(root_dir, 'linsys/direct/'),
        os.path.join(root_dir, 'linsys/direct/external/')
    ],
    library_dirs=library_dirs,
    libraries=libraries,
    extra_link_args=extra_link_args,
    extra_compile_args=extra_compile_args)

_scs_indirect = setuptools.Extension(
    name='_scs_indirect',
    sources=sources + glob.glob(os.path.join(root_dir, 'linsys/indirect/*.c')),
    define_macros=define_macros + extra_define_macros + [('INDIRECT', None)],
    include_dirs=include_dirs + [os.path.join(root_dir, 'linsys/indirect/')],
    library_dirs=library_dirs,
    libraries=libraries,
    extra_link_args=extra_link_args,
    extra_compile_args=extra_compile_args)

ext_modules = [_scs_direct, _scs_indirect]

setuptools.setup(
    name='scs',
    version='1.2.7',
    author='Brendan O\'Donoghue',
    author_email='bodonoghue85@gmail.com',
    url='http://github.com/cvxgrp/scs',
    description='scs: splitting conic solver',
    py_modules=['scs'],
    ext_modules=ext_modules,
    cmdclass={'build_ext': build_ext_scs},
    install_requires=['numpy >= 1.7', 'scipy >= 0.13.2'],
    license='MIT',
    zip_safe=False,
    long_description=
    """Solves convex cone programs via operator splitting. Can solve:
     linear programs (LPs), second-order cone programs (SOCPs),
     semidefinite programs (SDPs), exponential cone programs (ECPs),
     and power cone programs (PCPs), or problems with any combination of
     those cones. See http://github.com/cvxgrp/scs for more details.""")
