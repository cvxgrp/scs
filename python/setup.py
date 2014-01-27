from distutils.core import setup, Extension
from glob import glob
from platform import system
from numpy import get_include

lib = ['m']
if system() == 'Linux':
    lib += ['rt']

_scs_direct = Extension('_scs_direct', libraries = lib,
                    # define LDL and AMD to use long ints
                    # also define that we are building a python module
                    define_macros = [
                        ('PYTHON',None),
                        ('DLONG', None)],
                    include_dirs = ['../include', get_include(),
                        '../linsys/direct/',
                        '../linsys/direct/external/'],
                    sources = ['scsmodule.c',
                    ] + glob('../linsys/direct/*.c')
                      + glob('../linsys/direct/external/*.c')
                      + glob('../src/*.c'),
                    extra_compile_args=['-std=c99'])

_scs_indirect = Extension('_scs_indirect', libraries = lib,
                    # define LDL and AMD to use long ints
                    # also define that we are building a python module
                    define_macros = [
                        ('INDIRECT',None),
                        ('PYTHON',None),
                        ('DLONG', None)],
                    include_dirs = ['../include', get_include(),
                        '../linsys/indirect/'],
                    sources = ['scsmodule.c',
                    ] + glob('../linsys/indirect/*.c')
                      + glob('../src/*.c'),
                    extra_compile_args=['-std=c99'])

setup(  name = 'scs',
        version = '1.0',
        description = 'This is the Python package for scs: splittling cone solver. See Github page for more information.',
        py_modules = ['scs'],
        ext_modules = [_scs_direct, _scs_indirect],
        requires = ["numpy (>= 1.7)"])
