from distutils.core import setup, Extension
from glob import glob
from platform import system
from numpy import get_include

lib = []
if system() == 'Linux':
    lib += ['rt']

_ecos = Extension('_ecos', libraries = lib,
                    # define LDL and AMD to use long ints
                    # also define that we are building a python module
                    define_macros = [
                        ('PYTHON',None),
                        ('DLONG', None),
                        ('LDL_LONG', None)],
                    include_dirs = ['../include', get_include(),
                        '../external/amd/include',
                        '../external/ldl/include',
                        '../external/SuiteSparse_config'],
                    sources = ['ecosmodule.c',
                        '../external/ldl/src/ldl.c'
                    ] + glob('../external/amd/src/*.c')
                      + glob('../src/*.c'))


setup(  name = 'ecos',
        version = '1.0.1',
        description = 'This is the Python package for ECOS: Embedded Cone Solver. See Github page for more information.',
        py_modules = ['ecos'],
        ext_modules = [_ecos],
        requires = ["numpy (>= 1.7)"])
