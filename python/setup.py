from __future__ import print_function
from distutils.msvccompiler import MSVCCompiler
from glob import glob
from numpy import get_include
from numpy.distutils.system_info import get_info, BlasNotFoundError
from platform import system
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import shutil
import tempfile
import argparse
import os
import subprocess
import sys

SCS_ARG_MARK = '--scs'

parser = argparse.ArgumentParser(description='Compilation args for SCS.')
parser.add_argument(SCS_ARG_MARK, dest='scs', action='store_true',
                    default=False, help='Put this first to ensure following arguments are parsed correctly')
parser.add_argument('--gpu', dest='gpu', action='store_true',
                    default=False, help='Also compile the GPU CUDA version of SCS')
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
            env_vars['library_dirs'] = [env_lib_dirs]
        if env_libs:
            env_vars['libraries'] = env_libs.split(':')
        install_scs(blas_info=env_vars, lapack_info={})
        return

    # environment variables not set, using defaults instead
    blas_info = get_info('blas_opt')
    # ugly hack due to scipy bug
    if 'libraries' in blas_info:
        if 'mkl_intel_lp64' in blas_info['libraries']:
            blas_info = get_info('blas_mkl')
    if not blas_info:
        blas_info = get_info('blas')
    print(blas_info)

    lapack_info = get_info('lapack_opt')
    # ugly hack due to scipy bug
    if 'libraries' in lapack_info:
        if 'mkl_intel_lp64' in lapack_info['libraries']:
            lapack_info = get_info('lapack')
    if not lapack_info:
        lapack_info = get_info('lapack')
    print(lapack_info)

    try:
        install_scs(blas_info=blas_info, lapack_info=lapack_info)
    except:
        pass #TODO fix

def can_compile_with_openmp(cc, flags, gomp_paths):
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    # attempt to compile a test program with openmp
    test_filename_base = 'test_openmp'
    test_filename = test_filename_base + '.c'
    with open(test_filename, 'w') as f:
        f.write(
        '#include <omp.h>\n'
        '#include <stdio.h>\n'
        'int main(void) {\n'
        '  #pragma omp parallel\n'
        '  printf("thread: %d\\n", omp_get_thread_num());\n'
        '  return 0;\n'
        '}')

    compilation_cmd = [cc] + flags
    for path in gomp_paths:
        compilation_cmd  += ['-L' + path]
    compilation_cmd += [test_filename] + ['-o'] + [test_filename_base]
    run_cmd = ['./' + test_filename_base]

    try:
        with open(os.devnull, 'w') as fnull:
            subprocess.call(compilation_cmd,
                stdout=fnull, stderr=fnull)
            exit_code = subprocess.call(run_cmd,
                stdout=fnull, stderr=fnull)
    except OSError:
        exit_code = 1

    # clean up
    os.chdir(curdir)
    shutil.rmtree(tmpdir)
    return exit_code == 0


class build_ext_scs(build_ext):
    def build_extensions(self):
        # TODO: include paths for other systems
        gomp_paths = ['/usr/local/gfortran/lib']

        flags = []
        cc = None
        if isinstance(self.compiler, MSVCCompiler):
            # TODO: support Windows
            build_ext.build_extensions(self)
            return
        else:
            flags = ['-fopenmp']
            try:
                cc = self.compiler.compiler_so[0]
            except AttributeError:
                # we should only arrive here if the compiler is a BCPPCompiler
                print('compiler class does not contain a compiler_so object')
                build_ext.build_extensions(self)
                return

        if (can_compile_with_openmp(cc, flags, gomp_paths)):
            for e in self.extensions:
                e.extra_compile_args += flags
                e.library_dirs += gomp_paths
                # setuptools silently ignores extra_compile_args;
                # as a workaround, we include the flags here as well.
                # (see https://github.com/pypa/setuptools/issues/473)
                e.extra_link_args += flags + ['-lgomp']
        else:
            print('##################################################')
            print('# failed to compile with OpenMP;                 #')
            print('# some sections of SCS will not be parallelized. #')
            print('##################################################')
        build_ext.build_extensions(self)


def install_scs(**kwargs):
    blas_info = kwargs['blas_info']
    lapack_info = kwargs['lapack_info']

    extra_compile_args = ["-O3"]
    library_dirs = []
    extra_link_args = []
    libraries = []
    extra_define_macros = []
    sources = ['scsmodule.c', ] + glob(os.path.join(root_dir, 'src/*.c')) + glob(os.path.join(root_dir, 'linsys/*.c'))
    include_dirs = [root_dir, os.path.join(root_dir, 'include'), get_include(), os.path.join(root_dir, 'linsys')]
    define_macros = [('PYTHON', None), ('CTRLC', 1), ('COPYAMATRIX', None)]

    if system() == 'Linux':
        libraries += ['rt']
    if args.float32:
        define_macros += [('FLOAT', 1)] # single precision floating point
    if args.extraverbose:
        define_macros += [('EXTRAVERBOSE', 999)] # for debugging
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
                        sources=sources + glob(os.path.join(root_dir, 'linsys/direct/*.c')) + glob(os.path.join(root_dir, 'linsys/direct/external/*.c')),
                        define_macros=define_macros + extra_define_macros,
                        include_dirs=include_dirs + [os.path.join(root_dir, 'linsys/direct/'), os.path.join(root_dir, 'linsys/direct/external/')],
                        library_dirs=library_dirs,
                        libraries=libraries,
                        extra_link_args=extra_link_args,
                        extra_compile_args=extra_compile_args
                        )

    _scs_indirect = Extension(
                        name='_scs_indirect',
                        sources=sources + glob(os.path.join(root_dir, 'linsys/indirect/*.c')),
                        define_macros=define_macros + extra_define_macros + [('INDIRECT', None)],
                        include_dirs=include_dirs + [os.path.join(root_dir, 'linsys/indirect/')],
                        library_dirs=library_dirs,
                        libraries=libraries,
                        extra_link_args=extra_link_args,
                        extra_compile_args=extra_compile_args
                        )

    ext_modules = [_scs_direct, _scs_indirect]

    if args.gpu:
        _scs_gpu = Extension(
                        name='_scs_gpu',
                        sources=sources + glob(os.path.join(root_dir, 'linsys/gpu/*.c')),
                        define_macros=define_macros + [('GPU', None)],
                        include_dirs=include_dirs + [os.path.join(root_dir, 'linsys/gpu/'), '/usr/local/cuda/include'],
                        library_dirs=library_dirs + ['/usr/local/cuda/lib', '/usr/local/cuda/lib64'], # TODO probably not right for windows
                        libraries=libraries + ['cudart', 'cublas', 'cusparse'],
                        extra_link_args=extra_link_args,
                        extra_compile_args=extra_compile_args
                        )
        ext_modules += [_scs_gpu]

    setup(name='scs',
            version='1.2.6',
            author = 'Brendan O\'Donoghue',
            author_email = 'bodonoghue85@gmail.com',
            url = 'http://github.com/cvxgrp/scs',
            description='scs: splitting conic solver',
            py_modules=['scs'],
            ext_modules=ext_modules,
            cmdclass = {'build_ext': build_ext_scs},
            install_requires=["numpy >= 1.7","scipy >= 0.13.2"],
            license = "MIT",
            zip_safe=False,
            long_description="Solves convex cone programs via operator splitting. Can solve: linear programs (LPs), second-order cone programs (SOCPs), semidefinite programs (SDPs), exponential cone programs (ECPs), and power cone programs (PCPs), or problems with any combination of those cones. See http://github.com/cvxgrp/scs for more details."
            )

run_install()
