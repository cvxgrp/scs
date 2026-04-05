<h1 align="center" margin=0px>
<img src="https://github.com/cvxgrp/scs/blob/master/docs/src/_static/scs_logo.png" alt="Intersection of a cone and a polyhedron" width="450">
</h1>

[![Build Status](https://github.com/cvxgrp/scs/actions/workflows/build.yml/badge.svg)](https://github.com/cvxgrp/scs/actions/workflows/build.yml)
[![Documentation](https://img.shields.io/badge/docs-online-brightgreen?logo=read-the-docs&style=flat)](https://www.cvxgrp.org/scs/)
[![PyPI Downloads](https://img.shields.io/pypi/dm/scs.svg?label=PyPI%20downloads&cacheSeconds=86400)](https://pypistats.org/packages/scs)
[![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/scs.svg?label=Conda%20downloads)](https://anaconda.org/conda-forge/scs)

SCS (`splitting conic solver`) is a numerical optimization package for solving
large-scale convex cone problems. The current version is `3.2.11`.

The full documentation is available [here](https://www.cvxgrp.org/scs/).

If you wish to cite SCS please cite the papers listed [here](https://www.cvxgrp.org/scs/citing).

## Overview

SCS solves convex cone programs of the form:

```
minimize    (1/2) x'Px + c'x
subject to  Ax + s = b, s in K
```

where `K` is a product of convex cones. It uses an operator-splitting method
(ADMM) with Anderson acceleration to achieve fast convergence.

**Supported cones:** linear (LP), second-order (SOCP), semidefinite (SDP),
exponential, power, complex semidefinite, and spectral cones (log-determinant,
nuclear norm, sum-of-largest eigenvalues).

## Features

- Multiple linear solver backends: sparse direct (QDLDL), iterative (CG),
  Intel MKL Pardiso, NVIDIA cuDSS, GPU
- Anderson acceleration for faster convergence
- Problem data normalization/equilibration for numerical stability
- Warm-starting and incremental `b`/`c` updates via `scs_update`
- Ctrl-C signal handling for graceful interruption
- Configurable precision (`float`/`double`) and index type (`int`/`long long`)
- MIT licensed

## Building from Source

### Make (default)

```bash
make                    # Build direct + indirect solvers with BLAS/LAPACK
make test               # Build test binaries
./out/run_tests_direct  # Run direct solver tests
./out/run_tests_indirect
```

### CMake

```bash
mkdir build && cd build
cmake -DBUILD_TESTING=ON -DUSE_LAPACK=ON ..
cmake --build .
ctest
```

### Common Build Flags

| Flag | Default | Description |
|------|---------|-------------|
| `USE_LAPACK` | `1` | Enable BLAS/LAPACK (required for SDP support) |
| `DLONG` | `0` | Use 64-bit integer indexing |
| `SFLOAT` | `0` | Use single-precision floats |
| `USE_SPECTRAL_CONES` | `0` | Enable spectral cone support (requires LAPACK) |
| `USE_OPENMP` | `0` | Enable OpenMP parallelization |

Example: `make DLONG=1 USE_LAPACK=1`

## Language Interfaces

SCS has interfaces for several languages:

- **Python:** `pip install scs` ([PyPI](https://pypi.org/project/scs/))
- **Julia:** [SCS.jl](https://github.com/jump-dev/SCS.jl)
- **R:** [scs](https://cran.r-project.org/package=scs)
- **MATLAB:** See the [documentation](https://www.cvxgrp.org/scs/install/matlab.html)
- **Ruby:** [scs-ruby](https://github.com/ankane/scs-ruby)

SCS is the default solver in [CVXPY](https://www.cvxpy.org/).

## Project Structure

```
include/        Public API and internal headers
src/            Core solver implementation
linsys/         Linear solver backends (pluggable architecture)
  cpu/direct/     Sparse Cholesky (QDLDL, default)
  cpu/indirect/   Conjugate gradient
  mkl/direct/     Intel MKL Pardiso
  cudss/direct/   NVIDIA cuDSS
  gpu/indirect/   GPU iterative solver
  external/       Vendored dependencies (AMD, QDLDL)
test/           Test suite (minunit framework)
docs/src/       Sphinx documentation source
```

## License

MIT License. See [LICENSE.txt](LICENSE.txt).
