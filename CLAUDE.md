# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

SCS (Splitting Conic Solver) is a C library for solving large-scale convex cone optimization problems (LP, SOCP, SDP, ECP, PCP, complex SDP, spectral cones) using the ADMM splitting algorithm with optional Anderson acceleration.

## Build Commands

```bash
make                        # Build with BLAS/LAPACK (default)
make USE_LAPACK=0           # Build without LAPACK (no SDP support)
make DLONG=1                # Build with 64-bit integer indexing
make USE_SPECTRAL_CONES=1   # Build with spectral cone support (requires LAPACK)
make USE_OPENMP=1           # Build with OpenMP parallelization
make purge                  # Clean all build artifacts
```

CMake alternative:
```bash
mkdir build && cd build
cmake -DBUILD_TESTING=ON -DUSE_LAPACK=ON ..
cmake --build .
```

Outputs go to `out/`: static libs (`libscsdir.a`, `libscsindir.a`), shared libs, and demo/test binaries.

## Testing

```bash
make test                        # Build test binaries
./out/run_tests_direct           # Run direct solver tests
./out/run_tests_indirect         # Run indirect solver tests

# With options:
make test DLONG=1 USE_LAPACK=1
```

Tests use the `minunit.h` framework in `test/`. Individual test problems are in `test/problems/` and `test/spectral_cones_problems/`.

## Architecture

### Algorithm
SCS uses ADMM to split convex cone problems into two alternating steps:
1. **Linear system solve** — factorizes/solves a KKT-like system
2. **Cone projection** — projects onto the constraint cone

Anderson acceleration (`src/aa.c`) speeds convergence. Problem normalization/equilibration (`src/normalize.c`) improves conditioning before solving.

### Key Interfaces
- **`include/linsys.h`** — abstract interface that all linear system solvers must implement; `ScsLinSysWork` is the opaque solver-specific workspace
- **`include/cones.h`** — cone projection interface; `ScsConeWork` holds precomputed cone data
- **`include/scs.h`** — public API: `ScsMatrix`, `ScsData`, `ScsCone`, `ScsSettings`, `ScsSolution`

### Linear Solver Backends (in `linsys/`)
- `cpu/direct/` — QDLDL sparse Cholesky (default, bundled in `linsys/external/qdldl/`)
- `cpu/indirect/` — iterative (conjugate gradient / LSQR)
- `mkl/direct/` — Intel MKL Pardiso (enabled via `MKLROOT`)
- `accelerate/direct/` — Apple Accelerate sparse LDLt (macOS only, `make accelerate`)
- `cudss/direct/` — NVIDIA cuDSS GPU solver
- `gpu/indirect/` — GPU iterative solver

The Makefile selects which backend to compile based on flags. Only one direct or indirect backend is linked per build.

### Cone Types
Supported cones are registered in `src/cones.c`. The `ScsCone` struct fields map to: `z` (zero), `l` (non-negative), `q` (second-order), `s` (PSD), `cs` (complex PSD), `p` (power), `ep`/`ed` (exponential primal/dual). Spectral cones (`USE_SPECTRAL_CONES=1`) add log-determinant, nuclear norm, and sum-of-largest-eigenvalues support via `src/spectral_cones/`.

### Compiler Flags
Defined in `scs.mk`: `-Wall -Wwrite-strings -pedantic -funroll-loops -Werror=incompatible-pointer-types`. Type aliases `scs_float` and `scs_int` in `include/scs_types.h` change with `SFLOAT` and `DLONG` flags.

## Key Data Flow

The main solve lifecycle is:
1. **`scs_init`** — validates data, deep-copies and normalizes `A`/`P`/`b`/`c`, sets up linear system solver (factorization for direct), sets up cone workspace
2. **`scs_solve`** — ADMM loop:
   - Solve linear system (KKT) to get search direction (`project_lin_sys`)
   - Project onto dual cone via Moreau decomposition (`project_cones`)
   - Apply Anderson acceleration every `acceleration_interval` iterations
   - Check convergence / infeasibility / time limit
   - Optionally update scale parameter (adaptive scaling)
3. **`scs_finish`** — frees all workspace memory

`scs_update` allows changing `b`/`c` between solves without re-factorizing.

## Adding a New Linear Solver Backend

1. Create a new directory under `linsys/`, e.g. `linsys/mybackend/direct/`
2. Implement the 5 functions declared in `include/linsys.h`:
   - `scs_init_lin_sys_work()` — allocate workspace, factorize if direct
   - `scs_solve_lin_sys()` — solve the KKT system
   - `scs_update_lin_sys_diag_r()` — handle R diagonal updates
   - `scs_free_lin_sys_work()` — deallocate workspace
   - `scs_get_lin_sys_method()` — return a descriptive string
3. Create `private.h` defining `struct SCS_LIN_SYS_WORK` with your backend's fields
4. Add Makefile rules following the pattern of the existing backends (see `mkl`, `cudss`)
5. Add CMake rules following the same pattern in `CMakeLists.txt`
6. Add a CI workflow in `.github/workflows/`

See `linsys/cpu/direct/` for the reference implementation.

## Adding a New Cone Type

1. Add the new field(s) to `ScsCone` in `include/scs.h`
2. In `src/cones.c`:
   - Update `SCS(validate_cones)` to validate the new cone's parameters
   - Update `SCS(init_cone)` to compute cone boundaries and allocate workspace
   - Update `SCS(proj_dual_cone)` to project onto the new cone
   - Update `SCS(get_cone_header)` for printing
   - Update `SCS(set_r_y)` for the R diagonal entries
3. If the projection is complex, put it in a separate source file (like `src/exp_cone.c`)
4. Update the Makefile object lists and CMakeLists.txt source lists
5. Add a test in `test/problems/`

Cone ordering in the rows of `A` must match the order in `proj_dual_cone`.
