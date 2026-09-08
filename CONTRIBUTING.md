# Contributing to SCS

Thank you for your interest in contributing to SCS! For detailed guidelines,
see the [contribution docs](https://www.cvxgrp.org/scs/contributing).

## Quick Start

```bash
# Build
make

# Build and run tests
make test
./out/run_tests_direct
./out/run_tests_indirect

# Build with options
make DLONG=1 USE_LAPACK=1

# Clean everything
make purge
```

> **The spectral-cone tests are skipped by default.** `USE_SPECTRAL_CONES`
> defaults to `0` (see `scs.mk`), so a default `make test` reports nine tests
> as `skipped` and exercises none of `src/spectral_cones/` — while still
> printing `ALL TESTS PASSED`. If you touch the log-determinant, nuclear-norm,
> sum-of-largest or ell1 cone code, build and test with the flag enabled:
>
> ```bash
> make test USE_SPECTRAL_CONES=1
> ./out/run_tests_direct
> ```

## Development Workflow

1. Fork the repo and create a feature branch from `master`
2. Make your changes
3. Run the test suite to ensure nothing is broken
4. Submit a pull request against `master`

## Code Style

- C99 with `-Wall -Wwrite-strings -pedantic -Werror=incompatible-pointer-types`
- Internal functions use the `SCS(name)` macro (expands to `_scs_name`)
- Public API functions use the `scs_` prefix directly
- Types: `scs_float` (double/float), `scs_int` (int/long long)

## Project Layout

| Directory | Contents |
|-----------|----------|
| `include/` | Public API and internal headers |
| `src/` | Core solver implementation |
| `linsys/` | Linear solver backends (pluggable) |
| `test/` | Test suite (minunit framework) |
| `docs/src/` | Sphinx documentation source |

## Where to Read Next

| Topic | Reference |
|-------|-----------|
| Algorithm, termination criteria, scaling | [`docs/src/algorithm/`](docs/src/algorithm/) ([online](https://www.cvxgrp.org/scs/algorithm/)) |
| Linear solver backends and the KKT system | [`docs/src/linear_solver/`](docs/src/linear_solver/) ([online](https://www.cvxgrp.org/scs/linear_solver/)) |
| **Adding a new linear solver backend** | [Implementing a new linear solver](https://www.cvxgrp.org/scs/linear_solver/#implementing-a-new-linear-solver) — implement `ScsLinSysWork` and the functions in `include/linsys.h`; see `linsys/` for seven worked examples |
| Supported cones and their data layout | [`docs/src/api/cones.rst`](docs/src/api/cones.rst) ([online](https://www.cvxgrp.org/scs/api/cones.html)) |
| Compile-time flags | [`docs/src/api/compile_flags.rst`](docs/src/api/compile_flags.rst) ([online](https://www.cvxgrp.org/scs/api/compile_flags.html)) |
