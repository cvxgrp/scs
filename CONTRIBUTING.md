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

See `CLAUDE.md` for architecture details and guides on adding new solver
backends or cone types.
