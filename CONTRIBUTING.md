# Contributing

Thank you for your interest in contributing to ReMKiT1D. Please direct questions or suggestions to stefan.mijin@ukaea.uk or open an issue on GitHub.

## Workflow

- Perform feature development on separate branches named `dev-${userID}-${featureName}`, where `${userID}` is the primary developer's GitHub username and `${featureName}` is a short descriptive name.
- Open or reference a relevant GitHub issue before creating a feature branch, and link the issue in the pull request.
- Include unit tests (pFUnit) for all new functionality, especially low-level components. Pull requests that fail CI checks will not be merged.
- Document new public modules, types, and procedures using FORD docstrings.
- Add an entry to `CHANGELOG.md` for every user-facing change.
- Conform to the Fortran 2018 standard wherever possible, and avoid pre-Fortran 2003 constructs.

## Local checks before opening a PR

```bash
# Build and test in DEBUG mode (enables assertions and runtime checks)
cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=DEBUG
cmake --build build-debug -j"$(nproc)"
ctest --test-dir build-debug --output-on-failure

# Verify optimised RELEASE build
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=RELEASE
cmake --build build-release -j"$(nproc)"
ctest --test-dir build-release --output-on-failure

# Optionally verify coverage in TESTING mode and FORD docs generation
cmake -S . -B build-testing -DCMAKE_BUILD_TYPE=TESTING
cmake --build build-testing -j"$(nproc)"
ctest --test-dir build-testing --output-on-failure
ford ford-docs.md
```

## Coding style

- Variables and procedures: `camelCase`
- Derived types / classes: `PascalCase`
- Files, modules, and submodules: `snake_case`
- Separate declaration and implementation into modules and submodules, naming submodule files `*_procedures.f90` (small utility modules without derived types may remain in a single file).
- Prefer subroutine-based constructors (`call newObject%init(...)`) over function constructors (`newObject = Class(...)`).
- Use modern relational operators (`<`, `<=`, `==`, `/=`, `>`, `>=` rather than `.lt.`, `.eq.`, etc.).
