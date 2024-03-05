# CHANGELOG

## v1.2.0, 2024-03-05

- Bug fixes

### New Features

-
### Bug Fixes

- Fixed bugs with unary contract and expand. Added tests for non-trivial unary operators. 

## v1.1.0, 2024-02-02

- Solver and integrator improvements
- Debug mode bug fixes
- Multilinear interpolation support
- Support for some CRM features with non-default electron species ID
- Default support for release build assertions during startup

### New Features

- The BDE integrator internal controller now attempts to reduce the number of substeps back to 1 after every 50 successful integrations
- The BDE integrator now also has a hard maximum of attempted restarts regardless of whether consolidation happens. this avoids crashes due to allocating too much memory for the timestep size buffer
- The BDE integrator now prints the approximate convergence bottleneck variable (feature not completely reliable)
- Command line PETSc support for setting up the KSP solver object has been implemented. The command line ksp_type takes precedent over the integrator options. This now allows for command line customization of the KSP object. Note that tolerances are still set from config options. 
- New derivation type nDLinInterpDerivation - linear interpolation on n-dimensional data
- ModelboundCRMData can now have an associated non-default electron species ID (useful for multiple flux tube models). Note that Janev transitions still have hardcoded values. 
- Added the option to reset the value of the time variable upon loading from restart files
- Assertions are now run in release mode during simulation startup

### Bug Fixes

- Fixed a number of bugs in the code where divide-by-zero FPEs would be raised in debug mode. Some still remain.
- Fixed a bug with gfortran-11.4 which caused bound changes in an MPI buffer
- Some tests were failing when run on an M2 Mac due to FPE differences. This has been fixed. 

## v1.0.0, 2023-06-21

- Initial release

### Breaking Changes

- N/A

### Deprecations

- N/A

### New Features

- N/A

### Bug Fixes

- N/A
