# CHANGELOG

## v1.1.0, 2023-08-31

- Solver and integrator improvements
- Debug mode bug fixes
- Multilinear interpolation support

### New Features

- The BDE integrator internal controller now attempts to reduce the number of substeps back to 1 after every 50 successful integrations
- Command line PETSc support for setting up the KSP solver object has been implemented. The command line ksp_type takes precedent over the integrator options. This now allows for command line customization of the KSP object. Note that tolerances are still set from config options. 
- New derivation type nDLinInterpDerivation - linear interpolation on n-dimensional data

### Bug Fixes

- Fixed a number of bugs in the code where divide-by-zero FPEs would be raised in debug mode 

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