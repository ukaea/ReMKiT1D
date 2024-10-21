# CHANGELOG

## v1.2.1, 2024-10-21

- Added an option to remove the logLei discontinuity in NRL formulary fit
- Added more output to BDE integrator
- Added more unary transforms
- Added key for controlling BDE integrator consolidation

### Breaking Changes

- N/A

### New Features

- New standard textbook option to remoge logLei discontinuity in NRL formularly fit by moving branch cutoff to e^2 * Z^2 eV
- BDE integrator now outputs the current number of substeps and progress through those
- New "step" unary transform. 1 if the node evaluates > 0, 0 otherwise.
- New "filterWithin" unary transform. 0 if the value of the node is not within the first two passed realParams, otherwise no effect.
- New "BDEConsolidationInterval" key to set the number of timesteps for the integrator to attempt consolidation after. Defaults to 50 (old hardcoded value).

### Bug Fixes

- N/A

## v1.2.0, 2024-09-19

- Added CVODE integrator as an option
- Added new derivation-based explicit term
- Added new unary transformations
- New variable and manipulator features
- New timeloop mode - output-driven timeloop
- New restart option - initial output index
- Support for under- and over-relaxation in BDE integrator
- Updated documentation
- Bug fixes

### Breaking Changes

- Due to the need to now differentiate between term types in JSON input, pre v1.2.0 config files will not work with 1.2.0
- New explicit term interface incompatible with the previous one (relevant only for existing tests)

### New Features

- Basic CVODE integrator added
- Derivation-based explicit term added. This takes a derivation and an optional modelbound variable and evaluates to the product of the derivation result and the variable
- Added slope limiter related unary transformations 
- Variables can now be copied into/from arrays by passing a list of variable names to the container
- Manipulators now called before first time step
- The term evaluator manipulator can now accumulate values into the evaluation variable instead of overwriting it
- The term evaluator can now also explicitly request model and term updates (less fine-grained control than integrators)
- Variables can now be zeroed with a passed list of names
- Added option to explicitly change the maximum number of BDE integrator restarts (still hard-capped to 10)
- Another timeloop mode has been added where the output points are set, and the code makes sure they coincide with integrator steps. If the output point is sufficiently far away, the standard timestep behaviour is recovered. 
- It is now possible to set the initial output index. This is useful when restarting, allowing the user to avoid overwriting previous output files.
- New unary transform for flooring variables
- Timeloops now state which output index is written to
- BDE integrators can now have under- and over-relaxation set through a relaxation weight argument. A value <1 is under-relaxation, and >1 is over-relaxation.

### Bug Fixes

- Fixed bugs with unary contract and expand. Added tests for non-trivial unary operators. 
- Fixed segfault on finalize when no PETSc obj used
- Fixed weird segfault caused by having 0 implicit terms and adding a second general term to a group
- Fixed missing communication call in RK integrator
- Fixed wrong order of arguments in electron energy source term generator

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
