# Contributing

This file is under active development. Please send any suggestions to stefan.mijin@ukaea.uk.

- All feature development should be done in separate branches. Please use branch names of the format dev-\${userID}-\${featureName}, where \${userID} is best set to the main feature developer's GitHub username and \${featureName} is a short name for the feature.
- Raise a relevant issue (if it doesn't exist) before creating your feature branch, and refer to the issue in the PR
- All new features should include reasonable unit test coverage, especially if the feature is low level. Pull requests which fail CI tests will be automatically rejected. 
- Document all new features using FORD 
- Every PR must include corresponding changes in the CHANGELOG.md file
- Ideally, all new features should conform to the Fortran 2018 standard, and should not be using standards older than Fortran 2003 (where compilers allow).


- Coding style points: 
    - All variables and functions should be camelCase
    - Names of classes/types should be PascalCase
    - Names of files/modules/submodules should be snake_case
    - Separate declaration and implementation into modules and submodules, with files containing submodules ending in _procedures.f90. However, it is acceptable that minor modules not containing classes have their declaration and implementation all in one file.
    - Subroutine based constructors are preferred to function constructors (i.e. `call newObject%init()` over `newObject=Class()`)
    - Use modern logical operators (< over .lt. etc.)
    
