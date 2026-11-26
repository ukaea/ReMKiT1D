src_dir: ./src/modules/
output_dir: ./doc/ford
project: ReMKiT1D
project_github: https://github.com/ukaea/ReMKiT1D
version: 1.4.0
summary: Reactive Multifluid and Kinetic Transport in 1D
author: Stefan Mijin
github: https://github.com/SMijin
email: stefan.mijin@ukaea.uk
predocmark: >
docmark_alt: #
predocmark_alt: <
source: false
graph: false
search: true
license: by-nc
extra_filetypes: sh #
display: private
         protected
         public
extra_vartypes: Vec
                Mat
                KSP
                PC
                KSPType
                KSPNormType
                KSPConvergedReason
                PetscInt
                PetscCount
                PetscScalar
                PetscReal
                PetscBool
                PetscErrorCode
                PetscFortranAddr
                N_Vector
                SUNLinearSolver
                SUNMatrix
                MPI_Comm
exclude_dir: ./src/modules/testing_support
page_dir: ./docs-pages

## Brief description

ReMKiT1D is a framework for building 1D fluid and kinetic models with reactive multispecies capabilities, focused on Scrape-Off Layer plasmas.
