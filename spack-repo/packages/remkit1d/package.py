# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.cmake import CMakePackage

from spack.package import depends_on, license, variant, version


class Remkit1d(CMakePackage):
    """Framework for 1D multifluid and kinetic simulations geared towards Scrape-Off Layer plasmas"""

    homepage = "https://github.com/ukaea/ReMKiT1D"
    git = "https://github.com/ukaea/ReMKiT1D.git"
    supplier = "UK Atomic Energy Authority"

    license("GPL-3.0-or-later")

    version(
        "v1.3.0",
        tag="v1.3.0",
    )

    depends_on("cmake@3.18:", type="build")
    depends_on("c", type="build")
    depends_on("fortran")
    depends_on("mpi")
    depends_on("hdf5+fortran+mpi+hl")
    depends_on("petsc@3.17.5+fortran+mpi+hypre~debug")
    depends_on("hypre+fortran")
    depends_on("sundials@7.7.0+CVODE+mpi+f2003+lapack+shared")
    depends_on("json-fortran")

    variant("tests", default=False)
    depends_on("pfunit", when="+tests")

    def install(self, spec, prefix):
        pass

    def cmake_args(self):
        # FIXME: Add arguments other than
        # FIXME: CMAKE_INSTALL_PREFIX and CMAKE_BUILD_TYPE
        # FIXME: If not needed delete this function
        args = []
        return args
