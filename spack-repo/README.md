# Installing ReMKiT1D Using Spack

_Last tested on Ubuntu 24.04 LTS with Spack v1.2.1._

The [Spack](https://github.com/spack/spack) package manager can be used to manage the dependencies of ReMKiT1D on Linux workstations, HPC systems, container builds, and future multi-code packages. The Spack integration is intended to:

- Manage ReMKiT1D dependencies consistently on Linux and HPC systems.
- Allow dependency upgrades and testing through modifications to `package.py`.
- Support developer builds from local source trees using `spack develop`.

**The Spack build system for ReMKiT1D is currently under development.** It is strongly recommended the user only proceeds if they are familiar with Spack, and if they will be installing ReMKiT1D in their own development environment.

## One-Time Spack Installation

### Install Spack

Install the prerequisites listed in the [Spack documentation](https://spack.readthedocs.io/en/v1.2.0/installing_prerequisites.html).

Next, clone Spack into a permanent location:

```bash
git clone --depth=2 https://github.com/spack/spack.git
cd spack

git fetch origin --tags --depth=2
git switch --detach v1.2.1
```

### Add Spack to your shell environment

Add the following line to your shell startup file (e.g. `.bash_profile` or `.bashrc`), editing the path to the `spack/` repo as appropriate:

```bash
. /path_to_installation/spack/share/spack/setup-env.sh
```

Verify the installation by calling the version number:

```bash
spack --version
```

### Register available compilers

If Fortran/C compilers already exist on the user's system, let Spack recognise them:

```bash
spack compiler find
```

### Register the ReMKiT1D package repository

From within the `ReMKiT1D/` repository, run the following to let Spack see ReMKiT1D.

```bash
spack repo add spack-repo
```

## Using the ReMKiT1D Spack Environment

Always begin by activating the environment:

```bash
spack env activate .
```

To install ReMKiT1D dependencies (Note: this may take one to tens of minutes depending on available resources):

```bash
spack install --only dependencies
```

Optionally, `spack install` can be parallelised by adding `-j <num. cores>`.

One dependency of ReMKiT1D is a version of CMake, which you can use while in the active Spack environment by running:

```bash
spack add cmake
```

Then, build and test ReMKiT1D as normal:

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build
```

## Uninstalling Packages

List installed packages:

```bash
spack find
```

Remove the ReMKiT1D package:

```bash
spack uninstall remkit1d
```

Remove all dependents:

```bash
spack uninstall --all --dependents
```

Remove everything from the active environment:

```bash
spack uninstall --all
```

Clean Spack caches and build stages:

```bash
spack clean -a
```

## Developer Instructions

### Repository Layout

In general, a local Spack repository should have the following structure:

```text
ReMKiT1D/
├── spack.yaml
└── spack-repo/
    ├── repo.yaml
    └── packages/
        └── remkit1d/
            └── package.py
```

Once these files exist (and are valid syntax), verify that Spack can see the package:

```bash
spack list remkit1d
```

Expected output:

```text
==> 1 packages
remkit1d
```

### Modifying package.py

Dependency specifications are maintained in:

```text
spack-repo/packages/remkit1d/package.py
```

After modifying the package definition, force a reconcretization:

```bash
spack concretize -f
```

### Testing Without Spack

To verify that ReMKiT1D still builds using the older Docker-based workflow:

```bash
docker run -d -it \
    --name remkit1d-build \
    --mount type=bind,source=/home/user/ReMKiT1D,target=/home/ReMKiT1D \
    remkit1d:latest
```

This helps distinguish Spack integration issues from upstream build-system issues.

## Full Spack Build and Test

The following one-liner performs the following (in parallel with 16 cores):

1. Reconcretizes the environment,
2. Installs dependencies,
3. Builds ReMKiT1D into a folder named `build-spack`,
4. Runs the unit test suite.

```bash
spack concretize -f -j 16 && \
spack install -j 16 --only dependencies && \
cmake -S . -B build-spack && \
cmake --build build-spack -j 16 && \
ctest --test-dir build-spack
```
