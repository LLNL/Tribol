# Tribol: Contact Interface Physics Library

High fidelity simulations modeling complex interactions of moving bodies require
specialized contact algorithms to enforce constraints between surfaces that come
into contact in order to prevent penetration and to compute the associated contact
response forces. Tribol aims to provide a unified interface for various contact
algorithms, specifically, contact detection and enforcement, and serve as a
common infrastructure enabling the research and development of advanced contact
algorithms.


## Quick Start Guide

### Clone the repository

```
git clone --recursive git@github.com:LLNL/Tribol.git
```

### Setup for Development

Development tools can optionally be installed through the Spack package manager.
The command to do this is
```
python3 scripts/uberenv/uberenv.py --project-json=scripts/spack/devtools.json --spack-config-dir=scripts/spack/configs/<platform> --prefix=../tribol_devtools
```
Please verify the `compilers.yaml` and `packages.yaml` files in
`scripts/spack/configs/<platform>` match your system configuration.  In
particular, local versions of `pyshroud` and `pysphinx` are assumed present.

### Installing dependencies

Tribol dependency installation is managed through uberenv, which invokes a local
instance of the spack package manager to install and manage dependencies.  To
install dependencies, run

```
python3 scripts/uberenv/uberenv.py --spack-config-dir=scripts/spack/configs/<platform> --prefix=../tribol_libs
```

See additional options by running

```
python3 scripts/uberenv/uberenv.py --help
```

Tribol is tested on three platforms: 
- Ubuntu 20.04 LTS (via Windows WSL 2)
- TOSS 3
- BlueOS

Note, running on Ubuntu 20.04 LTS requires inspecting and tailoring the
`compilers.yaml` and `packages.yaml` files for your specific system.

See `scripts/spack/packages/tribol/package.py` for possible variants in the
spack spec. The file `scripts/spack/specs.json` lists spack specs which are
known to build successfully on different platforms.  Note the development tools
can be built with dependencies using the `+devtools` variant.

### Build the code

After running uberenv, a host config file is created in the tribol repo root
directory.  Use the `config-build.py` script to create build and install
directories and invoke CMake.

```
python3 ./config-build.py -hc <host-config>
```

Enter the build directory and run

```
make -j
```

to build Tribol.


## Dependencies

The Tribol contact physics library requires:
- CMake 3.14 or higher
- C++11 compiler
- MPI
- mfem
- axom

## License

Tribol is distributed under the terms of the MIT license. All new contributions must be 
made under this license.

See [LICENSE](LICENSE) and [NOTICE](NOTICE) for details.

SPDX-License-Identifier: MIT

LLNL-CODE-846697

## SPDX usage

Individual files contain SPDX tags instead of the full license text.
This enables machine processing of license information based on the SPDX
License Identifiers that are available here: https://spdx.org/licenses/

Files that are licensed as MIT contain the following
text in the license header:

    SPDX-License-Identifier: (MIT)

## External Packages

Tribol bundles some of its external dependencies in its repository.  These
packages are covered by various permissive licenses.  A summary listing
follows.  See the license included with each package for full details.


[//]: # (Note: The spaces at the end of each line below add line breaks)

PackageName: BLT  
PackageHomePage: https://github.com/LLNL/blt  
PackageLicenseDeclared: BSD-3-Clause  

PackageName: uberenv  
PackageHomePage: https://github.com/LLNL/uberenv  
PackageLicenseDeclared: BSD-3-Clause  
