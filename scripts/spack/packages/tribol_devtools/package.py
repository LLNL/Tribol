# Copyright 2013-2019 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *

class TribolDevtools(BundlePackage):
    """This is a set of tools necessary for the developers of Tribol"""

    version('fakeversion')

    depends_on("doxygen")
    depends_on("python")
    depends_on("py-shroud")
    depends_on("py-sphinx")
    depends_on("llvm+clang@10.0.0")
