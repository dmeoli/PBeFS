# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.cmake import CMakePackage

from spack.package import *


class Pbefs(CMakePackage):
    """Parallel best-first search in C++20: the sequential and parallel A*
    of an implicit state-space graph and branch-and-bound of a tree with an
    incumbent, on C++ threads, OpenMP and MPI."""

    homepage = "https://github.com/dmeoli/PBeFS"
    url = "https://github.com/dmeoli/PBeFS/releases/download/v0.1.0/PBeFS-0.1.0.tar.gz"
    git = "https://github.com/dmeoli/PBeFS.git"

    maintainers("dmeoli")

    license("LGPL-3.0-only", checked_by="dmeoli")

    version("main", branch="main")
    version("0.1.0", sha256="0000000000000000000000000000000000000000000000000000000000000000")

    variant("shared", default=True, description="Build the shared library")
    variant("openmp", default=True, description="The target PBeFS::OpenMP")
    variant("mpi", default=False, description="The target PBeFS::MPI")

    depends_on("cxx", type="build")
    depends_on("cmake@3.21:", type="build")
    depends_on("ads")
    depends_on("mpi", when="+mpi")

    def cmake_args(self):
        # FastFlow, which no package provides, would be downloaded at
        # configure time: it is left out
        return [
            self.define_from_variant("BUILD_SHARED_LIBS", "shared"),
            self.define_from_variant("PBEFS_WITH_OPENMP", "openmp"),
            self.define_from_variant("PBEFS_WITH_MPI", "mpi"),
            self.define("PBEFS_WITH_FASTFLOW", False),
            self.define("PBEFS_BUILD_EXAMPLES", False),
            self.define("PBEFS_BUILD_TOOLS", False),
            self.define("PBEFS_BUILD_TESTS", self.run_tests),
            self.define("BUILD_TESTING", self.run_tests),
        ]
