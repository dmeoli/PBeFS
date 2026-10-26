import os

from conan import ConanFile
from conan.tools.build import check_min_cppstd
from conan.tools.cmake import CMake, CMakeToolchain, cmake_layout
from conan.tools.files import copy, get, rmdir

required_conan_version = ">=2.0"


class PbefsConan(ConanFile):
    name = "pbefs"
    description = ("Parallel best-first search in C++20: A* and "
                   "branch-and-bound on threads, OpenMP and MPI")
    license = "LGPL-3.0-only"
    url = "https://github.com/conan-io/conan-center-index"
    homepage = "https://github.com/dmeoli/PBeFS"
    topics = ("search", "a-star", "branch-and-bound", "parallel")
    package_type = "library"
    settings = "os", "arch", "compiler", "build_type"
    options = {"shared": [True, False], "fPIC": [True, False]}
    default_options = {"shared": False, "fPIC": True}

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def configure(self):
        if self.options.shared:
            self.options.rm_safe("fPIC")

    def layout(self):
        cmake_layout(self, src_folder="src")

    def validate(self):
        check_min_cppstd(self, 20)

    def source(self):
        get(self, **self.conan_data["sources"][self.version], strip_root=True)

    def generate(self):
        # the library alone: the targets of the runtimes are those of the
        # CMake package of the project, which CMakeDeps does not use
        tc = CMakeToolchain(self)
        tc.cache_variables["PBEFS_WITH_OPENMP"] = False
        tc.cache_variables["PBEFS_WITH_MPI"] = False
        tc.cache_variables["PBEFS_WITH_FASTFLOW"] = False
        tc.cache_variables["PBEFS_BUILD_TESTS"] = False
        tc.cache_variables["PBEFS_BUILD_EXAMPLES"] = False
        tc.cache_variables["PBEFS_BUILD_TOOLS"] = False
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        copy(self, "LICENSE", self.source_folder,
             os.path.join(self.package_folder, "licenses"))
        cmake = CMake(self)
        cmake.install()
        rmdir(self, os.path.join(self.package_folder, "lib", "cmake"))
        rmdir(self, os.path.join(self.package_folder, "lib", "pkgconfig"))
        rmdir(self, os.path.join(self.package_folder, "share"))

    def package_info(self):
        self.cpp_info.libs = ["PBeFS"]
        self.cpp_info.set_property("cmake_file_name", "PBeFS")
        self.cpp_info.set_property("cmake_target_name", "PBeFS::PBeFS")
        self.cpp_info.set_property("pkg_config_name", "pbefs")
        if self.settings.os in ("Linux", "FreeBSD"):
            self.cpp_info.system_libs = ["pthread", "m"]
