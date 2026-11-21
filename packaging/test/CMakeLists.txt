# --------------------------------------------------------------------------- #
#    CMake file for PBeFS/packaging/test                                      #
#                                                                             #
#    A project that uses an installed PBeFS through find_package(PBeFS), and, #
#    when pkg-config is there, also through pbefs.pc, and the symbolic        #
#    searches through PBeFS::Symbolic if the package has them: the recipes of #
#    packaging/ build it against the package they made and run it.            #
#                                                                             #
#        $ cmake -S packaging/test -B build-use -DCMAKE_PREFIX_PATH=<prefix>  #
#        $ cmake --build build-use && ctest --test-dir build-use              #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Università di Pisa                              #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.21)
project(PBeFSUse LANGUAGES CXX)

find_package(PBeFS 0.1 REQUIRED)
add_executable(use-cmake main.cpp)
target_link_libraries(use-cmake PRIVATE PBeFS::PBeFS)

enable_testing()
add_test(NAME use-cmake COMMAND use-cmake)

if (TARGET PBeFS::Symbolic)
    add_executable(use-symbolic symbolic.cpp)
    target_link_libraries(use-symbolic PRIVATE PBeFS::Symbolic)
    add_test(NAME use-symbolic COMMAND use-symbolic)
endif ()

# pkg-config gives the flags of the library, not the standard of C++
find_package(PkgConfig QUIET)
if (PKG_CONFIG_FOUND)
    pkg_check_modules(PBEFSPC REQUIRED IMPORTED_TARGET pbefs)
    add_executable(use-pkgconfig main.cpp)
    target_link_libraries(use-pkgconfig PRIVATE PkgConfig::PBEFSPC)
    target_compile_features(use-pkgconfig PRIVATE cxx_std_20)
    add_test(NAME use-pkgconfig COMMAND use-pkgconfig)
endif ()

# --------------------------------------------------------------------------- #
