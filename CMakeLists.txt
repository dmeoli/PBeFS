# --------------------------------------------------------------------------- #
#    Main CMake file for PBeFS                                                #
#                                                                             #
#    This file builds the library of PBeFS (the searches, templates           #
#    over the Domain, are in its headers), its interface targets for the      #
#    parallel runtimes, the examples, the benchmark tools and the tests, and  #
#    installs them with a package configuration file, so that other CMake     #
#    projects can use the library with find_package(PBeFS), and with the      #
#    pkg-config file pbefs.pc for the others.                                 #
#                                                                             #
#        $ cmake -S <source-path> -B <build-path>                             #
#        $ cmake --build <build-path>                                         #
#        $ ctest --test-dir <build-path>                                      #
#        $ cmake --build <build-path> --target install                        #
#                                                                             #
#    ADS, whose executors run the workers of the parallel searches, is        #
#    found by find_package(ADS) or else downloaded at a fixed commit. The     #
#    other parallel runtimes are optional: OpenMP and MPI are used if found,  #
#    FastFlow if its headers are in FASTFLOW_ROOT (/opt/fastflow by default)  #
#    or on the include path, else it is downloaded at a fixed commit; each    #
#    can be turned off with the PBEFS_WITH_* options. PBEFS_STRUCTURES        #
#    chooses the data structures of A*: those of ADS (ads, the default), of   #
#    the standard library (std) or of Boost (boost). Sylvan, the library of   #
#    binary decision diagrams of the symbolic searches, is optional too: it   #
#    is found by find_package(sylvan), or else downloaded at a fixed commit   #
#    (with CMake 3.23 or later, which it needs), and without it the library   #
#    PBeFS::Symbolic and its tests are left out.                              #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Università di Pisa                              #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.21)
cmake_policy(VERSION 3.15)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/DeriveVersion.cmake)
pbefs_derive_version(PBEFS_VERSION)

project(PBeFS
        VERSION ${PBEFS_VERSION}
        DESCRIPTION "Parallel best-first search: A* and branch-and-bound"
        HOMEPAGE_URL https://github.com/dmeoli/PBeFS
        LANGUAGES CXX)

# These variables make the code harder to read but easier to change.
set(modName ${PROJECT_NAME})
set(modNamespace ${PROJECT_NAME})

# Find out if it's being called by another project.
get_directory_property(hasParent PARENT_DIRECTORY)

# This adds the cmake directory to the module search paths,
# allowing us to use our modules.
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# ----- Settings ------------------------------------------------------------ #
# The tests, the examples and the tools are built only when PBeFS is
# the top-level project, not when another project adds it with
# add_subdirectory() or FetchContent.
option(PBEFS_BUILD_TESTS "Build the tests" ${PROJECT_IS_TOP_LEVEL})
option(PBEFS_BUILD_EXAMPLES "Build the examples" ${PROJECT_IS_TOP_LEVEL})
option(PBEFS_BUILD_TOOLS "Build the benchmark tools" ${PROJECT_IS_TOP_LEVEL})

# The options of the compiler for the tests and the tools, such as the
# warnings, which depend on the compiler and are therefore given by the
# presets of CMakePresets.json rather than here.
set(PBEFS_CXX_OPTIONS "" CACHE STRING
        "Options of the compiler for the tests and the tools (a list)")

# This creates the BUILD_TESTING option (default: ON) and enables the
# testing using the ctest executable.
# See: https://cmake.org/cmake/help/latest/manual/ctest.1.html
if (PBEFS_BUILD_TESTS)
    include(CTest)
endif ()

# Sets the default build type (if none was specified).
# See: https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html
include(BuildType)

# The standard installation directories, used by the subdirectories too.
include(GNUInstallDirs)

option(PBEFS_WITH_OPENMP "Build the OpenMP runtime if OpenMP is found" ON)
option(PBEFS_WITH_MPI "Build the MPI search if MPI is found" ON)
option(PBEFS_WITH_FASTFLOW "Build the FastFlow runtime if found" ON)

option(PBEFS_WITH_SYLVAN "Build the symbolic searches if Sylvan is found"
       ON)

option(PBEFS_FETCH_FASTFLOW "Download FastFlow if it is not installed" ON)
option(PBEFS_FETCH_SYLVAN "Download Sylvan if it is not installed" ON)
option(PBEFS_WITH_HIGHWAY
       "Prefetch through Google Highway if it is found, in a downloaded ADS"
       ON)

set(FASTFLOW_ROOT /opt/fastflow CACHE PATH
    "Directory holding the ff/ headers of FastFlow")

# a commit of a dependency to download, in the cache: one given by hand
# stays, while the default of an earlier version of PBeFS is replaced by
# that of this one, so that a build directory follows the new pins
function(pbefs_pin var value doc)
    if (NOT DEFINED CACHE{${var}} OR
        "$CACHE{${var}}" STREQUAL "$CACHE{${var}_DEFAULT}")
        set(${var} ${value} CACHE STRING "${doc}" FORCE)
    endif ()
    set(${var}_DEFAULT ${value} CACHE INTERNAL "The default of ${var}")
endfunction()

# the commit of FastFlow downloaded when it is not installed, the one the
# CI tests
pbefs_pin(PBEFS_FASTFLOW_COMMIT d476f66ab924d8d122f54b4b90aee00ef979aea8
          "The commit of FastFlow to download")

# the data structures of A* [see include/PBeFS/astar/Structures.h], and the
# commit of ADS downloaded when it is not installed
set(PBEFS_STRUCTURES ads CACHE STRING
    "The data structures of A*: ads, std or boost")
set_property(CACHE PBEFS_STRUCTURES PROPERTY STRINGS ads std boost)
set(ADS_ROOT /opt/ADS CACHE PATH "Prefix where ADS is installed")
pbefs_pin(PBEFS_ADS_COMMIT d826038dd17f5b79058f8b895433adae25a64bfd
          "The commit of ADS to download")

# the commit of Sylvan downloaded when it is not installed (the release
# 1.12.6), which downloads the release of Lace it needs
pbefs_pin(PBEFS_SYLVAN_COMMIT b5e236c39ef79a0b481dfdfb6b215745e3f7b852
          "The commit of Sylvan to download")

# ----- Requirements -------------------------------------------------------- #
find_package(Threads REQUIRED)

# ADS is installed (ADS_ROOT is a hint), or else the sources of the commit
# PBEFS_ADS_COMMIT are added as a subproject, with the backends and the
# interoperations that PBeFS does not use off, since PBeFS asks for OpenMP,
# FastFlow and MPI by itself, Highway, which gives ADS::prefetch, as
# PBEFS_WITH_HIGHWAY says (without it a prefetch does nothing), and Boost
# (header-only here, i.e., Boost::headers) and hwloc where they are found.
# With any ADS, CLOSED is the PointerTable of ADS, which the searches
# prefetch, and not the boost::unordered_flat_set of an ADS with Boost
# [see the README, "The data structures"].
set(PBEFS_ADS_DEFINITIONS ADS_POINTERTABLE_BOOST=0)
find_package(ADS CONFIG QUIET HINTS ${ADS_ROOT})
if (NOT ADS_FOUND)
    foreach (backend GMP BLAS MPI LEMON PYTHON STDPAR STDEXEC TASKFLOW KOKKOS
                     BOOST_MPI CEREAL EIGEN METIS KAHIP GRAPHBLAS NETWORKIT
                     CROARING SDSL)
        set(ADS_WITH_${backend} OFF CACHE STRING
            "The backend ${backend} of ADS")
    endforeach ()
    set(ADS_WITH_BOOST AUTO CACHE STRING "The backend BOOST of ADS")
    set(ADS_WITH_HWLOC AUTO CACHE STRING "The backend HWLOC of ADS")
    set(highway OFF)
    if (PBEFS_WITH_HIGHWAY)
        set(highway AUTO)
    endif ()
    set(ADS_WITH_HIGHWAY ${highway} CACHE STRING
        "The backend HIGHWAY of ADS, from PBEFS_WITH_HIGHWAY" FORCE)
    include(FetchContent)
    FetchContent_Declare(
            ADS
            GIT_REPOSITORY https://github.com/dmeoli/ADS.git
            GIT_TAG ${PBEFS_ADS_COMMIT})
    FetchContent_MakeAvailable(ADS)
endif ()
if (NOT TARGET ADS::ADS)
    message(FATAL_ERROR "PBeFS needs ADS: install it or set ADS_ROOT")
endif ()
if (PBEFS_STRUCTURES STREQUAL "boost")
    find_package(Boost 1.81 REQUIRED)
elseif (NOT PBEFS_STRUCTURES MATCHES "^(ads|std)$")
    message(FATAL_ERROR "PBEFS_STRUCTURES is ads, std or boost, not "
            "${PBEFS_STRUCTURES}")
endif ()

if (PBEFS_WITH_OPENMP)
    find_package(OpenMP COMPONENTS CXX)
endif ()

if (PBEFS_WITH_MPI)
    find_package(MPI COMPONENTS CXX)
endif ()

# FastFlow is header-only: its headers are looked for in FASTFLOW_ROOT and
# on the include path, and if they are not there the sources of the commit
# PBEFS_FASTFLOW_COMMIT are downloaded (its CMakeLists.txt is not used).
set(PBEFS_HAS_FASTFLOW OFF)
set(runtimes "")
if (PBEFS_WITH_FASTFLOW)
    find_path(FASTFLOW_INCLUDE_DIR ff/ff.hpp
              HINTS ${FASTFLOW_ROOT} PATH_SUFFIXES include)
    if (NOT FASTFLOW_INCLUDE_DIR AND PBEFS_FETCH_FASTFLOW)
        include(FetchContent)
        FetchContent_Declare(
                fastflow
                GIT_REPOSITORY https://github.com/fastflow/fastflow.git
                GIT_TAG ${PBEFS_FASTFLOW_COMMIT}
                SOURCE_SUBDIR none)
        FetchContent_MakeAvailable(fastflow)
        set(FASTFLOW_INCLUDE_DIR ${fastflow_SOURCE_DIR} CACHE PATH
            "Directory holding the ff/ headers of FastFlow" FORCE)
    endif ()
    if (FASTFLOW_INCLUDE_DIR AND EXISTS ${FASTFLOW_INCLUDE_DIR}/ff/ff.hpp)
        set(PBEFS_HAS_FASTFLOW ON)
    endif ()
endif ()

# Sylvan is installed, or else the sources of the commit PBEFS_SYLVAN_COMMIT
# are added as a subproject, in C, without its tests, examples and C++
# bindings, position independent so that a shared PBeFS::Symbolic can link
# it; the sanitizers of the tests reach it through CMAKE_C_FLAGS.
set(PBEFS_HAS_SYLVAN OFF)
if (PBEFS_WITH_SYLVAN)
    find_package(sylvan CONFIG QUIET)
    if (NOT TARGET sylvan::sylvan AND PBEFS_FETCH_SYLVAN)
        if (CMAKE_VERSION VERSION_LESS 3.23)
            message(STATUS "Sylvan needs CMake 3.23: the symbolic searches "
                    "are left out")
        else ()
            enable_language(C)
            set(SYLVAN_BUILD_TESTS OFF CACHE BOOL "Build Sylvan tests")
            set(SYLVAN_BUILD_EXAMPLES OFF CACHE BOOL "Build Sylvan examples")
            set(SYLVAN_BUILD_CPP OFF CACHE BOOL "Build Sylvan C++ bindings")
            set(SYLVAN_ENABLE_PIC ON CACHE BOOL "Build Sylvan with -fPIC")
            include(FetchContent)
            FetchContent_Declare(
                    sylvan
                    GIT_REPOSITORY https://github.com/trolando/sylvan.git
                    GIT_TAG ${PBEFS_SYLVAN_COMMIT})
            FetchContent_MakeAvailable(sylvan)
        endif ()
    endif ()
    if (TARGET sylvan::sylvan)
        set(PBEFS_HAS_SYLVAN ON)
    endif ()
endif ()

# ----- Configuration header ------------------------------------------------ #
# This will generate the PBeFSConfig.h header in the build directory.
string(MAKE_C_IDENTIFIER "${modName}" modMacro)
configure_file(cmake/${modName}Config.h.in ${modName}Config.h)

# ----- Library ------------------------------------------------------------- #
# The searches are templates over the Domain and live in the headers; the
# library holds what does not depend on it. With the BUILD_SHARED_LIBS
# variable we can specify if the library will be STATIC or SHARED, and one
# more INTERFACE target for each optional runtime carries what that runtime
# needs.
add_library(${modName})
target_compile_features(${modName} PUBLIC cxx_std_20)

# The name of the shared library carries the version of its binary interface,
# which any minor release may change while the major version is 0, and then
# any major release; MSVC exports the symbols only if asked to. The static
# library is position independent too, so that a shared library of another
# project can link it.
if (PROJECT_VERSION_MAJOR EQUAL 0)
    set(soVersion ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
else ()
    set(soVersion ${PROJECT_VERSION_MAJOR})
endif ()
set_target_properties(${modName} PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION ${soVersion}
        WINDOWS_EXPORT_ALL_SYMBOLS ON
        POSITION_INDEPENDENT_CODE ON)

# Note: do not GLOB files here.
target_sources(${modName} PRIVATE
        src/Grid.cpp
        src/Planning.cpp
        src/SASHeuristics.cpp
        src/SASTask.cpp
        src/Solution.cpp)

target_include_directories(
        ${modName} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:include>)
target_link_libraries(${modName} PUBLIC Threads::Threads ADS::ADS)
target_compile_definitions(${modName} PUBLIC ${PBEFS_ADS_DEFINITIONS})
if (NOT PBEFS_STRUCTURES STREQUAL "ads")
    string(TOUPPER ${PBEFS_STRUCTURES} structures)
    target_compile_definitions(${modName} PUBLIC
            PBEFS_STRUCTURES_${structures})
endif ()
if (PBEFS_STRUCTURES STREQUAL "boost")
    target_link_libraries(${modName} PUBLIC Boost::headers)
endif ()
add_library(${modNamespace}::${modName} ALIAS ${modName})

if (OpenMP_CXX_FOUND)
    add_library(${modName}_OpenMP INTERFACE)
    target_link_libraries(${modName}_OpenMP INTERFACE
            ${modName} OpenMP::OpenMP_CXX)
    set_target_properties(${modName}_OpenMP PROPERTIES EXPORT_NAME OpenMP)
    add_library(${modNamespace}::OpenMP ALIAS ${modName}_OpenMP)
    list(APPEND runtimes OpenMP)
endif ()

if (MPI_CXX_FOUND)
    add_library(${modName}_MPI INTERFACE)
    target_link_libraries(${modName}_MPI INTERFACE ${modName} MPI::MPI_CXX)
    target_compile_definitions(${modName}_MPI INTERFACE ADS_WITH_MPI)
    set_target_properties(${modName}_MPI PROPERTIES EXPORT_NAME MPI)
    add_library(${modNamespace}::MPI ALIAS ${modName}_MPI)
    list(APPEND runtimes MPI)
endif ()

# The headers of FastFlow are included as system headers, so that their
# warnings do not drown ours; installed, they go next to ours.
if (PBEFS_HAS_FASTFLOW)
    add_library(${modName}_FastFlow INTERFACE)
    target_include_directories(
            ${modName}_FastFlow SYSTEM INTERFACE
            $<BUILD_INTERFACE:${FASTFLOW_INCLUDE_DIR}>
            $<INSTALL_INTERFACE:include/${modNamespace}/fastflow>)
    target_link_libraries(${modName}_FastFlow INTERFACE ${modName})
    # FastFlow would pin its threads from the first core on, whatever the
    # CPUs given to the process: the workers are placed by pin_worker()
    target_compile_definitions(${modName}_FastFlow INTERFACE
            NO_DEFAULT_MAPPING ADS_WITH_FASTFLOW)
    set_target_properties(${modName}_FastFlow PROPERTIES
            EXPORT_NAME FastFlow)
    add_library(${modNamespace}::FastFlow ALIAS ${modName}_FastFlow)
    list(APPEND runtimes FastFlow)
endif ()

# The symbolic searches are in the headers but for the diagrams of Sylvan,
# which src/Bdd.cpp wraps, so that their users do not include the headers of
# Sylvan; the library is a component of its own, as the runtimes.
if (PBEFS_HAS_SYLVAN)
    add_library(${modName}_Symbolic src/Bdd.cpp)
    target_link_libraries(${modName}_Symbolic PUBLIC ${modName}
            PRIVATE sylvan::sylvan)
    target_compile_definitions(${modName}_Symbolic PUBLIC PBEFS_HAS_SYLVAN)
    set_target_properties(${modName}_Symbolic PROPERTIES
            EXPORT_NAME Symbolic
            OUTPUT_NAME ${modName}Symbolic
            VERSION ${PROJECT_VERSION}
            SOVERSION ${soVersion}
            WINDOWS_EXPORT_ALL_SYMBOLS ON
            POSITION_INDEPENDENT_CODE ON)
    add_library(${modNamespace}::Symbolic ALIAS ${modName}_Symbolic)
    list(APPEND runtimes Symbolic)
endif ()

message(STATUS "${modName} ${PROJECT_VERSION}: OpenMP ${OpenMP_CXX_FOUND}, "
        "MPI ${MPI_CXX_FOUND}, FastFlow ${PBEFS_HAS_FASTFLOW}, "
        "Sylvan ${PBEFS_HAS_SYLVAN}, structures ${PBEFS_STRUCTURES}")

# ----- Subdirectories ------------------------------------------------------ #
if (PBEFS_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif ()

if (PBEFS_BUILD_TOOLS)
    add_subdirectory(tools)
endif ()

if (PBEFS_BUILD_TESTS AND BUILD_TESTING)
    add_subdirectory(test)
endif ()

# ----- Install instructions ------------------------------------------------ #
# The following commands are used when installing the library
# and its CMake configuration files on the system. The target of each
# runtime has an import file of its own, which the package configuration
# file reads only if the consumer finds that runtime too, so that a
# runtime missing there removes only its target.
install(TARGETS ${modName} EXPORT ${modName}Targets)
foreach (runtime IN LISTS runtimes)
    install(TARGETS ${modName}_${runtime}
            EXPORT ${modName}${runtime}Targets)
    install(EXPORT ${modName}${runtime}Targets
            NAMESPACE ${modNamespace}::
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})
endforeach ()

# Install the headers and the auto-generated configuration header.
install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${PROJECT_BINARY_DIR}/${modName}Config.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace})
if (PBEFS_HAS_FASTFLOW)
    install(DIRECTORY ${FASTFLOW_INCLUDE_DIR}/ff
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace}/fastflow)
endif ()

# (Generate and) install the target import file, that allows other
# CMake projects to import the targets.
install(EXPORT ${modName}Targets
        NAMESPACE ${modNamespace}::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Generate the package version file, that allows other
# CMake projects to know the version.
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}ConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion)

# Generate the package configuration file, that allows other
# CMake projects to find the library with find_package().
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/${modName}Config.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}Config.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Install the package version and configuration files.
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}ConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Generate and install the pkg-config file pbefs.pc, which gives the
# library alone (the runtimes are for CMake), with the macro of the data
# structures it was built with. The prefix is found from the directory of
# the file, so that a package installed elsewhere than CMAKE_INSTALL_PREFIX
# (e.g., by cmake --install --prefix) gives its own paths.
if (IS_ABSOLUTE "${CMAKE_INSTALL_LIBDIR}")
    set(pcPrefix "${CMAKE_INSTALL_PREFIX}")
else ()
    set(pcPrefix "\${pcfiledir}/..")
    string(REGEX MATCHALL "[^/]+" pcLevels "${CMAKE_INSTALL_LIBDIR}")
    foreach (pcLevel IN LISTS pcLevels)
        string(APPEND pcPrefix "/..")
    endforeach ()
endif ()
if (IS_ABSOLUTE "${CMAKE_INSTALL_INCLUDEDIR}")
    set(pcIncludeDir "${CMAKE_INSTALL_INCLUDEDIR}")
else ()
    set(pcIncludeDir "\${prefix}/${CMAKE_INSTALL_INCLUDEDIR}")
endif ()
if (IS_ABSOLUTE "${CMAKE_INSTALL_LIBDIR}")
    set(pcLibDir "${CMAKE_INSTALL_LIBDIR}")
else ()
    set(pcLibDir "\${prefix}/${CMAKE_INSTALL_LIBDIR}")
endif ()
list(JOIN PBEFS_ADS_DEFINITIONS " -D" pcCflags)
set(pcCflags " -D${pcCflags}")
set(pcRequires "ads")
if (NOT PBEFS_STRUCTURES STREQUAL "ads")
    string(APPEND pcCflags " -DPBEFS_STRUCTURES_${structures}")
endif ()
configure_file(cmake/pbefs.pc.in ${CMAKE_CURRENT_BINARY_DIR}/pbefs.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pbefs.pc
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

# Install the README and LICENSE files.
install(FILES
        ${CMAKE_CURRENT_LIST_DIR}/README.md
        ${CMAKE_CURRENT_LIST_DIR}/LICENSE
        DESTINATION ${CMAKE_INSTALL_DATADIR}/${modName})

# ----- Packaging support --------------------------------------------------- #
# The following commands allow one to build a proper package using CPack.
# See: https://cmake.org/cmake/help/latest/manual/cpack.1.html
if (NOT hasParent)
    set(CPACK_PACKAGE_VENDOR "Donato Meoli")
    set(CPACK_PACKAGE_CONTACT "Donato Meoli")
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY ${PROJECT_DESCRIPTION})
    set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_LIST_DIR}/LICENSE")
    set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_LIST_DIR}/README.md")
    set(CPACK_SOURCE_IGNORE_FILES
            \\\\.git/
            \\\\.gitignore$
            \\\\.gitattributes$)
    include(CPack)
endif ()

# --------------------------------------------------------------------------- #
# Remove from the search paths the cmake directory we added at the beginning.
list(REMOVE_AT CMAKE_MODULE_PATH -1)

# --------------------------------------------------------------------------- #
