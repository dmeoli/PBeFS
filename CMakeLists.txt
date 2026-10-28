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
#    The parallel runtimes are optional: OpenMP and MPI are used if found,    #
#    FastFlow if its headers are in FASTFLOW_ROOT (/opt/fastflow by default)  #
#    or on the include path, else it is downloaded at a fixed commit; each    #
#    can be turned off with the PBEFS_WITH_* options. PBEFS_STRUCTURES        #
#    chooses the data structures of A*: those of the library (own, the        #
#    default), of the standard library (std), of Boost (boost) or of ADS      #
#    (ads), the latter found by find_package(ADS) or else downloaded at a     #
#    fixed commit.                                                            #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
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

option(PBEFS_FETCH_FASTFLOW "Download FastFlow if it is not installed" ON)

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
# commit of ADS downloaded when they are those of ADS and it is not installed
set(PBEFS_STRUCTURES own CACHE STRING
    "The data structures of A*: own, std, boost or ads")
set_property(CACHE PBEFS_STRUCTURES PROPERTY STRINGS own std boost ads)
set(ADS_ROOT /opt/ADS CACHE PATH "Prefix where ADS is installed")
pbefs_pin(PBEFS_ADS_COMMIT 45312c537e947f4d9c8ba35bcc4c3deccfb39c4f
          "The commit of ADS to download")

# ----- Requirements -------------------------------------------------------- #
find_package(Threads REQUIRED)

# Boost is header-only here (Boost::headers); ADS is installed (ADS_ROOT is
# a hint), or else the sources of the commit PBEFS_ADS_COMMIT are added as
# a subproject.
if (PBEFS_STRUCTURES STREQUAL "boost")
    find_package(Boost 1.81 REQUIRED)
elseif (PBEFS_STRUCTURES STREQUAL "ads")
    find_package(ADS CONFIG QUIET HINTS ${ADS_ROOT})
    if (NOT ADS_FOUND)
        include(FetchContent)
        FetchContent_Declare(
                ADS
                GIT_REPOSITORY https://github.com/dmeoli/ADS.git
                GIT_TAG ${PBEFS_ADS_COMMIT})
        FetchContent_MakeAvailable(ADS)
    endif ()
    if (NOT TARGET ADS::ADS)
        message(FATAL_ERROR "PBEFS_STRUCTURES ads needs ADS: install it or "
                "set ADS_ROOT")
    endif ()
elseif (NOT PBEFS_STRUCTURES MATCHES "^(own|std)$")
    message(FATAL_ERROR "PBEFS_STRUCTURES is own, std, boost or ads, not "
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
        src/Distribution.cpp
        src/Grid.cpp
        src/SearchUtils.cpp
        src/Solution.cpp)

target_include_directories(
        ${modName} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:include>)
target_link_libraries(${modName} PUBLIC Threads::Threads)
if (NOT PBEFS_STRUCTURES STREQUAL "own")
    string(TOUPPER ${PBEFS_STRUCTURES} structures)
    target_compile_definitions(${modName} PUBLIC
            PBEFS_STRUCTURES_${structures})
endif ()
if (PBEFS_STRUCTURES STREQUAL "boost")
    target_link_libraries(${modName} PUBLIC Boost::headers)
elseif (PBEFS_STRUCTURES STREQUAL "ads")
    target_link_libraries(${modName} PUBLIC ADS::ADS)
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
            NO_DEFAULT_MAPPING)
    set_target_properties(${modName}_FastFlow PROPERTIES
            EXPORT_NAME FastFlow)
    add_library(${modNamespace}::FastFlow ALIAS ${modName}_FastFlow)
    list(APPEND runtimes FastFlow)
endif ()

message(STATUS "${modName} ${PROJECT_VERSION}: OpenMP ${OpenMP_CXX_FOUND}, "
        "MPI ${MPI_CXX_FOUND}, FastFlow ${PBEFS_HAS_FASTFLOW}, "
        "structures ${PBEFS_STRUCTURES}")

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
set(pcCflags "")
set(pcRequires "")
if (NOT PBEFS_STRUCTURES STREQUAL "own")
    set(pcCflags " -DPBEFS_STRUCTURES_${structures}")
endif ()
if (PBEFS_STRUCTURES STREQUAL "ads")
    set(pcRequires "ads")
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
