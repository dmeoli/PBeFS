# --------------------------------------------------------------------------- #
# DeriveVersion.cmake                                                         #
#                                                                             #
# Derive the module version from git, so that cutting a release is just       #
# `git tag x.y.z` with NO version number written by hand anywhere: the        #
# version in project(), in the generated <Module>Config.h and in the tools'   #
# --version all follow the tag.                                               #
#                                                                             #
# Resolution order (the first that yields a numeric x[.y[.z]] wins):          #
#   1. a source tarball (e.g. the one that release.yml attaches to a GitHub   #
#      release, NOT a git repo): the VERSION.txt file, which git fills in at  #
#      `git archive` time via the `VERSION.txt export-subst` .gitattributes   #
#      entry (its content is "$Format:%(describe:tags,abbrev=0)$",            #
#      substituted with the tag by git);                                      #
#   2. a working git clone whose top level is the source directory:           #
#      `git describe --tags --abbrev=0`; a tarball unpacked inside another    #
#      git repository (the build trees of vcpkg, a FetchContent of the        #
#      tarball under the build directory of a project) is not asked, since    #
#      the tags would be those of that repository;                            #
#   3. otherwise the sentinel "0.0.0" (no version info available).            #
#                                                                             #
# NOTE: the file is VERSION.txt, NOT VERSION: libc++'s C++20 <version> header #
# would otherwise be shadowed by it on case-insensitive filesystems (macOS,   #
# Windows) whenever the repo root is on a -I include path.                    #
#                                                                             #
# The leading x[.y[.z]] is extracted from whatever git returns, so a describe #
# with a distance suffix ("1.2.3-4-gdeadbee", e.g. an archive of a non-tag    #
# commit) still yields "1.2.3" rather than the sentinel.                      #
#                                                                             #
# Call it BEFORE project():                                                   #
#     include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/DeriveVersion.cmake)          #
#     pbefs_derive_version(MOD_VERSION)                                       #
#     project(<Module> VERSION ${MOD_VERSION} ...)                            #
# --------------------------------------------------------------------------- #

function(pbefs_derive_version out_var)
    set(_raw "")

    # 1) tarball: the VERSION.txt file filled in by `git archive`
    #    (export-subst); in a git clone it still holds the literal
    #    "$Format:...$", which is left to git
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION.txt")
        file(READ "${CMAKE_CURRENT_SOURCE_DIR}/VERSION.txt" _raw)
        string(STRIP "${_raw}" _raw)
        if(_raw MATCHES "^\\$Format:")
            set(_raw "")
        endif()
    endif()

    # 2) a git clone of this repository: the prefix of the source directory
    #    in its work tree is empty only at the top level
    find_package(Git QUIET)
    if(NOT _raw AND GIT_FOUND)
        execute_process(
                COMMAND "${GIT_EXECUTABLE}" rev-parse --show-prefix
                WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                OUTPUT_VARIABLE _prefix
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
                RESULT_VARIABLE _result)
        if(_result EQUAL 0 AND _prefix STREQUAL "")
            execute_process(
                    COMMAND "${GIT_EXECUTABLE}" describe --tags --abbrev=0
                    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                    OUTPUT_VARIABLE _raw
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET
                    RESULT_VARIABLE _result)
            if(NOT _result EQUAL 0)
                set(_raw "")
            endif()
        endif()
    endif()

    # extract the leading x[.y[.z]] (tolerating a leading 'v' and a trailing
    # "-<n>-g<sha>" describe suffix); 3) sentinel if nothing matched
    if(_raw MATCHES "^v?([0-9]+(\\.[0-9]+)*)")
        set(_version "${CMAKE_MATCH_1}")
    else()
        set(_version "0.0.0")
    endif()

    set(${out_var} "${_version}" PARENT_SCOPE)
endfunction()
