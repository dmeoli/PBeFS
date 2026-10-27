# --------------------------------------------------------------------------- #
#    CMake file for PBeFS/test                                                #
#                                                                             #
#    In this file we define the tests of the library: one tester, test.cpp,   #
#    run once for each case, every case of a parallel runtime only if the     #
#    runtime is there, and the MPI case under mpiexec with 1, 2 and 4         #
#    processes.                                                               #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

# ----- the tester ---------------------------------------------------------- #
add_executable(${modName}_test test.cpp)
target_link_libraries(${modName}_test PRIVATE ${modNamespace}::${modName})
if (MSVC)
    # the instantiations of every search on every domain exceed the sections
    # of an object file of MSVC
    target_compile_options(${modName}_test PRIVATE /W4 /bigobj)
else ()
    target_compile_options(${modName}_test PRIVATE -Wall -Wextra -Wpedantic)
endif ()

if (TARGET ${modName}_OpenMP)
    target_link_libraries(${modName}_test PRIVATE ${modNamespace}::OpenMP)
    target_compile_definitions(${modName}_test PRIVATE PBEFS_HAS_OPENMP)
endif ()
if (TARGET ${modName}_FastFlow)
    target_link_libraries(${modName}_test PRIVATE ${modNamespace}::FastFlow)
    target_compile_definitions(${modName}_test PRIVATE PBEFS_HAS_FASTFLOW)
endif ()
if (TARGET ${modName}_MPI)
    target_link_libraries(${modName}_test PRIVATE ${modNamespace}::MPI)
    target_compile_definitions(${modName}_test PRIVATE PBEFS_HAS_MPI)
endif ()

set(dataDir ${PROJECT_SOURCE_DIR}/data)

# ----- the cases ----------------------------------------------------------- #
# The cases are those of tests.txt whose runtime is built, and the case
# registry checks that tests.txt and the cases of the tester are the same.
# The MPI case runs on 1, 2 and 4 processes; the environment lets Open MPI
# (4 and 5) start more processes than the cores of a small runner, and the
# other implementations ignore it.
set(mpiOversubscribe
        "OMPI_MCA_rmaps_base_oversubscribe=1"
        "PRTE_MCA_rmaps_default_mapping_policy=:oversubscribe")
file(STRINGS tests.txt lines REGEX "^[a-z]")
foreach (line IN LISTS lines)
    string(REGEX MATCH "^([^ ]+) +([^ ]+)" _ "${line}")
    set(case ${CMAKE_MATCH_1})
    set(needs ${CMAKE_MATCH_2})
    if ((needs STREQUAL "openmp" AND NOT TARGET ${modName}_OpenMP) OR
            (needs STREQUAL "fastflow" AND NOT TARGET ${modName}_FastFlow) OR
            (needs STREQUAL "mpi" AND NOT TARGET ${modName}_MPI))
        continue ()
    endif ()
    if (needs STREQUAL "mpi")
        foreach (np 1 2 4)
            add_test(NAME ${modName}.${case}-${np}
                     COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG}
                             ${np} ${MPIEXEC_PREFLAGS}
                             $<TARGET_FILE:${modName}_test>
                             ${MPIEXEC_POSTFLAGS} ${case} ${dataDir})
            set_tests_properties(${modName}.${case}-${np} PROPERTIES
                    ENVIRONMENT "${mpiOversubscribe}")
        endforeach ()
    else ()
        add_test(NAME ${modName}.${case}
                 COMMAND ${modName}_test ${case} ${dataDir})
    endif ()
endforeach ()
add_test(NAME ${modName}.registry
         COMMAND ${modName}_test registry
                 ${CMAKE_CURRENT_SOURCE_DIR}/tests.txt)

# ----- the label of this library's tests ----------------------------------- #
# The label lets the CI select exactly these tests with `ctest -L PBeFS`.
get_property(_tests DIRECTORY PROPERTY TESTS)
if (_tests)
    set_tests_properties(${_tests} PROPERTIES LABELS "${modName}"
                         TIMEOUT 900)
endif ()

# --------------------------------------------------------------------------- #
