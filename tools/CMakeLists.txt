# --------------------------------------------------------------------------- #
#    CMake file for PBeFS/tools                                               #
#                                                                             #
#    AStarBench, the benchmark of the shared-memory searches, and, if MPI is  #
#    there, AStarBenchMPI, the same source compiled for HDAStarMPI; both are  #
#    compiled with the alignment of the branches against the erratum JCC of   #
#    Skylake, where the assembler takes it.                                   #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

# On x86 the branches are kept off the 32-byte boundaries, where the
# erratum JCC of Skylake slows a loop by the alignment alone (about 10%),
# if the assembler takes the option; the searches and the data structures
# they are compared on (std, Boost, ADS) are all in the headers compiled
# here, and hence all of them get it alike.
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wa,-mbranches-within-32B-boundaries PBEFS_HAS_JCC)

add_executable(AStarBench AStarBench.cpp)
target_link_libraries(AStarBench PRIVATE ${modNamespace}::${modName})
if (TARGET ${modName}_OpenMP)
    target_link_libraries(AStarBench PRIVATE ${modNamespace}::OpenMP)
    target_compile_definitions(AStarBench PRIVATE PBEFS_HAS_OPENMP)
endif ()
if (TARGET ${modName}_FastFlow)
    target_link_libraries(AStarBench PRIVATE ${modNamespace}::FastFlow)
    target_compile_definitions(AStarBench PRIVATE PBEFS_HAS_FASTFLOW)
endif ()

if (TARGET ${modName}_MPI)
    add_executable(AStarBenchMPI AStarBench.cpp)
    target_link_libraries(AStarBenchMPI PRIVATE ${modNamespace}::MPI)
    target_compile_definitions(AStarBenchMPI PRIVATE PBEFS_BENCH_MPI)
endif ()

if (PBEFS_HAS_JCC)
    foreach (bench AStarBench AStarBenchMPI)
        if (TARGET ${bench})
            target_compile_options(${bench} PRIVATE
                    -Wa,-mbranches-within-32B-boundaries)
        endif ()
    endforeach ()
endif ()

install(TARGETS AStarBench RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
if (TARGET AStarBenchMPI)
    install(TARGETS AStarBenchMPI RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif ()

# --------------------------------------------------------------------------- #
