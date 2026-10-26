# --------------------------------------------------------------------------- #
#    CMake file for PBeFS/tools                                               #
#                                                                             #
#    AStarBench, the benchmark of the shared-memory searches, and, if MPI is  #
#    there, AStarBenchMPI, the same source compiled for HDAStarMPI;           #
#    BranchAndBoundBench and BranchAndBoundBenchMPI are the same for the      #
#    Branch-and-Bound; PDBBench is the benchmark of the pattern databases.    #
#    All are compiled with the alignment of the branches against the erratum  #
#    JCC of Skylake, where the assembler takes it.                            #
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
add_executable(BranchAndBoundBench BranchAndBoundBench.cpp)
target_link_libraries(BranchAndBoundBench PRIVATE ${modNamespace}::${modName})
add_executable(PDBBench PDBBench.cpp)
target_link_libraries(PDBBench PRIVATE ${modNamespace}::${modName})
foreach (bench AStarBench BranchAndBoundBench PDBBench)
    # the instantiations of every search on every domain exceed the sections
    # of an object file of MSVC
    target_compile_options(${bench} PRIVATE
            $<$<CXX_COMPILER_ID:MSVC>:/bigobj>)
    if (TARGET ${modName}_OpenMP)
        target_link_libraries(${bench} PRIVATE ${modNamespace}::OpenMP)
        target_compile_definitions(${bench} PRIVATE PBEFS_HAS_OPENMP)
    endif ()
    if (TARGET ${modName}_FastFlow)
        target_link_libraries(${bench} PRIVATE ${modNamespace}::FastFlow)
        target_compile_definitions(${bench} PRIVATE PBEFS_HAS_FASTFLOW)
    endif ()
endforeach ()

if (TARGET ${modName}_MPI)
    add_executable(AStarBenchMPI AStarBench.cpp)
    target_link_libraries(AStarBenchMPI PRIVATE ${modNamespace}::MPI)
    target_compile_definitions(AStarBenchMPI PRIVATE PBEFS_BENCH_MPI)
    add_executable(BranchAndBoundBenchMPI BranchAndBoundBench.cpp)
    target_link_libraries(BranchAndBoundBenchMPI PRIVATE
            ${modNamespace}::MPI)
    target_compile_definitions(BranchAndBoundBenchMPI PRIVATE
            PBEFS_BENCH_MPI)
endif ()

if (PBEFS_HAS_JCC)
    foreach (bench AStarBench AStarBenchMPI BranchAndBoundBench
            BranchAndBoundBenchMPI PDBBench)
        if (TARGET ${bench})
            target_compile_options(${bench} PRIVATE
                    -Wa,-mbranches-within-32B-boundaries)
        endif ()
    endforeach ()
endif ()

install(TARGETS AStarBench BranchAndBoundBench PDBBench
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
foreach (bench AStarBenchMPI BranchAndBoundBenchMPI)
    if (TARGET ${bench})
        install(TARGETS ${bench} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    endif ()
endforeach ()

# --------------------------------------------------------------------------- #
