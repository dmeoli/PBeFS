# --------------------------------------------------------------------------- #
#    CMake file for PBeFS/tools                                               #
#                                                                             #
#    AStarBench, the benchmark of the shared-memory searches, and, if MPI is  #
#    there, AStarBenchMPI, the same sources compiled for HDAStarMPI;          #
#    BranchAndBoundBench and BranchAndBoundBenchMPI are the same for the      #
#    Branch-and-Bound; PDBBench is the benchmark of the pattern databases.    #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Università di Pisa                              #
# --------------------------------------------------------------------------- #

# the searches of AStarBench are compiled by AStarBenchPart.cpp once for
# each Domain of AStarBenchDomains.h and each of its three variants, through
# a unit generated here for each pair, so that no unit compiles the searches
# of more than one Domain type
file(STRINGS AStarBenchDomains.h benchDomains REGEX "^PBEFS_BENCH_DOMAIN")
list(LENGTH benchDomains benchDomainCount)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
        AStarBenchDomains.h)
set(PBEFS_ASTARBENCH_SOURCES AStarBench.cpp)
foreach (i RANGE 1 ${benchDomainCount})
    foreach (part 0 1 2)
        set(unit ${CMAKE_CURRENT_BINARY_DIR}/AStarBenchPart-${i}-${part}.cpp)
        file(CONFIGURE OUTPUT ${unit} @ONLY CONTENT
                "#define PBEFS_BENCH_DOMAIN_INDEX ${i}
#define PBEFS_BENCH_PART ${part}
#include \"${CMAKE_CURRENT_SOURCE_DIR}/AStarBenchPart.cpp\"
")
        list(APPEND PBEFS_ASTARBENCH_SOURCES ${unit})
    endforeach ()
endforeach ()
add_executable(AStarBench ${PBEFS_ASTARBENCH_SOURCES})
target_link_libraries(AStarBench PRIVATE ${modNamespace}::${modName})
add_executable(BranchAndBoundBench BranchAndBoundBench.cpp)
target_link_libraries(BranchAndBoundBench PRIVATE ${modNamespace}::${modName})
add_executable(PDBBench PDBBench.cpp)
target_link_libraries(PDBBench PRIVATE ${modNamespace}::${modName})
foreach (bench AStarBench BranchAndBoundBench PDBBench)
    target_compile_options(${bench} PRIVATE ${PBEFS_CXX_OPTIONS})
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
    add_executable(AStarBenchMPI ${PBEFS_ASTARBENCH_SOURCES})
    target_link_libraries(AStarBenchMPI PRIVATE ${modNamespace}::MPI)
    target_compile_definitions(AStarBenchMPI PRIVATE PBEFS_BENCH_MPI)
    add_executable(BranchAndBoundBenchMPI BranchAndBoundBench.cpp)
    target_link_libraries(BranchAndBoundBenchMPI PRIVATE
            ${modNamespace}::MPI)
    target_compile_definitions(BranchAndBoundBenchMPI PRIVATE
            PBEFS_BENCH_MPI)
    foreach (bench AStarBenchMPI BranchAndBoundBenchMPI)
        target_compile_options(${bench} PRIVATE ${PBEFS_CXX_OPTIONS})
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
