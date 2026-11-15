# --------------------------------------------------------------------------- #
#    CMake file for PBeFS/tools                                               #
#                                                                             #
#    AStarBench, the benchmark of the shared-memory searches, and, if MPI is  #
#    there, AStarBenchMPI, the same sources compiled for HDAStarMPI;          #
#    BranchAndBoundBench and BranchAndBoundBenchMPI are the same for the      #
#    Branch-and-Bound; PDBBench is the benchmark of the pattern databases,    #
#    and SymbolicBench, with Sylvan, that of the symbolic searches. All are   #
#    compiled with the alignment of the branches against the erratum JCC of   #
#    Skylake, where the assembler takes it.                                   #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Università di Pisa                              #
# --------------------------------------------------------------------------- #

# On x86 the branches are kept off the 32-byte boundaries, where the
# erratum JCC of Skylake slows a loop by the alignment alone (about 10%),
# if the assembler takes the option; the searches and the data structures
# they are compared on (std, Boost, ADS) are all in the headers compiled
# here, and hence all of them get it alike.
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-Wa,-mbranches-within-32B-boundaries PBEFS_HAS_JCC)

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

if (TARGET ${modName}_Symbolic)
    add_executable(SymbolicBench SymbolicBench.cpp)
    target_link_libraries(SymbolicBench PRIVATE ${modNamespace}::Symbolic)
    foreach (runtime OpenMP FastFlow)
        if (TARGET ${modName}_${runtime})
            string(TOUPPER ${runtime} macro)
            target_link_libraries(SymbolicBench PRIVATE
                    ${modNamespace}::${runtime})
            target_compile_definitions(SymbolicBench PRIVATE
                    PBEFS_HAS_${macro})
        endif ()
    endforeach ()
    install(TARGETS SymbolicBench RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif ()

if (TARGET ${modName}_MPI)
    add_executable(AStarBenchMPI ${PBEFS_ASTARBENCH_SOURCES})
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
            BranchAndBoundBenchMPI PDBBench SymbolicBench)
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
