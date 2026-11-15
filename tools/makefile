###############################################################################
################################## makefile ###################################
###############################################################################
#                                                                             #
#   makefile of the tools of PBeFS: AStarBench and BranchAndBoundBench,       #
#   with the threads, OpenMP and FastFlow, AStarBenchMPI and                  #
#   BranchAndBoundBenchMPI, the same sources compiled with mpicxx for         #
#   HDAStarMPI and MPIBranchAndBoundSearch, and PDBBench, the benchmark of    #
#   the pattern databases                                                     #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Università di Pisa                              #
#                                                                             #
###############################################################################

# basic directory
DIR = .

# common flags
COMMON_SW = -std=c++20 -Wall -Wextra

# debug switches
SW_DEBUG = -g3 $(PBeFSSWDBG) -O0 $(COMMON_SW) $(PBeFSSWERR)
# production switches, with JCCFLAGS [see below]
SW_RELEASE = -O3 $(COMMON_SW) -DNDEBUG $(JCCFLAGS)

# compilers
CC = g++
MPICC = mpicxx

# on x86 the branches are kept off the 32-byte boundaries, where the
# erratum JCC of Skylake slows a loop by the alignment alone (about 10%),
# if the assembler takes the option, for the searches and for the data
# structures they are compared on alike
JCCFLAGS := $(shell echo 'int main(){}' | $(CC) -x c++ -c -o /dev/null \
            -Wa,-mbranches-within-32B-boundaries - 2>/dev/null && \
            echo -Wa,-mbranches-within-32B-boundaries)

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: release

# debug target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

debug: SW = $(SW_DEBUG)
debug: $(DIR)/AStarBench $(DIR)/BranchAndBoundBench $(DIR)/PDBBench

# release target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

release: SW = $(SW_RELEASE)
release: $(DIR)/AStarBench $(DIR)/BranchAndBoundBench $(DIR)/PDBBench

# mpi target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

mpi: SW = $(SW_RELEASE)
mpi: $(DIR)/AStarBenchMPI $(DIR)/BranchAndBoundBenchMPI

# clean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)/*.o $(DIR)/*~ $(DIR)/AStarBench $(DIR)/AStarBenchMPI \
	      $(DIR)/BranchAndBoundBench $(DIR)/BranchAndBoundBenchMPI \
	      $(DIR)/PDBBench

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -

# PBeFS
PBeFSSDR = ..
include $(PBeFSSDR)/makefile-c

# the programs- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# AStarBench and the files that instantiate its searches by Domain
ASTARBENCHSRC = $(DIR)/AStarBench.cpp \
                $(DIR)/AStarBenchTiles.cpp \
                $(DIR)/AStarBenchGrids.cpp \
                $(DIR)/AStarBenchGraphs.cpp \
                $(DIR)/AStarBenchPancakes.cpp \
                $(DIR)/AStarBenchBigPancakes.cpp \
                $(DIR)/AStarBenchBurned.cpp \
                $(DIR)/AStarBenchBigBurned.cpp \
                $(DIR)/AStarBenchTopSpin.cpp \
                $(DIR)/AStarBenchBigTopSpin.cpp \
                $(DIR)/AStarBenchHanoi.cpp \
                $(DIR)/AStarBenchExamples.cpp
ASTARBENCHH = $(DIR)/AStarBench.h $(DIR)/AStarBenchSolve.h
ASTARBENCHOBJ = $(ASTARBENCHSRC:.cpp=.o)
ASTARBENCHMPIOBJ = $(ASTARBENCHSRC:.cpp=-mpi.o)

$(ASTARBENCHOBJ): %.o: %.cpp $(ASTARBENCHH) $(PBeFSH)
	$(CC) -c -o $@ $< $(PBeFSINC) $(PBeFSDEF) $(SW)

$(ASTARBENCHMPIOBJ): %-mpi.o: %.cpp $(ASTARBENCHH) $(PBeFSH) $(PBeFSMPIH)
	$(MPICC) -c -o $@ $< $(PBeFSINC) $(PBeFSMPIDEF) -DPBEFS_BENCH_MPI \
	-pthread $(SW)

$(DIR)/AStarBench: $(ASTARBENCHOBJ) $(PBeFSOBJ)
	$(CC) -o $@ $(ASTARBENCHOBJ) $(PBeFSOBJ) $(PBeFSLIB) $(SW)

$(DIR)/AStarBenchMPI: $(ASTARBENCHMPIOBJ) $(PBeFSOBJ)
	$(MPICC) -o $@ $(ASTARBENCHMPIOBJ) $(PBeFSOBJ) -pthread $(SW)

$(DIR)/BranchAndBoundBench: $(DIR)/BranchAndBoundBench.cpp $(PBeFSH) \
                            $(PBeFSOBJ)
	$(CC) -o $@ $< $(PBeFSOBJ) $(PBeFSINC) $(PBeFSDEF) $(PBeFSLIB) $(SW)

$(DIR)/BranchAndBoundBenchMPI: $(DIR)/BranchAndBoundBench.cpp $(PBeFSH) \
                               $(PBeFSMPIH) $(PBeFSOBJ)
	$(MPICC) -o $@ $< $(PBeFSOBJ) $(PBeFSINC) $(PBeFSMPIDEF) -DPBEFS_BENCH_MPI \
	-pthread $(SW)

$(DIR)/PDBBench: $(DIR)/PDBBench.cpp $(PBeFSH) $(PBeFSOBJ)
	$(CC) -o $@ $< $(PBeFSOBJ) $(PBeFSINC) $(PBeFSDEF) $(PBeFSLIB) $(SW)

# distclean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

distclean: clean

# phony targets - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.PHONY: debug release mpi clean distclean

############################### End of makefile ###############################
