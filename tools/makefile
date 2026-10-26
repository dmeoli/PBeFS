##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of the tools of PBeFS: AStarBench, with the threads,            #
#   OpenMP and FastFlow, and AStarBenchMPI, the same source compiled with    #
#   mpicxx for HDAStarMPI                                                    #
#                                                                            #
#                                Donato Meoli                                #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# basic directory
DIR = .

# common flags
COMMON_SW = -std=c++20 -Wall -Wextra

# debug switches
SW_DEBUG = -g -O0 $(COMMON_SW)
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
debug: $(DIR)/AStarBench

# release target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

release: SW = $(SW_RELEASE)
release: $(DIR)/AStarBench

# mpi target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

mpi: SW = $(SW_RELEASE)
mpi: $(DIR)/AStarBenchMPI

# clean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)/*.o $(DIR)/*~ $(DIR)/AStarBench $(DIR)/AStarBenchMPI

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -

# PBeFS
PBeFSSDR = ..
include $(PBeFSSDR)/makefile-c

# the programs- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

$(DIR)/AStarBench: $(DIR)/AStarBench.cpp $(PBeFSH) $(PBeFSOBJ)
	$(CC) -o $@ $< $(PBeFSOBJ) $(PBeFSINC) $(PBeFSDEF) $(PBeFSLIB) $(SW)

$(DIR)/AStarBenchMPI: $(DIR)/AStarBench.cpp $(PBeFSH) $(PBeFSMPIH) $(PBeFSOBJ)
	$(MPICC) -o $@ $< $(PBeFSOBJ) $(PBeFSINC) -DPBEFS_BENCH_MPI -pthread $(SW)

# distclean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

distclean: clean

# phony targets - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.PHONY: debug release mpi clean distclean

############################ End of makefile #################################
