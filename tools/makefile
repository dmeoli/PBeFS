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
SW_DEBUG = -g3 -O0 $(COMMON_SW)
# production switches
SW_RELEASE = -O3 $(COMMON_SW) -DNDEBUG

# compilers
CC = g++
MPICC = mpicxx

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

# AStarBench, with its searches compiled by AStarBenchPart.cpp once for each
# Domain of AStarBenchDomains.h and each of its three variants, the object
# AStarBenchPart-<Domain>-<variant>.o, so that no unit compiles the searches
# of more than one Domain type
ASTARBENCHH = $(DIR)/AStarBench.h $(DIR)/AStarBenchSolve.h \
              $(DIR)/AStarBenchDomains.h
ASTARBENCHN := $(shell seq $$(grep -c ^PBEFS_BENCH_DOMAIN \
               $(DIR)/AStarBenchDomains.h))
ASTARBENCHPARTS = $(foreach i,$(ASTARBENCHN),$(foreach p,0 1 2,$(i)-$(p)))
ASTARBENCHPARTOBJ = $(ASTARBENCHPARTS:%=$(DIR)/AStarBenchPart-%.o)
ASTARBENCHMPIPARTOBJ = $(ASTARBENCHPARTS:%=$(DIR)/AStarBenchMPIPart-%.o)
ASTARBENCHOBJ = $(DIR)/AStarBench.o $(ASTARBENCHPARTOBJ)
ASTARBENCHMPIOBJ = $(DIR)/AStarBench-mpi.o $(ASTARBENCHMPIPARTOBJ)
# the Domain and the variant of an object, from the stem <Domain>-<variant>
ASTARBENCHPARTDEF = -DPBEFS_BENCH_DOMAIN_INDEX=$(word 1,$(subst -, ,$*)) \
                    -DPBEFS_BENCH_PART=$(word 2,$(subst -, ,$*))

$(DIR)/AStarBench.o: $(DIR)/AStarBench.cpp $(ASTARBENCHH) $(PBeFSH)
	$(CC) -c -o $@ $< $(PBeFSINC) $(PBeFSDEF) $(SW)

$(ASTARBENCHPARTOBJ): $(DIR)/AStarBenchPart-%.o: $(DIR)/AStarBenchPart.cpp \
                      $(ASTARBENCHH) $(PBeFSH)
	$(CC) -c -o $@ $< $(ASTARBENCHPARTDEF) $(PBeFSINC) $(PBeFSDEF) $(SW)

$(DIR)/AStarBench-mpi.o: $(DIR)/AStarBench.cpp $(ASTARBENCHH) $(PBeFSH) \
                         $(PBeFSMPIH)
	$(MPICC) -c -o $@ $< $(PBeFSINC) $(PBeFSMPIDEF) -DPBEFS_BENCH_MPI \
	-pthread $(SW)

$(ASTARBENCHMPIPARTOBJ): $(DIR)/AStarBenchMPIPart-%.o: \
                         $(DIR)/AStarBenchPart.cpp $(ASTARBENCHH) \
                         $(PBeFSH) $(PBeFSMPIH)
	$(MPICC) -c -o $@ $< $(ASTARBENCHPARTDEF) $(PBeFSINC) $(PBeFSMPIDEF) \
	-DPBEFS_BENCH_MPI -pthread $(SW)

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
