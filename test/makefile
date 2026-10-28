###############################################################################
################################## makefile ###################################
###############################################################################
#                                                                             #
#   makefile of test (for PBeFS)                                              #
#                                                                             #
#   make builds the tester with the threads, OpenMP and FastFlow, make mpi    #
#   the one with MPI too (compiled with mpicxx), make check runs all the      #
#   cases as ctest does.                                                      #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
#                                                                             #
###############################################################################

# module name
NAME = PBeFS_test

# basic directory
DIR = .

# common flags
COMMON_SW = -std=c++20 -Wall -Wextra -Wpedantic

# debug switches
SW_DEBUG = -g3 $(PBeFSSWDBG) -O0 $(COMMON_SW) $(PBeFSSWERR)
# debug switches with address sanitizer
SW_DEBUG_ASAN = -g3 -O1 $(COMMON_SW) -fsanitize=undefined -fsanitize=address \
                -fno-omit-frame-pointer
# production switches
SW_RELEASE = -O3 $(COMMON_SW) -DNDEBUG

# compiler
CC = g++
MPICC = mpicxx

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: release

# debug target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

debug: SW = $(SW_DEBUG)
debug: $(DIR)/$(NAME)

# debug_asan target - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

debug_asan: SW = $(SW_DEBUG_ASAN)
debug_asan: $(DIR)/$(NAME)

# release target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

release: SW = $(SW_RELEASE)
release: $(DIR)/$(NAME)

# mpi target: the tester with MPI too - - - - - - - - - - - - - - - - - - - -

mpi: SW = $(SW_RELEASE)
mpi: $(DIR)/$(NAME)_mpi

# check target: all the cases - - - - - - - - - - - - - - - - - - - - - - - -

# the cases of tests.txt, those with MPI apart, and the check of the list

CASES = $(shell awk '/^[a-z]/ && $$2 != "mpi" { print $$1 }' tests.txt)
MPICASES = $(shell awk '/^[a-z]/ && $$2 == "mpi" { print $$1 }' tests.txt)

check: release
	@for c in $(CASES); do $(DIR)/$(NAME) $$c ../data | tail -1 || exit 1; done
	@$(DIR)/$(NAME) registry tests.txt | tail -1

check_mpi: mpi
	@for c in $(MPICASES); do for n in 1 2 4; do \
	mpiexec -n $$n $(DIR)/$(NAME)_mpi $$c ../data | tail -1 || exit 1; \
	done; done

# clean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)/*.o $(DIR)/*~ $(DIR)/$(NAME) $(DIR)/$(NAME)_mpi

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# each module outputs some macros to be used here:
# *OBJ is the final object(s) / library
# *LIB is the external libraries + -L< libdirs >
# *H   is the list of all include files
# *INC is the -I< include directories >

# PBeFS
PBeFSSDR = ..
include $(PBeFSSDR)/makefile-c

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# object files
MOBJ = $(PBeFSOBJ)

# libraries
MLIB = $(PBeFSLIB)

$(DIR)/$(NAME): $(DIR)/test.cpp $(PBeFSH) $(MOBJ)
	$(CC) -o $@ $(DIR)/test.cpp $(MOBJ) $(PBeFSINC) -I$(DIR) $(PBeFSDEF) \
	$(MLIB) $(SW)

$(DIR)/$(NAME)_mpi: $(DIR)/test.cpp $(PBeFSH) $(PBeFSMPIH) $(MOBJ)
	$(MPICC) -o $@ $(DIR)/test.cpp $(MOBJ) $(PBeFSINC) -I$(DIR) $(PBeFSDEF) \
	$(PBeFSMPIDEF) $(MLIB) $(SW)

# the configuration header that CMake writes in the build directory
$(DIR)/PBeFSConfig.h:
	v=$$(git -C $(PBeFSSDR) describe --tags --abbrev=0 2>/dev/null) ; \
	printf '#define PBeFS_VERSION "%s"\n' "$${v:-0.0.0}" > $@

$(DIR)/$(NAME) $(DIR)/$(NAME)_mpi: $(DIR)/PBeFSConfig.h

# distclean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

distclean: clean
	rm -f $(DIR)/PBeFSConfig.h

# phony targets - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.PHONY: debug debug_asan release mpi check check_mpi clean distclean

############################### End of makefile ###############################
