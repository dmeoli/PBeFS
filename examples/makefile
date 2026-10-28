###############################################################################
################################## makefile ###################################
###############################################################################
#                                                                             #
#   makefile of the examples of PBeFS: 8Puzzle, 15Puzzle, FindPath,           #
#   RouteFinding, NQueens, WaterJugs, MissionariesCannibals, VacuumWorld      #
#   and TravelingSalesman                                                     #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
#                                                                             #
###############################################################################

# basic directory
DIR = .

# the programs
NAMES = 8Puzzle 15Puzzle FindPath RouteFinding NQueens WaterJugs \
        MissionariesCannibals VacuumWorld \
        TravelingSalesman

# common flags
COMMON_SW = -std=c++20 -Wall -Wextra

# debug switches
SW_DEBUG = -g3 $(PBeFSSWDBG) -O0 $(COMMON_SW) $(PBeFSSWERR)
# production switches
SW_RELEASE = -O3 $(COMMON_SW) -DNDEBUG

# compiler
CC = g++

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: release

# debug target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

debug: SW = $(SW_DEBUG)
debug: $(addprefix $(DIR)/,$(NAMES))

# release target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

release: SW = $(SW_RELEASE)
release: $(addprefix $(DIR)/,$(NAMES))

# clean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)/*.o $(DIR)/*~ $(addprefix $(DIR)/,$(NAMES))

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -

# PBeFS (without makefile-c: the examples use no parallel runtime
# save the C++ threads)
PBeFSSDR = ..
include $(PBeFSSDR)/makefile

# the programs- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

$(DIR)/%: $(DIR)/%.cpp $(PBeFSH) $(PBeFSOBJ)
	$(CC) -o $@ $< $(PBeFSOBJ) $(PBeFSINC) -pthread $(SW)

# distclean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

distclean: clean

# phony targets - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.PHONY: debug release clean distclean

############################### End of makefile ###############################
