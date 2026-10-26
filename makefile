##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of PBeFS                                                        #
#                                                                            #
#   The searches of PBeFS are templates over the Domain, and they are        #
#   in the headers; the objects hold the parts that do not depend on it.     #
#                                                                            #
#   Input:  $(CC)      = compiler command                                    #
#           $(SW)      = compiler options                                    #
#           $(PBeFSSDR)  = the directory where PBeFS is                      #
#           $(PBeFSSTRUCT) = the data structures of A*: std, boost or ads,   #
#                          none for those of PBeFS                           #
#           $(ADSSDR)    = the directory where ADS is, with ads              #
#                                                                            #
#   Output: $(PBeFSOBJ)  = the final object(s) / library                     #
#           $(PBeFSH)    = the .h files to include                           #
#           $(PBeFSINC)  = the -I$( include directory )                      #
#                                                                            #
#                                Donato Meoli                                #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

PBeFSOBJ = $(PBeFSSDR)/obj/Distribution.o \
           $(PBeFSSDR)/obj/Grid.o \
           $(PBeFSSDR)/obj/SearchUtils.o \
           $(PBeFSSDR)/obj/Solution.o

PBeFSINC = -I$(PBeFSSDR)/include

# the data structures of A* [see include/PBeFS/astar/Structures.h]
ifeq ($(PBeFSSTRUCT),std)
PBeFSINC += -DPBEFS_STRUCTURES_STD
endif
ifeq ($(PBeFSSTRUCT),boost)
PBeFSINC += -DPBEFS_STRUCTURES_BOOST
endif
ifeq ($(PBeFSSTRUCT),ads)
PBeFSINC += -DPBEFS_STRUCTURES_ADS -I$(ADSSDR)/include
endif

PBeFSH   = $(PBeFSSDR)/include/PBeFS/PBeFS.h \
           $(PBeFSSDR)/include/PBeFS/astar.h \
           $(PBeFSSDR)/include/PBeFS/bnb.h \
           $(PBeFSSDR)/include/PBeFS/core.h \
           $(PBeFSSDR)/include/PBeFS/core/SearchUtils.h \
           $(PBeFSSDR)/include/PBeFS/core/ThreadMailbox.h \
           $(PBeFSSDR)/include/PBeFS/core/ThreadRuntime.h \
           $(PBeFSSDR)/include/PBeFS/astar/AStar.h \
           $(PBeFSSDR)/include/PBeFS/astar/Distribution.h \
           $(PBeFSSDR)/include/PBeFS/astar/Domain.h \
           $(PBeFSSDR)/include/PBeFS/astar/ExplicitGraph.h \
           $(PBeFSSDR)/include/PBeFS/astar/Grid.h \
           $(PBeFSSDR)/include/PBeFS/astar/HDAStar.h \
           $(PBeFSSDR)/include/PBeFS/astar/MissionariesCannibals.h \
           $(PBeFSSDR)/include/PBeFS/astar/NQueens.h \
           $(PBeFSSDR)/include/PBeFS/astar/Node.h \
           $(PBeFSSDR)/include/PBeFS/astar/NodeTable.h \
           $(PBeFSSDR)/include/PBeFS/astar/OpenList.h \
           $(PBeFSSDR)/include/PBeFS/astar/PASE.h \
           $(PBeFSSDR)/include/PBeFS/astar/Romania.h \
           $(PBeFSSDR)/include/PBeFS/astar/SPAStar.h \
           $(PBeFSSDR)/include/PBeFS/astar/SafePBNF.h \
           $(PBeFSSDR)/include/PBeFS/astar/SequenceScore.h \
           $(PBeFSSDR)/include/PBeFS/astar/SlidingTile.h \
           $(PBeFSSDR)/include/PBeFS/astar/SlowDomain.h \
           $(PBeFSSDR)/include/PBeFS/astar/Solution.h \
           $(PBeFSSDR)/include/PBeFS/astar/Structures.h \
           $(PBeFSSDR)/include/PBeFS/astar/TravelingSalesman.h \
           $(PBeFSSDR)/include/PBeFS/astar/VacuumWorld.h \
           $(PBeFSSDR)/include/PBeFS/astar/WaterJugs.h \
           $(wildcard $(PBeFSSDR)/include/PBeFS/astar/structures/*/*.h) \
           $(PBeFSSDR)/include/PBeFS/bnb/BranchAndBound.h \
           $(PBeFSSDR)/include/PBeFS/bnb/BranchAndBoundScheme.h \
           $(PBeFSSDR)/include/PBeFS/bnb/DeterministicBranchAndBound.h \
           $(PBeFSSDR)/include/PBeFS/bnb/Knapsack.h \
           $(PBeFSSDR)/include/PBeFS/bnb/ParallelBranchAndBound.h \
           $(PBeFSSDR)/include/PBeFS/bnb/SelfSplitBranchAndBound.h \
           $(PBeFSSDR)/include/PBeFS/bnb/StealingBranchAndBound.h

# the headers of the optional runtimes- - - - - - - - - - - - - - - - - - - -

PBeFSOMPH = $(PBeFSSDR)/include/PBeFS/core/OpenMPRuntime.h
PBeFSFFH  = $(PBeFSSDR)/include/PBeFS/core/FastFlowRuntime.h \
            $(PBeFSSDR)/include/PBeFS/astar/FarmAStar.h \
            $(PBeFSSDR)/include/PBeFS/bnb/FarmBranchAndBound.h
PBeFSMPIH = $(PBeFSSDR)/include/PBeFS/astar/HDAStarMPI.h \
            $(PBeFSSDR)/include/PBeFS/bnb/MPIBranchAndBound.h \
            $(PBeFSSDR)/include/PBeFS/core/MPIMailbox.h

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(PBeFSOBJ) $(PBeFSSDR)/*~ $(PBeFSSDR)/include/PBeFS/*~ \
	      $(PBeFSSDR)/include/PBeFS/*/*~ $(PBeFSSDR)/src/*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(PBeFSSDR)/obj/Distribution.o: $(PBeFSSDR)/src/Distribution.cpp $(PBeFSH)
	$(CC) -c $(PBeFSSDR)/src/Distribution.cpp -o $@ $(PBeFSINC) $(SW)

$(PBeFSSDR)/obj/Grid.o: $(PBeFSSDR)/src/Grid.cpp $(PBeFSH)
	$(CC) -c $(PBeFSSDR)/src/Grid.cpp -o $@ $(PBeFSINC) $(SW)

$(PBeFSSDR)/obj/SearchUtils.o: $(PBeFSSDR)/src/SearchUtils.cpp $(PBeFSH)
	$(CC) -c $(PBeFSSDR)/src/SearchUtils.cpp -o $@ $(PBeFSINC) $(SW)

$(PBeFSSDR)/obj/Solution.o: $(PBeFSSDR)/src/Solution.cpp $(PBeFSH)
	$(CC) -c $(PBeFSSDR)/src/Solution.cpp -o $@ $(PBeFSINC) $(SW)

############################ End of makefile #################################
