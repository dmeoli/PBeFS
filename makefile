###############################################################################
################################## makefile ###################################
###############################################################################
#                                                                             #
#   makefile of PBeFS                                                         #
#                                                                             #
#   The searches of PBeFS are templates over the Domain, and they are         #
#   in the headers; the objects hold the parts that do not depend on it.      #
#                                                                             #
#   Input:  $(CC)      = compiler command                                     #
#           $(SW)      = compiler options                                     #
#           $(PBeFSSDR)  = the directory where PBeFS is                       #
#           $(PBeFSSTRUCT) = the data structures of A*: std or boost, none    #
#                          (or ads) for those of ADS                          #
#           $(ADSSDR)    = the directory where ADS is                         #
#                                                                             #
#   Output: $(PBeFSOBJ)  = the final object(s) / library                      #
#           $(PBeFSH)    = the .h files to include                            #
#           $(PBeFSINC)  = the -I$( include directory )                       #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Università di Pisa                              #
#                                                                             #
###############################################################################

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

PBeFSOBJ = $(PBeFSSDR)/obj/Grid.o \
           $(PBeFSSDR)/obj/Planning.o \
           $(PBeFSSDR)/obj/SASHeuristics.o \
           $(PBeFSSDR)/obj/SASTask.o \
           $(PBeFSSDR)/obj/Solution.o

PBeFSINC = -I$(PBeFSSDR)/include -isystem $(ADSSDR)/include

# the data structures of A* [see include/PBeFS/astar/Structures.h]
ifeq ($(PBeFSSTRUCT),std)
PBeFSINC += -DPBEFS_STRUCTURES_STD
endif
ifeq ($(PBeFSSTRUCT),boost)
PBeFSINC += -DPBEFS_STRUCTURES_BOOST
endif

PBeFSH   = $(PBeFSSDR)/include/PBeFS/PBeFS.h \
           $(PBeFSSDR)/include/PBeFS/astar.h \
           $(PBeFSSDR)/include/PBeFS/bnb.h \
           $(PBeFSSDR)/include/PBeFS/core.h \
           $(PBeFSSDR)/include/PBeFS/core/SearchUtils.h \
           $(PBeFSSDR)/include/PBeFS/astar/AbstractionPDB.h \
           $(PBeFSSDR)/include/PBeFS/astar/AnytimeAStar.h \
           $(PBeFSSDR)/include/PBeFS/astar/AStar.h \
           $(PBeFSSDR)/include/PBeFS/astar/AutoDistribution.h \
           $(PBeFSSDR)/include/PBeFS/astar/BatchAStar.h \
           $(PBeFSSDR)/include/PBeFS/astar/BeamSearch.h \
           $(PBeFSSDR)/include/PBeFS/astar/BreadthFirst.h \
           $(PBeFSSDR)/include/PBeFS/astar/ConcurrentClosed.h \
           $(PBeFSSDR)/include/PBeFS/astar/Deepening.h \
           $(PBeFSSDR)/include/PBeFS/astar/DimacsGraph.h \
           $(PBeFSSDR)/include/PBeFS/astar/Distribution.h \
           $(PBeFSSDR)/include/PBeFS/astar/Domain.h \
           $(PBeFSSDR)/include/PBeFS/astar/DynamicGraph.h \
           $(PBeFSSDR)/include/PBeFS/astar/EnforcedHillClimbing.h \
           $(PBeFSSDR)/include/PBeFS/astar/EPASE.h \
           $(PBeFSSDR)/include/PBeFS/astar/ExplicitGraph.h \
           $(PBeFSSDR)/include/PBeFS/astar/FarmAStar.h \
           $(PBeFSSDR)/include/PBeFS/astar/FocalSearch.h \
           $(PBeFSSDR)/include/PBeFS/astar/FringeSearch.h \
           $(PBeFSSDR)/include/PBeFS/astar/GreedySearch.h \
           $(PBeFSSDR)/include/PBeFS/astar/Grid.h \
           $(PBeFSSDR)/include/PBeFS/astar/HDAStar.h \
           $(PBeFSSDR)/include/PBeFS/astar/HDLPAStar.h \
           $(PBeFSSDR)/include/PBeFS/astar/HPAStar.h \
           $(PBeFSSDR)/include/PBeFS/astar/ImprovedOptimisticSearch.h \
           $(PBeFSSDR)/include/PBeFS/astar/IteratedSearch.h \
           $(PBeFSSDR)/include/PBeFS/astar/LPAStar.h \
           $(PBeFSSDR)/include/PBeFS/astar/MissionariesCannibals.h \
           $(PBeFSSDR)/include/PBeFS/astar/MOAStar.h \
           $(PBeFSSDR)/include/PBeFS/astar/MOStructures.h \
           $(PBeFSSDR)/include/PBeFS/astar/MPLP.h \
           $(PBeFSSDR)/include/PBeFS/astar/MQAStar.h \
           $(PBeFSSDR)/include/PBeFS/astar/MultiAStar.h \
           $(PBeFSSDR)/include/PBeFS/astar/MultiObjective.h \
           $(PBeFSSDR)/include/PBeFS/astar/MultiObjectiveGraph.h \
           $(PBeFSSDR)/include/PBeFS/astar/MultiObjectiveGrid.h \
           $(PBeFSSDR)/include/PBeFS/astar/NBlockMap.h \
           $(PBeFSSDR)/include/PBeFS/astar/Node.h \
           $(PBeFSSDR)/include/PBeFS/astar/NodeTable.h \
           $(PBeFSSDR)/include/PBeFS/astar/NQueens.h \
           $(PBeFSSDR)/include/PBeFS/astar/OpenList.h \
           $(PBeFSSDR)/include/PBeFS/astar/OptimisticSearch.h \
           $(PBeFSSDR)/include/PBeFS/astar/Pancake.h \
           $(PBeFSSDR)/include/PBeFS/astar/ParallelGreedy.h \
           $(PBeFSSDR)/include/PBeFS/astar/ParallelIDAStar.h \
           $(PBeFSSDR)/include/PBeFS/astar/ParallelMOAStar.h \
           $(PBeFSSDR)/include/PBeFS/astar/PartialExpansion.h \
           $(PBeFSSDR)/include/PBeFS/astar/PASE.h \
           $(PBeFSSDR)/include/PBeFS/astar/PSDD.h \
           $(PBeFSSDR)/include/PBeFS/astar/PatternDatabase.h \
           $(PBeFSSDR)/include/PBeFS/astar/Permutation.h \
           $(PBeFSSDR)/include/PBeFS/astar/PermutationPDB.h \
           $(PBeFSSDR)/include/PBeFS/astar/Planning.h \
           $(PBeFSSDR)/include/PBeFS/astar/PotentialSearch.h \
           $(PBeFSSDR)/include/PBeFS/astar/PWAStar.h \
           $(PBeFSSDR)/include/PBeFS/astar/RBFS.h \
           $(PBeFSSDR)/include/PBeFS/astar/RealTimeSearch.h \
           $(PBeFSSDR)/include/PBeFS/astar/RobotNav2D.h \
           $(PBeFSSDR)/include/PBeFS/astar/Romania.h \
           $(PBeFSSDR)/include/PBeFS/astar/RubiksCube.h \
           $(PBeFSSDR)/include/PBeFS/astar/SafePBNF.h \
           $(PBeFSSDR)/include/PBeFS/astar/SASHeuristics.h \
           $(PBeFSSDR)/include/PBeFS/astar/SASTask.h \
           $(PBeFSSDR)/include/PBeFS/astar/SequenceAlignment.h \
           $(PBeFSSDR)/include/PBeFS/astar/SequenceScore.h \
           $(PBeFSSDR)/include/PBeFS/astar/SlidingTile.h \
           $(PBeFSSDR)/include/PBeFS/astar/SlowDomain.h \
           $(PBeFSSDR)/include/PBeFS/astar/Solution.h \
           $(PBeFSSDR)/include/PBeFS/astar/SPAStar.h \
           $(PBeFSSDR)/include/PBeFS/astar/Structures.h \
           $(PBeFSSDR)/include/PBeFS/astar/TDS.h \
           $(PBeFSSDR)/include/PBeFS/astar/TopSpin.h \
           $(PBeFSSDR)/include/PBeFS/astar/TowersOfHanoi.h \
           $(PBeFSSDR)/include/PBeFS/astar/TravelingSalesman.h \
           $(PBeFSSDR)/include/PBeFS/astar/VacuumWorld.h \
           $(PBeFSSDR)/include/PBeFS/astar/WaterJugs.h \
           $(PBeFSSDR)/include/PBeFS/astar/Weight.h \
           $(PBeFSSDR)/include/PBeFS/astar/WidthSearch.h \
           $(wildcard $(PBeFSSDR)/include/PBeFS/astar/structures/*/*.h) \
           $(PBeFSSDR)/include/PBeFS/bnb/BranchAndBound.h \
           $(PBeFSSDR)/include/PBeFS/bnb/BranchAndBoundScheme.h \
           $(PBeFSSDR)/include/PBeFS/bnb/DeterministicBranchAndBound.h \
           $(PBeFSSDR)/include/PBeFS/bnb/FarmBranchAndBound.h \
           $(PBeFSSDR)/include/PBeFS/bnb/Knapsack.h \
           $(PBeFSSDR)/include/PBeFS/bnb/MultiQueueBranchAndBound.h \
           $(PBeFSSDR)/include/PBeFS/bnb/ParallelBranchAndBound.h \
           $(PBeFSSDR)/include/PBeFS/bnb/SelfSplitBranchAndBound.h \
           $(PBeFSSDR)/include/PBeFS/bnb/StealingBranchAndBound.h

# the headers of the optional runtimes- - - - - - - - - - - - - - - - - - - -

# what is compiled with mpicxx takes the C API of MPI only, without its C++
# bindings, which are deprecated and whose headers give warnings
PBeFSMPIDEF = -DPBEFS_HAS_MPI -DADS_WITH_MPI -DOMPI_SKIP_MPICXX \
              -DMPICH_SKIP_MPICXX

PBeFSMPIH = $(PBeFSSDR)/include/PBeFS/astar/HDAStarMPI.h \
            $(PBeFSSDR)/include/PBeFS/astar/ParallelIDAStarMPI.h \
            $(PBeFSSDR)/include/PBeFS/bnb/MPIBranchAndBound.h

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(PBeFSOBJ) $(PBeFSSDR)/*~ $(PBeFSSDR)/include/PBeFS/*~ \
	      $(PBeFSSDR)/include/PBeFS/*/*~ $(PBeFSSDR)/src/*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(PBeFSSDR)/obj/Grid.o: $(PBeFSSDR)/src/Grid.cpp $(PBeFSH)
	$(CC) -c $(PBeFSSDR)/src/Grid.cpp -o $@ $(PBeFSINC) $(SW)

$(PBeFSSDR)/obj/Planning.o: $(PBeFSSDR)/src/Planning.cpp $(PBeFSH)
	$(CC) -c $(PBeFSSDR)/src/Planning.cpp -o $@ $(PBeFSINC) $(SW)

$(PBeFSSDR)/obj/SASHeuristics.o: $(PBeFSSDR)/src/SASHeuristics.cpp $(PBeFSH)
	$(CC) -c $(PBeFSSDR)/src/SASHeuristics.cpp -o $@ $(PBeFSINC) $(SW)

$(PBeFSSDR)/obj/SASTask.o: $(PBeFSSDR)/src/SASTask.cpp $(PBeFSH)
	$(CC) -c $(PBeFSSDR)/src/SASTask.cpp -o $@ $(PBeFSINC) $(SW)

$(PBeFSSDR)/obj/Solution.o: $(PBeFSSDR)/src/Solution.cpp $(PBeFSH)
	$(CC) -c $(PBeFSSDR)/src/Solution.cpp -o $@ $(PBeFSINC) $(SW)

############################### End of makefile ###############################
