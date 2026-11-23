/*--------------------------------------------------------------------------*/
/*---------------------------- File PBeFS.cppm -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Interface of the C++20 module pbefs, which exports every name that the
 * headers of include/PBeFS declare in the namespace PBeFS_di_unipi_it and in
 * those nested in it, but for the namespaces detail; a translation unit can
 * then write
 *
 *     import pbefs;
 *
 * in place of the #include of the headers. The names that only some sets of
 * data structures of A* declare [see PBeFS/astar/Structures.h] are exported
 * with the set that the module is compiled with, those of the searches on MPI
 * only if it is compiled with ADS_WITH_MPI, as the target PBeFS::MPI does,
 * and those of the symbolic searches only if it is compiled with
 * PBEFS_HAS_SYLVAN, as the target PBeFS::Symbolic does. The macros are not
 * exported by a module, and those that configure the headers take the value
 * they have where the module is compiled.
 *
 * This file is written by cmake/pbefs_module_exports.py, which reads the
 * names from the AST of the headers: it is not to be edited by hand.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Università di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------- GLOBAL MODULE FRAGMENT -------------------------*/
/*--------------------------------------------------------------------------*/

module;

#include "PBeFS/PBeFS.h"

#if defined( ADS_WITH_MPI )
 #include "PBeFS/astar/HDAStarMPI.h"
 #include "PBeFS/astar/ParallelIDAStarMPI.h"
 #include "PBeFS/bnb/MPIHubBranchAndBound.h"
#endif

#if defined( PBEFS_HAS_SYLVAN )
 #include "PBeFS/symbolic.h"
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- MODULE PURVIEW -----------------------------*/
/*--------------------------------------------------------------------------*/

export module pbefs;

export namespace PBeFS_di_unipi_it
{

// astar/AStar.h
using PBeFS_di_unipi_it::AStar;
using PBeFS_di_unipi_it::astar;

// astar/AbstractionPDB.h
using PBeFS_di_unipi_it::AbstractSpace;
using PBeFS_di_unipi_it::AbstractionTable;
using PBeFS_di_unipi_it::abstract_distances;

// astar/AnytimeAStar.h
using PBeFS_di_unipi_it::AnytimeAStar;
using PBeFS_di_unipi_it::AnytimeOptions;
using PBeFS_di_unipi_it::anytime_astar;

// astar/AutoDistribution.h
using PBeFS_di_unipi_it::CutObjective;
using PBeFS_di_unipi_it::DistributionEstimate;
using PBeFS_di_unipi_it::FeatureGraph;
using PBeFS_di_unipi_it::FeatureModel;
using PBeFS_di_unipi_it::Gf2Basis;
using PBeFS_di_unipi_it::auto_distribution;
using PBeFS_di_unipi_it::auto_distribution_names;
using PBeFS_di_unipi_it::dynamic_state_abstraction;
using PBeFS_di_unipi_it::fluency_abstraction;
using PBeFS_di_unipi_it::greedy_abstraction;
using PBeFS_di_unipi_it::kernel_selection_table;
using PBeFS_di_unipi_it::kernel_table;
using PBeFS_di_unipi_it::make_distribution;
using PBeFS_di_unipi_it::operator_table;
using PBeFS_di_unipi_it::sparsest_cut_abstraction;
using PBeFS_di_unipi_it::state_abstraction;

// astar/BatchAStar.h
using PBeFS_di_unipi_it::BatchAStar;
using PBeFS_di_unipi_it::BatchOptions;
using PBeFS_di_unipi_it::batch_astar;

// astar/BeamSearch.h
using PBeFS_di_unipi_it::BeamOptions;
using PBeFS_di_unipi_it::BeamSearch;
using PBeFS_di_unipi_it::beam_search;

// astar/BiBestFirst.h
using PBeFS_di_unipi_it::BiBestFirst;
using PBeFS_di_unipi_it::BiPriority;
using PBeFS_di_unipi_it::bae;
using PBeFS_di_unipi_it::bidijkstra;
using PBeFS_di_unipi_it::fmm;
using PBeFS_di_unipi_it::mm;
using PBeFS_di_unipi_it::mme;

// astar/Bidirectional.h
using PBeFS_di_unipi_it::BiFronts;
using PBeFS_di_unipi_it::BiNode;
using PBeFS_di_unipi_it::BiQueue;
using PBeFS_di_unipi_it::BidirectionalOptions;
using PBeFS_di_unipi_it::BlindDomain;
using PBeFS_di_unipi_it::ReverseDomain;
using PBeFS_di_unipi_it::Side;
using PBeFS_di_unipi_it::SideChoice;
using PBeFS_di_unipi_it::backward_side;
using PBeFS_di_unipi_it::forward_side;

// astar/BreadthFirst.h
using PBeFS_di_unipi_it::BreadthFirst;
using PBeFS_di_unipi_it::breadth_first;

// astar/DDDAStar.h
using PBeFS_di_unipi_it::DDDAStar;
using PBeFS_di_unipi_it::DDDOptions;
using PBeFS_di_unipi_it::DiskDomain;
using PBeFS_di_unipi_it::ddd_astar;

// astar/DVCBS.h
using PBeFS_di_unipi_it::DVCBS;
using PBeFS_di_unipi_it::dvcbs;

// astar/Deepening.h
using PBeFS_di_unipi_it::CutoffHistogram;
using PBeFS_di_unipi_it::DeepeningOptions;

// astar/DimacsGraph.h
using PBeFS_di_unipi_it::DimacsCoordinates;
using PBeFS_di_unipi_it::dimacs_distance;
using PBeFS_di_unipi_it::dimacs_scale;
using PBeFS_di_unipi_it::read_dimacs_coordinates;
using PBeFS_di_unipi_it::read_dimacs_graph;
using PBeFS_di_unipi_it::set_dimacs_heuristic;

// astar/Distribution.h
using PBeFS_di_unipi_it::AbstractStateDistribution;
using PBeFS_di_unipi_it::AbstractZobristDistribution;
using PBeFS_di_unipi_it::FeatureAbstraction;
using PBeFS_di_unipi_it::HyperplaneDistribution;
using PBeFS_di_unipi_it::IncrementalDistribution;
using PBeFS_di_unipi_it::PerfectHashDistribution;
using PBeFS_di_unipi_it::RandomDistribution;
using PBeFS_di_unipi_it::StateHashDistribution;
using PBeFS_di_unipi_it::TableZobristDistribution;
using PBeFS_di_unipi_it::ZobristDistribution;
using PBeFS_di_unipi_it::distribution_key;
using PBeFS_di_unipi_it::key_owner;
using PBeFS_di_unipi_it::reduce;
using PBeFS_di_unipi_it::successor_owners;

// astar/Domain.h
using PBeFS_di_unipi_it::AbstractFeatureDomain;
using PBeFS_di_unipi_it::AbstractStateDomain;
using PBeFS_di_unipi_it::BackwardOperatorDomain;
using PBeFS_di_unipi_it::BatchHeuristicDomain;
using PBeFS_di_unipi_it::BidirectionalDomain;
using PBeFS_di_unipi_it::CostChanges;
using PBeFS_di_unipi_it::DeferredDomain;
using PBeFS_di_unipi_it::DistanceDomain;
using PBeFS_di_unipi_it::Domain;
using PBeFS_di_unipi_it::EdgeCostChanged;
using PBeFS_di_unipi_it::EdgeDomain;
using PBeFS_di_unipi_it::EventDomain;
using PBeFS_di_unipi_it::FeatureDomain;
using PBeFS_di_unipi_it::FeatureGroupDomain;
using PBeFS_di_unipi_it::FeatureTransitionDomain;
using PBeFS_di_unipi_it::GoalChanged;
using PBeFS_di_unipi_it::GoalFeatureDomain;
using PBeFS_di_unipi_it::HeterogeneousEdgeDomain;
using PBeFS_di_unipi_it::IncrementalFeatureDomain;
using PBeFS_di_unipi_it::LazyEdgeDomain;
using PBeFS_di_unipi_it::LevelDomain;
using PBeFS_di_unipi_it::MinCostDomain;
using PBeFS_di_unipi_it::MultiHeuristicDomain;
using PBeFS_di_unipi_it::NBlockDomain;
using PBeFS_di_unipi_it::OperatorDomain;
using PBeFS_di_unipi_it::OperatorOf;
using PBeFS_di_unipi_it::PairwiseDomain;
using PBeFS_di_unipi_it::PartialExpansionDomain;
using PBeFS_di_unipi_it::PerfectHashDomain;
using PBeFS_di_unipi_it::PredecessorDomain;
using PBeFS_di_unipi_it::PreferredDomain;
using PBeFS_di_unipi_it::RelevantFeatureDomain;
using PBeFS_di_unipi_it::TreeOperatorDomain;
using PBeFS_di_unipi_it::cost_changes_of;
using PBeFS_di_unipi_it::cost_t;
using PBeFS_di_unipi_it::for_each_child;
using PBeFS_di_unipi_it::for_each_parent;
using PBeFS_di_unipi_it::for_each_tree_child;
using PBeFS_di_unipi_it::heuristic_of;
using PBeFS_di_unipi_it::heuristics_of;
using PBeFS_di_unipi_it::no_abstract_feature;
using PBeFS_di_unipi_it::no_op;
using PBeFS_di_unipi_it::path_cost;
using PBeFS_di_unipi_it::state_t;

// astar/DynamicGraph.h
using PBeFS_di_unipi_it::DynamicGraph;

// astar/EPASE.h
using PBeFS_di_unipi_it::EPASE;
using PBeFS_di_unipi_it::EPASEOptions;

// astar/EnforcedHillClimbing.h
using PBeFS_di_unipi_it::EHCOptions;
using PBeFS_di_unipi_it::EnforcedHillClimbing;
using PBeFS_di_unipi_it::PreferredUsage;
using PBeFS_di_unipi_it::ehc;

// astar/ExplicitGraph.h
using PBeFS_di_unipi_it::ExplicitGraph;

// astar/FarmAStar.h
using PBeFS_di_unipi_it::FarmAStar;
using PBeFS_di_unipi_it::FarmOptions;

// astar/FocalSearch.h
using PBeFS_di_unipi_it::FocalOptions;
using PBeFS_di_unipi_it::FocalOrder;
using PBeFS_di_unipi_it::FocalSearch;
using PBeFS_di_unipi_it::LazyQueue;

// astar/FringeSearch.h
using PBeFS_di_unipi_it::FringeSearch;
using PBeFS_di_unipi_it::fringe_search;

// astar/GBFHS.h
using PBeFS_di_unipi_it::GBFHS;
using PBeFS_di_unipi_it::bibfs;
using PBeFS_di_unipi_it::gbfhs;

// astar/GreedySearch.h
using PBeFS_di_unipi_it::GreedyEvaluation;
using PBeFS_di_unipi_it::GreedyNode;
using PBeFS_di_unipi_it::GreedyOpen;
using PBeFS_di_unipi_it::GreedyOptions;
using PBeFS_di_unipi_it::GreedySearch;
using PBeFS_di_unipi_it::HeuristicCache;
using PBeFS_di_unipi_it::QueueKind;
using PBeFS_di_unipi_it::QueueSpec;
using PBeFS_di_unipi_it::TieBreak;
using PBeFS_di_unipi_it::gbfs;

// astar/Grid.h
using PBeFS_di_unipi_it::Grid;
using PBeFS_di_unipi_it::GridBoard;
using PBeFS_di_unipi_it::GridCost;
using PBeFS_di_unipi_it::GridMap;
using PBeFS_di_unipi_it::GridNBlocks;
using PBeFS_di_unipi_it::grid_sqrt2;

// astar/HDAStar.h
using PBeFS_di_unipi_it::HDAOptions;
using PBeFS_di_unipi_it::HDAStar;
using PBeFS_di_unipi_it::LockedChannelsOf;
using PBeFS_di_unipi_it::ParallelOptions;

// astar/HDAStarMPI.h
#if defined( ADS_WITH_MPI )
using PBeFS_di_unipi_it::HDAStarMPI;
#endif

// astar/HDLPAStar.h
using PBeFS_di_unipi_it::HDDStarLite;
using PBeFS_di_unipi_it::HDLPAOptions;
using PBeFS_di_unipi_it::HDLPAStar;
using PBeFS_di_unipi_it::IncrementalOptions;
using PBeFS_di_unipi_it::LPADefaultQueue;

// astar/HPAStar.h
using PBeFS_di_unipi_it::HPAOptions;
using PBeFS_di_unipi_it::HPAStar;

// astar/IDAStar.h
using PBeFS_di_unipi_it::IDAOptions;
using PBeFS_di_unipi_it::IDAStar;
using PBeFS_di_unipi_it::cost_bounded_dfs;
using PBeFS_di_unipi_it::idastar;

// astar/IDDAStar.h
using PBeFS_di_unipi_it::IDDAStar;
using PBeFS_di_unipi_it::IDDOptions;
using PBeFS_di_unipi_it::idd_astar;

// astar/ImprovedOptimisticSearch.h
using PBeFS_di_unipi_it::ImprovedOptimisticSearch;

// astar/IteratedSearch.h
using PBeFS_di_unipi_it::IteratedOptions;
using PBeFS_di_unipi_it::IteratedSearch;
using PBeFS_di_unipi_it::lama_phases;

// astar/LPAStar.h
using PBeFS_di_unipi_it::DStarLite;
using PBeFS_di_unipi_it::LPABuckets;
using PBeFS_di_unipi_it::LPAHeap;
using PBeFS_di_unipi_it::LPAKey;
using PBeFS_di_unipi_it::LPAStar;
using PBeFS_di_unipi_it::RobotDomain;
using PBeFS_di_unipi_it::add_cost;
using PBeFS_di_unipi_it::erase_heap_key;
using PBeFS_di_unipi_it::lpa_below;
using PBeFS_di_unipi_it::set_heap_key;

// astar/MOAStar.h
using PBeFS_di_unipi_it::AutoFrontier;
using PBeFS_di_unipi_it::EMOAFrontier;
using PBeFS_di_unipi_it::MOAStar;
using PBeFS_di_unipi_it::MOAlgorithm;
using PBeFS_di_unipi_it::NAMOAStarDR;
using PBeFS_di_unipi_it::add_costs;
using PBeFS_di_unipi_it::finish_front;
using PBeFS_di_unipi_it::mo_algorithm;
using PBeFS_di_unipi_it::mo_search;
using PBeFS_di_unipi_it::trace_labels;

// astar/MPLP.h
using PBeFS_di_unipi_it::MPLP;
using PBeFS_di_unipi_it::MPLPOptions;

// astar/MQAStar.h
using PBeFS_di_unipi_it::MQAStar;

// astar/MissionariesCannibals.h
using PBeFS_di_unipi_it::MissionariesCannibals;

// astar/MultiAStar.h
using PBeFS_di_unipi_it::MultiAStar;

// astar/MultiObjective.h
using PBeFS_di_unipi_it::ArrayFrontier;
using PBeFS_di_unipi_it::BlockIndex;
using PBeFS_di_unipi_it::LayerIndex;
using PBeFS_di_unipi_it::MOBucketQueue;
using PBeFS_di_unipi_it::MOSearchOptions;
using PBeFS_di_unipi_it::MOSearchStats;
using PBeFS_di_unipi_it::MultiObjectiveDomain;
using PBeFS_di_unipi_it::ParetoSolution;
using PBeFS_di_unipi_it::ScalarFrontier;
using PBeFS_di_unipi_it::TreeFrontier;
using PBeFS_di_unipi_it::dominates;
using PBeFS_di_unipi_it::mo_cost_t;
using PBeFS_di_unipi_it::mo_dead_end;
using PBeFS_di_unipi_it::mo_scalar_t;
using PBeFS_di_unipi_it::mo_state_t;
using PBeFS_di_unipi_it::objective_order;
using PBeFS_di_unipi_it::objectives_of;
using PBeFS_di_unipi_it::pareto_front;
using PBeFS_di_unipi_it::pareto_solution_t;
using PBeFS_di_unipi_it::permute;
using PBeFS_di_unipi_it::truncate;
using PBeFS_di_unipi_it::unpermute;
using PBeFS_di_unipi_it::weakly_dominates;

// astar/MultiObjectiveGraph.h
using PBeFS_di_unipi_it::ExtraCost;
using PBeFS_di_unipi_it::MultiObjectiveGraph;
using PBeFS_di_unipi_it::distances_to;
using PBeFS_di_unipi_it::extra_cost;
using PBeFS_di_unipi_it::read_dimacs_objectives;

// astar/MultiObjectiveGrid.h
using PBeFS_di_unipi_it::MultiObjectiveGrid;

// astar/NBS.h
using PBeFS_di_unipi_it::NBS;
using PBeFS_di_unipi_it::nbs;

// astar/NBlockMap.h
using PBeFS_di_unipi_it::NBlockMap;

// astar/NQueens.h
using PBeFS_di_unipi_it::NQueens;

// astar/Node.h
using PBeFS_di_unipi_it::Node;
using PBeFS_di_unipi_it::cheaper;
using PBeFS_di_unipi_it::cost_slack;
using PBeFS_di_unipi_it::infinite_cost;
using PBeFS_di_unipi_it::packed_fg;
using PBeFS_di_unipi_it::real_cost_tolerance;
using PBeFS_di_unipi_it::trace_path;

// astar/OpenList.h
using PBeFS_di_unipi_it::LowGOpenList;
using PBeFS_di_unipi_it::OpenEntry;
using PBeFS_di_unipi_it::OpenKind;

// astar/OptimisticSearch.h
using PBeFS_di_unipi_it::OptimisticOptions;
using PBeFS_di_unipi_it::OptimisticSearch;
using PBeFS_di_unipi_it::optimistic_search;

// astar/PASE.h
using PBeFS_di_unipi_it::PASE;
using PBeFS_di_unipi_it::PASEOptions;

// astar/PEMSearch.h
using PBeFS_di_unipi_it::PEMKind;
using PBeFS_di_unipi_it::PEMOptions;
using PBeFS_di_unipi_it::PEMSearch;

// astar/PSDD.h
using PBeFS_di_unipi_it::BFPSDD;
using PBeFS_di_unipi_it::DynamicBoundedPSDD;
using PBeFS_di_unipi_it::IDPSDD;
using PBeFS_di_unipi_it::LayeredNBlockSearch;
using PBeFS_di_unipi_it::PSDD;
using PBeFS_di_unipi_it::PSDDOptions;

// astar/PWAStar.h
using PBeFS_di_unipi_it::PWAStar;
using PBeFS_di_unipi_it::PWAStarOptions;

// astar/Pancake.h
using PBeFS_di_unipi_it::BurnedPancakePuzzle;
using PBeFS_di_unipi_it::PancakePuzzle;

// astar/ParallelBidirectional.h
using PBeFS_di_unipi_it::ParallelBiKind;
using PBeFS_di_unipi_it::ParallelBidirectional;
using PBeFS_di_unipi_it::pnba;

// astar/ParallelGreedy.h
using PBeFS_di_unipi_it::GreedyRule;
using PBeFS_di_unipi_it::GreedyTrace;
using PBeFS_di_unipi_it::HDGBFS;
using PBeFS_di_unipi_it::ParallelGreedy;
using PBeFS_di_unipi_it::ParallelGreedyOptions;

// astar/ParallelIDAStar.h
using PBeFS_di_unipi_it::FrontierRule;
using PBeFS_di_unipi_it::PIDAOptions;
using PBeFS_di_unipi_it::ParallelIDAStar;
using PBeFS_di_unipi_it::UnitSchedule;
using PBeFS_di_unipi_it::parallel_idastar;

// astar/ParallelIDAStarMPI.h
#if defined( ADS_WITH_MPI )
using PBeFS_di_unipi_it::ParallelIDAStarMPI;
#endif

// astar/ParallelMOAStar.h
using PBeFS_di_unipi_it::HDMOAStar;
using PBeFS_di_unipi_it::MOAStarPar;
using PBeFS_di_unipi_it::ParallelMOOptions;
using PBeFS_di_unipi_it::objective_orders;

// astar/PartialExpansion.h
using PBeFS_di_unipi_it::EPEAStar;
using PBeFS_di_unipi_it::PEAStar;
using PBeFS_di_unipi_it::PartialDelta;
using PBeFS_di_unipi_it::PartialExpansionBase;
using PBeFS_di_unipi_it::PartialExtra;
using PBeFS_di_unipi_it::PartialNoExtra;
using PBeFS_di_unipi_it::PartialOptions;
using PBeFS_di_unipi_it::epea_star;
using PBeFS_di_unipi_it::pea_star;

// astar/PatternDatabase.h
using PBeFS_di_unipi_it::AdditivePDB;
using PBeFS_di_unipi_it::BlankRule;
using PBeFS_di_unipi_it::ByteEntries;
using PBeFS_di_unipi_it::MinNibbleEntries;
using PBeFS_di_unipi_it::Mod3Entries;
using PBeFS_di_unipi_it::NibbleEntries;
using PBeFS_di_unipi_it::PDBOptions;
using PBeFS_di_unipi_it::ParityBitEntries;
using PBeFS_di_unipi_it::Partition;
using PBeFS_di_unipi_it::PatternTable;
using PBeFS_di_unipi_it::SlidingTilePDB;
using PBeFS_di_unipi_it::SparseEntries;
using PBeFS_di_unipi_it::TilePattern;
using PBeFS_di_unipi_it::partition_15_5_5_5;
using PBeFS_di_unipi_it::partition_15_6_6_3;
using PBeFS_di_unipi_it::partition_15_7_8;
using PBeFS_di_unipi_it::partition_24_6_6_6_6;

// astar/Permutation.h
using PBeFS_di_unipi_it::LexRanking;
using PBeFS_di_unipi_it::MyrvoldRuskeyRanking;
using PBeFS_di_unipi_it::Ranking;
using PBeFS_di_unipi_it::TreeRanking;
using PBeFS_di_unipi_it::bytes_hash;
using PBeFS_di_unipi_it::k_permutations;

// astar/PermutationPDB.h
using PBeFS_di_unipi_it::PDBDomain;
using PBeFS_di_unipi_it::PermutationDomain;
using PBeFS_di_unipi_it::PermutationGroup;
using PBeFS_di_unipi_it::PermutationHeuristic;
using PBeFS_di_unipi_it::PermutationMove;
using PBeFS_di_unipi_it::PermutationPDB;
using PBeFS_di_unipi_it::PermutationSpace;

// astar/Planning.h
using PBeFS_di_unipi_it::Planning;
using PBeFS_di_unipi_it::PlanningBase;
using PBeFS_di_unipi_it::PlanningHeuristic;
using PBeFS_di_unipi_it::PlanningOptions;
using PBeFS_di_unipi_it::planning_heuristic;

// astar/PotentialSearch.h
using PBeFS_di_unipi_it::PotentialSearch;

// astar/RBFS.h
using PBeFS_di_unipi_it::RBFS;
using PBeFS_di_unipi_it::rbfs;

// astar/RealTimeSearch.h
using PBeFS_di_unipi_it::RealTimeOptions;
using PBeFS_di_unipi_it::RealTimeSearch;
using PBeFS_di_unipi_it::lss_lrta_star;

// astar/RobotNav2D.h
using PBeFS_di_unipi_it::NavMap;
using PBeFS_di_unipi_it::NavOptions;
using PBeFS_di_unipi_it::RobotNav2D;

// astar/Romania.h
using PBeFS_di_unipi_it::Romania;

// astar/RubiksCube.h
using PBeFS_di_unipi_it::RubikPDBHeuristic;
using PBeFS_di_unipi_it::RubiksCube;

// astar/SASHeuristics.h
using PBeFS_di_unipi_it::LMCut;
using PBeFS_di_unipi_it::RelaxedHeuristic;
using PBeFS_di_unipi_it::SASAdditiveCartesian;
using PBeFS_di_unipi_it::SASCartesianAbstraction;
using PBeFS_di_unipi_it::SASMergeAndShrink;
using PBeFS_di_unipi_it::SASPatternDatabase;
using PBeFS_di_unipi_it::dead_end;

// astar/SASTask.h
using PBeFS_di_unipi_it::SASTask;
using PBeFS_di_unipi_it::StripsTask;

// astar/SPAStar.h
using PBeFS_di_unipi_it::SPAStar;

// astar/SafePBNF.h
using PBeFS_di_unipi_it::PBNFOptions;
using PBeFS_di_unipi_it::SafePBNF;

// astar/SequenceAlignment.h
using PBeFS_di_unipi_it::ScoreMatrix;
using PBeFS_di_unipi_it::SequenceAlignment;
using PBeFS_di_unipi_it::amino_count;
using PBeFS_di_unipi_it::amino_gap;
using PBeFS_di_unipi_it::amino_symbols;
using PBeFS_di_unipi_it::pam250;

// astar/SequenceScore.h
using PBeFS_di_unipi_it::SequenceScore;

// astar/SlidingTile.h
using PBeFS_di_unipi_it::EightPuzzle;
using PBeFS_di_unipi_it::FifteenPuzzle;
using PBeFS_di_unipi_it::SlidingTile;
using PBeFS_di_unipi_it::TileRegions;
using PBeFS_di_unipi_it::TwentyFourPuzzle;
using PBeFS_di_unipi_it::Word128;

// astar/SlowDomain.h
using PBeFS_di_unipi_it::EdgeWork;
using PBeFS_di_unipi_it::SlowDomain;
using PBeFS_di_unipi_it::SlowEdgeDomain;

// astar/Solution.h
using PBeFS_di_unipi_it::SearchOptions;
using PBeFS_di_unipi_it::SearchStats;
using PBeFS_di_unipi_it::Solution;
using PBeFS_di_unipi_it::Status;

// astar/TDS.h
using PBeFS_di_unipi_it::TDS;
using PBeFS_di_unipi_it::TDSOptions;

// astar/TopSpin.h
using PBeFS_di_unipi_it::TopSpin;

// astar/TowersOfHanoi.h
using PBeFS_di_unipi_it::HanoiPDB;
using PBeFS_di_unipi_it::HanoiSpace;
using PBeFS_di_unipi_it::TowersOfHanoi;
using PBeFS_di_unipi_it::hanoi_tops;

// astar/TravelingSalesman.h
using PBeFS_di_unipi_it::TravelingSalesman;

// astar/VacuumWorld.h
using PBeFS_di_unipi_it::VacuumWorld;

// astar/WaterJugs.h
using PBeFS_di_unipi_it::WaterJugs;

// astar/Weight.h
using PBeFS_di_unipi_it::UnitWeight;
using PBeFS_di_unipi_it::Weight;
using PBeFS_di_unipi_it::finish_solution;
using PBeFS_di_unipi_it::settle_improvements;

// astar/WidthSearch.h
using PBeFS_di_unipi_it::BestFirstWidth;
using PBeFS_di_unipi_it::IteratedWidth;
using PBeFS_di_unipi_it::NoveltyTable;
using PBeFS_di_unipi_it::SerializedWidth;
using PBeFS_di_unipi_it::WidthOptions;
using PBeFS_di_unipi_it::added_features;
using PBeFS_di_unipi_it::bfws;
using PBeFS_di_unipi_it::features_of;
using PBeFS_di_unipi_it::iw;
using PBeFS_di_unipi_it::siw;

// astar/structures/ads/ConcurrentClosed.h
using PBeFS_di_unipi_it::ConcurrentClosed;

// astar/structures/ads/ExplicitAdjacency.h
using PBeFS_di_unipi_it::ExplicitAdjacency;

// astar/structures/ads/GridAdjacency.h
using PBeFS_di_unipi_it::GridAdjacency;

// astar/structures/ads/MOStructures.h
using PBeFS_di_unipi_it::LexSet;
using PBeFS_di_unipi_it::MOHeap;

// astar/structures/ads/NodePool.h
using PBeFS_di_unipi_it::NodePool;

// astar/structures/ads/NodeTable.h
using PBeFS_di_unipi_it::NodeTable;

// astar/structures/ads/OpenList.h
using PBeFS_di_unipi_it::BucketOpenList;
using PBeFS_di_unipi_it::HeapOpenList;

// bnb/BBCheckpoint.h
using PBeFS_di_unipi_it::BBCheckpoint;
using PBeFS_di_unipi_it::SerializableBranchAndBound;
using PBeFS_di_unipi_it::bb_read_checkpoint;
using PBeFS_di_unipi_it::bb_write_checkpoint;
using PBeFS_di_unipi_it::checkpointed_branch_and_bound;

// bnb/BranchAndBound.h
using PBeFS_di_unipi_it::BBChild;
using PBeFS_di_unipi_it::BBCursor;
using PBeFS_di_unipi_it::BBEvaluator;
using PBeFS_di_unipi_it::BBLeaf;
using PBeFS_di_unipi_it::BBLeafKind;
using PBeFS_di_unipi_it::BBMode;
using PBeFS_di_unipi_it::BBOptions;
using PBeFS_di_unipi_it::BBQueue;
using PBeFS_di_unipi_it::BBSolution;
using PBeFS_di_unipi_it::BBStats;
using PBeFS_di_unipi_it::BBStatus;
using PBeFS_di_unipi_it::BBValue;
using PBeFS_di_unipi_it::BBValueTraits;
using PBeFS_di_unipi_it::BasicBBOptions;
using PBeFS_di_unipi_it::BasicEvaluation;
using PBeFS_di_unipi_it::BranchAndBound;
using PBeFS_di_unipi_it::BranchAndBoundSearch;
using PBeFS_di_unipi_it::EstimatingBranchAndBound;
using PBeFS_di_unipi_it::Evaluation;
using PBeFS_di_unipi_it::HashedBranchAndBound;
using PBeFS_di_unipi_it::IndexedBBEvaluator;
using PBeFS_di_unipi_it::LazyBBEvaluator;
using PBeFS_di_unipi_it::Selection;
using PBeFS_di_unipi_it::Sense;
using PBeFS_di_unipi_it::bb_value_t;
using PBeFS_di_unipi_it::branch_and_bound;

// bnb/BranchAndBoundScheme.h
using PBeFS_di_unipi_it::BBScheme;
using PBeFS_di_unipi_it::bb_scheme;
using PBeFS_di_unipi_it::bb_scheme_name;
using PBeFS_di_unipi_it::parallel_branch_and_bound;

// bnb/DeterministicBranchAndBound.h
using PBeFS_di_unipi_it::DeterministicBranchAndBound;

// bnb/FarmBranchAndBound.h
using PBeFS_di_unipi_it::FarmBranchAndBound;

// bnb/Flowshop.h
using PBeFS_di_unipi_it::Flowshop;

// bnb/IntervalBranchAndBound.h
using PBeFS_di_unipi_it::IntervalBranchAndBound;
using PBeFS_di_unipi_it::RadixBranchAndBound;

// bnb/Knapsack.h
using PBeFS_di_unipi_it::BasicKnapsack;
using PBeFS_di_unipi_it::Knapsack;

// bnb/MPIBranchAndBound.h
#if defined( ADS_WITH_MPI )
using PBeFS_di_unipi_it::MPIBranchAndBound;
using PBeFS_di_unipi_it::MPIBranchAndBoundSearch;
#endif

// bnb/MPIHubBranchAndBound.h
#if defined( ADS_WITH_MPI )
using PBeFS_di_unipi_it::BBHubParams;
using PBeFS_di_unipi_it::MPIHubBranchAndBoundSearch;
#endif

// bnb/MaxClique.h
using PBeFS_di_unipi_it::BitGraph;
using PBeFS_di_unipi_it::BitSet;
using PBeFS_di_unipi_it::MaxClique;

// bnb/MaxCommonSubgraph.h
using PBeFS_di_unipi_it::MaxCommonSubgraph;

// bnb/MultiQueueBranchAndBound.h
using PBeFS_di_unipi_it::MultiQueueBranchAndBound;
using PBeFS_di_unipi_it::MultiQueueValue;

// bnb/ParallelBranchAndBound.h
using PBeFS_di_unipi_it::BBIncumbent;
using PBeFS_di_unipi_it::BasicParallelBBOptions;
using PBeFS_di_unipi_it::DistributedBranchAndBound;
using PBeFS_di_unipi_it::ParallelBBOptions;
using PBeFS_di_unipi_it::SharedBranchAndBound;

// bnb/QAP.h
using PBeFS_di_unipi_it::QAP;
using PBeFS_di_unipi_it::linear_assignment;

// bnb/RacingBranchAndBound.h
using PBeFS_di_unipi_it::racing_branch_and_bound;

// bnb/SelfSplitBranchAndBound.h
using PBeFS_di_unipi_it::SelfSplitBranchAndBound;

// bnb/SkeletonBranchAndBound.h
using PBeFS_di_unipi_it::BBRule;
using PBeFS_di_unipi_it::SkeletonBranchAndBound;

// bnb/StealingBranchAndBound.h
using PBeFS_di_unipi_it::StealingBranchAndBound;

// bnb/SubgraphIsomorphism.h
using PBeFS_di_unipi_it::SubgraphIsomorphism;

// bnb/UTS.h
using PBeFS_di_unipi_it::SHA1;
using PBeFS_di_unipi_it::UTS;

// core/SearchUtils.h
using PBeFS_di_unipi_it::ParallelParams;
using PBeFS_di_unipi_it::Stopwatch;
using PBeFS_di_unipi_it::cache_line;
using PBeFS_di_unipi_it::pin_worker;
using PBeFS_di_unipi_it::prefetch;
using PBeFS_di_unipi_it::random_below;
using PBeFS_di_unipi_it::random_unit;
using PBeFS_di_unipi_it::run_workers;

// symbolic/Bdd.h
#if defined( PBEFS_HAS_SYLVAN )
using PBeFS_di_unipi_it::Bdd;
using PBeFS_di_unipi_it::BddRuntime;
using PBeFS_di_unipi_it::disjoin;
using PBeFS_di_unipi_it::image;
using PBeFS_di_unipi_it::image_into;
using PBeFS_di_unipi_it::images;
using PBeFS_di_unipi_it::nodes;
using PBeFS_di_unipi_it::preimage;
#endif

// symbolic/SymbolicDomain.h
#if defined( PBEFS_HAS_SYLVAN )
using PBeFS_di_unipi_it::DecodableSymbolicDomain;
using PBeFS_di_unipi_it::SymbolicAxiom;
using PBeFS_di_unipi_it::SymbolicAxiomDomain;
using PBeFS_di_unipi_it::SymbolicBackwardHeuristicDomain;
using PBeFS_di_unipi_it::SymbolicDomain;
using PBeFS_di_unipi_it::SymbolicEffect;
using PBeFS_di_unipi_it::SymbolicFact;
using PBeFS_di_unipi_it::SymbolicHeuristicDomain;
using PBeFS_di_unipi_it::SymbolicOperator;
using PBeFS_di_unipi_it::SymbolicTerm;
using PBeFS_di_unipi_it::infinite_term;
#endif

// symbolic/SymbolicModel.h
#if defined( PBEFS_HAS_SYLVAN )
using PBeFS_di_unipi_it::SymbolicLayout;
using PBeFS_di_unipi_it::SymbolicModel;
#endif

// symbolic/SymbolicPuzzles.h
#if defined( PBEFS_HAS_SYLVAN )
using PBeFS_di_unipi_it::SymbolicPancake;
using PBeFS_di_unipi_it::SymbolicSlidingTile;
#endif

// symbolic/SymbolicSAS.h
#if defined( PBEFS_HAS_SYLVAN )
using PBeFS_di_unipi_it::SymbolicSAS;
#endif

// symbolic/SymbolicSearch.h
#if defined( PBEFS_HAS_SYLVAN )
using PBeFS_di_unipi_it::SymbolicDirection;
using PBeFS_di_unipi_it::SymbolicImages;
using PBeFS_di_unipi_it::SymbolicOptions;
using PBeFS_di_unipi_it::SymbolicReachable;
using PBeFS_di_unipi_it::SymbolicSearch;
using PBeFS_di_unipi_it::SymbolicSolution;
using PBeFS_di_unipi_it::SymbolicStats;
using PBeFS_di_unipi_it::symbolic_bfs;
using PBeFS_di_unipi_it::symbolic_reachable;
using PBeFS_di_unipi_it::symbolic_search;
#endif

}

/*--------------------------------------------------------------------------*/
/*-------------------------- End File PBeFS.cppm ---------------------------*/
/*--------------------------------------------------------------------------*/
