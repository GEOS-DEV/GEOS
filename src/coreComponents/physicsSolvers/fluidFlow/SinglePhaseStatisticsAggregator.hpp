/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhaseStatisticsAggregator.hpp
 * @details Region statistics data is stored as follow:

 * Problem : ProblemManager
 * |-> domain : DomainPartition
 *     |-> MeshBodies : Group
 *         |-> "cartesianMesh" : MeshBody
 *             |-> meshLevels : Group
 *                 |-> Level0 : MeshLevel
 *                 |   |-> ElementRegions : ElementRegionManager
 *                 |   |   |-> elementRegionsGroup : Group
 *                 |   |       |-> "Channel" : CellElementRegion
 *                 |   |       |   |-> elementSubRegions : Group
 *                 |   |       |       |-> "cb-0_0_0" : CellElementSubRegion
 *                 |   |       |       |   | * pressure : Wrapper< real64 array >
 *                 |   |       |       |   | * temperature : Wrapper< real64 array >
 *                 |   |       |       |   |-> sets : Group
 *                 |   |       |       |   |   | * all : Wrapper< index array >
 *                 |   |       |       |   |   | * myBox : Wrapper< index array >
 *                 |   |       |       |   |   [...] (other element sets)
 *                 |   |       |       |   |
 *                 |   |       |       |   [...] (other fields)
 *                 |   |       |       |
 *                 |   |       |       |-> "cb-0_0_1" : CellElementSubRegion
 *                 |   |       |       |   | * pressure : Wrapper< real64 array >
 *                 |   |       |       |   | * temperature : Wrapper< real64 array >
 *                 |   |       |       |   |-> sets : Group
 *                 |   |       |       |   |   | * all : Wrapper< index array >
 *                 |   |       |       |   |   | * myBox : Wrapper< index array >
 *                 |   |       |       |   |   [...] (other element sets)
 *                 |   |       |       |   |
 *                 |   |       |       |   [...] (other fields)
 *                 |   |       |       |
 *                 |   |       |       [...] (other sub-regions)
 *                 |   |       |
 *                 |   |       |-> "Barrier" : CellElementRegion
 *                 |   |           |-> "cb-1_0_0" : CellElementSubRegion
 *                 |   |           |-> "cb-1_0_1" : CellElementSubRegion
 *                 |   |           [...] (other sub-regions)
 *                 |   |
 *                 |   |-> nodeManager : NodeManager
 *                 |   |   |-> sets : Group
 *                 |   |       | * all : Wrapper< index array >
 *                 |   |       | * "myBox" : Wrapper< index array >
 *                 |   |       [...] (other element sets)
 *                 |   |
 *                 |   [...] (other element managers)
 *          ____   |   |
 *          |      |   |-> statistics : Group (storage for all stats)
 *          |      |       |-> "flowStats" : Group (storage for this instance stats)
 *          |      |       |   |-> setsStatistics : RegionStatistics (selected sets aggregate, mpi reduced)
 *          |      |       |       |-> "..." : RegionStatistics (set aggregate, mpi reduced)
 *          |      |       |       |   |       - all requested sets
 *          |      |       |       |   |       - "all" if no set restriction or if requested
 *          |      |       |       |   |       - "__compound" if a compoud set is needed to compute mesh-level stats
 *          |      |       |       |   |
 *          |      |       |       |   |-> "Channel" : RegionStatistics (region aggregate, mpi reduced)
 *          |      |       |       |   |   |-> "cb-0_0_0" : RegionStatistics (sub-region compute read-back)
 *  stats   |      |       |       |   |   |-> "cb-0_0_1" : RegionStatistics (sub-region compute read-back)
 *  data -> |      |       |       |   |   [...] (other sub-regions stats)
 *          |      |       |       |   |
 *          |      |       |       |   |-> "Barrier" : RegionStatistics (region aggregate, mpi reduced)
 *          |      |       |       |       |-> cb-1_0_0 : RegionStatistics (sub-region compute read-back)
 *          |      |       |       |       |-> cb-1_0_1 : RegionStatistics (sub-region compute read-back)
 *          |      |       |       |       [...] (other sub-regions stats)
 *          |      |       |       |
 *          |      |       |       [...] (other set aggregates)
 *          |      |       |
 *          |___   |       [...] (other stats storages)
 *                 |
 *                 [...] (other discretizations)
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASESTATISTICSAGGREGATOR_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASESTATISTICSAGGREGATOR_HPP_

#include "physicsSolvers/StatisticsAggregatorBase.hpp"

namespace geos
{

class SinglePhaseBase;

namespace singlePhaseStatistics
{
class StatsAggregator;
class RegionStatistics;
}

template<>
struct StatsAggregatorTraits< singlePhaseStatistics::StatsAggregator >
{
  using SolverType = SinglePhaseBase;
  using StatsGroupType = singlePhaseStatistics::RegionStatistics;
};

namespace singlePhaseStatistics
{

/**
 * @brief Output data group to contain the result of a given stat aggregator on the dataRepository.
 *        Attributes are public since the class is a POD.
 * @todo repair 1D HDF5 outputs to enable stats HDF5 outputs
 */
class RegionStatistics : public RegionStatisticsBase
{
public:

  /// Element Count
  globalIndex m_elemCount;

  /// average region pressure (numerator value before postAggregateCompute())
  real64 m_averagePressure;
  /// minimum region pressure
  real64 m_minPressure;
  /// maximum region pressure
  real64 m_maxPressure;

  /// minimum region delta pressure
  real64 m_minDeltaPressure;
  /// maximum region delta pressure
  real64 m_maxDeltaPressure;

  /// average region temperature (numerator value before postAggregateCompute())
  real64 m_averageTemperature;
  /// minimum region temperature
  real64 m_minTemperature;
  /// maximum region temperature
  real64 m_maxTemperature;

  /// fluid mass
  real64 m_totalMass;

  /// total region pore volume
  real64 m_totalDynamicPoreVolume;
  /// total region uncompacted pore volume (not necessarily output, useful for weighting cell pressure data)
  real64 m_totalUncompactedPoreVolume;

  // TODO? -> split to struct PressureStats...MassStats:
  // - optional computation of each stats
  // - VKS for struct name ("pressureStats"..."massStats")
  // - current RegionStatistics struct bits

  /**
   * @brief Construct a new Region Statistics object
   * @param targetName name of the data-repository object that is targeted by the statistics
   *                   (mesh level / region / sub-region).
   * @param parent the instance parent in data-repository
   */
  RegionStatistics( string const & targetName,
                    dataRepository::Group * const parent,
                    string_view setName,
                    bool statsOutputEnabled );

  RegionStatistics( RegionStatistics && ) = default;

};

/**
 * @brief Reponsible of computing physical statistics over the grid, registering the result in the
 *        data repository, but not storing / outputing it by itself. It does not have mutable state
 *        except the encountered issues.
 * @todo repair 1D HDF5 outputs to enable stats HDF5 outputs
 */
class StatsAggregator : public StatsAggregatorBase< StatsAggregator >
{
public:

  using Base = StatsAggregatorBase< StatsAggregator >;

  /**
   * @brief Construct a new Stats Aggregator object
   * @param ownerName the unique name of the entity requesting the statistics.
   *                  An error is thrown if not unique in this context.
   * @param meshBodies The Group containing the MeshBody objects
   * @param statsOutputEnabled If true, the stats are saved in the output HDF5
   *                           (through dataRepository::RestartFlags, but not functional for this output for now).
   */
  StatsAggregator( dataRepository::DataContext const & ownerDataContext,
                   dataRepository::Group & meshBodies,
                   bool statsOutputEnabled );

  /**
   * @brief Enable the computation of region statistics, initialize data structure to collect them.
   *        Register the resulting data wrappers so they will be targeted by TimeHistory output
   * @note Must be called in or after the "registerDataOnMesh" initialization phase
   * @param setNames The list of mesh element sets to restrict the statistics to.
   *                 If empty, the whole discretization is processed.
   */
  void enableRegionStatisticsAggregation( string_array const & setNames = string_array() );

  // template implementations
  /// @cond DO_NOT_DOCUMENT

  void initStats( RegionStatistics & stats, real64 time ) const;
  void computeSubRegionRankStats( CellElementSubRegion & subRegion,
                                  RegionStatistics & subRegionStats,
                                  SetType const & targetSet ) const;
  void aggregateStats( RegionStatistics & stats, RegionStatistics const & other ) const;
  void mpiAggregateStats( RegionStatistics & stats ) const;
  void postAggregateStats( RegionStatistics & stats );

  /// @endcond

};

} /* namespace singlePhaseStatistics */

} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASESTATISTICSAGGREGATOR_HPP_ */
