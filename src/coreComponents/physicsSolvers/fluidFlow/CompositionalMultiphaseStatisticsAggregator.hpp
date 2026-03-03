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
 * @file CompositionalMultiphaseStatisticsAggregator.hpp
 * @details Region statistics data is stored as follow:

 * Problem : ProblemManager
 * |-> domain : DomainPartition
 *     |-> MeshBodies : Group
 *         |-> cartesianMesh : MeshBody
 *             |-> meshLevels : Group
 *                 |-> Level0 : MeshLevel
 *                 |   |-> nodeManager : NodeManager
 *                 |   |   |-> sets : Group
 *                 |   |       | * all : Wrapper< index array >
 *                 |   |       | * xneg : Wrapper< index array >
 *                 |   |       [...] (other element sets)
 *                 |   |
 *                 |   |-> ElementRegions : ElementRegionManager
 *                 |   |   |-> Channel : CellElementRegion
 *                 |   |   |   |-> cb-0_0_0 : CellElementSubRegion
 *                 |   |   |   |   | * pressure : Wrapper< real64 array >
 *                 |   |   |   |   | * temperature : Wrapper< real64 array >
 *                 |   |   |   |   [...] (other fields)
 *                 |   |   |   |
 *                 |   |   |   |-> cb-0_0_1 : CellElementSubRegion
 *                 |   |   |   |   | * pressure : Wrapper< real64 array >
 *                 |   |   |   |   | * temperature : Wrapper< real64 array >
 *                 |   |   |   |   [...] (other fields)
 *                 |   |   |   |
 *                 |   |   |   [...] (other sub-regions)
 *                 |   |   |
 *                 |   |   |-> Barrier : CellElementRegion
 *                 |   |       |-> cb-1_0_0 : CellElementSubRegion
 *                 |   |       |-> cb-1_0_1 : CellElementSubRegion
 *                 |   |       [...] (other sub-regions)
 *                 |   |
 *                 |   [...] (other element managers)
 *          ____   |   |
 *          |      |   |-> statistics : Group (storage for all stats)
 *          |      |       |-> compFlowStats : Group (storage for this instance stats)
 *          |      |       |   |-> cflStatistics : CFLStatistics
 *          |      |       |   |-> regionsStatistics : RegionStatistics (aggregate)
 *          |      |       |       |-> Channel : RegionStatistics (aggregate, mpi reduced)
 *          |      |       |       |   |-> cb-0_0_0 : RegionStatistics (compute read-back)
 *  stats   |      |       |       |   |-> cb-0_0_1 : RegionStatistics (compute read-back)
 *  data -> |      |       |       |   [...] (other sub-regions stats)
 *          |      |       |       |
 *          |      |       |       |-> Barrier : RegionStatistics (aggregate, mpi reduced)
 *          |      |       |           |-> cb-1_0_0 : RegionStatistics (compute read-back)
 *          |      |       |           |-> cb-1_0_1 : RegionStatistics (compute read-back)
 *          |      |       |           [...] (other sub-regions stats)
 *          |      |       |
 *          |___   |       [...] (other stats storages)
 *                 |
 *                 [...] (other discretizations)
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASESTATISTICSAGGREGATOR_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASESTATISTICSAGGREGATOR_HPP_

#include "common/DataTypes.hpp"
#include "physicsSolvers/StatisticsAggregatorBase.hpp"

namespace geos
{

class CompositionalMultiphaseBase;

namespace compositionalMultiphaseStatistics
{
class StatsAggregator;
class RegionStatistics;
}

template<>
struct StatsAggregatorTraits< compositionalMultiphaseStatistics::StatsAggregator >
{
  using SolverType = CompositionalMultiphaseBase;
  using StatsGroupType = compositionalMultiphaseStatistics::RegionStatistics;
};

namespace compositionalMultiphaseStatistics
{

struct AggregatorParameters
{
  /// Threshold to decide whether a phase is considered "mobile" or not
  real64 m_relpermThreshold;

  // TODO: add other params like views and stuff
};

/**
 * @brief Output data group to contain the result of a given stat aggregator on the dataRepository.
 *        Attributes are public since the class is a POD.
 * @todo repair 1D HDF5 outputs to enable stats HDF5 outputs
 */
class RegionStatistics : public RegionStatisticsBase
{
public:

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

  /// total region pore volume
  real64 m_totalPoreVolume;
  /// total region uncompacted pore volume (not necessarily output, useful for weighting cell pressure data)
  real64 m_totalUncompactedPoreVolume;
  /// phase region dynamic pore volume
  array1d< real64 > const m_phaseDynamicPoreVolume;

  /// region phase mass (trapped and non-trapped, immobile and mobile)
  array1d< real64 > const m_phaseMass;
  /// trapped region phase mass
  array1d< real64 > const m_trappedPhaseMass;
  /// non-trapped region phase mass (available after postAggregateCompute())
  array1d< real64 > const m_nonTrappedPhaseMass;
  /// immobile region phase mass
  array1d< real64 > const m_immobilePhaseMass;
  /// mobile region phase mass (available after postAggregateCompute())
  array1d< real64 > const m_mobilePhaseMass;
  /// region component mass
  array2d< real64 > const m_componentMass;

  // TODO? -> split to struct PressureStats...MassStats:
  // - optional computation of each stats
  // - VKS for struct name ("pressureStats"..."massStats")
  // - current RegionStatistics struct bits

  /**
   * @brief Construct a new Region Statistics object
   * @param targetName name of the data-repository object that is targeted by the statistics
   *                   (mesh level / region / sub-region).
   * @param parent the instance parent in data-repository
   * @param numPhases Fluid phase count
   * @param numComponents Fluid component count
   */
  RegionStatistics( string const & targetName, dataRepository::Group * const parent,
                    integer numPhases, integer numComponents );

  RegionStatistics( RegionStatistics && ) = default;

};

/**
 * @brief Output data group to contain the result of a given stat aggregator on the dataRepository.
 *        Attributes are public since the class is a POD. Can it be replaced by a wrapped-struct?
 */
class CFLStatistics : public RegionStatisticsBase
{
public:

  /// Maximum Courant Friedrichs Lewy number in the grid for each phase
  real64 m_maxPhaseCFL;

  /// Maximum Courant-Friedrichs-Lewy number in the grid for each component
  real64 m_maxCompCFL;

  /**
   * @brief Construct a new CFLStatistics object
   * @param name instance name in data-repository
   * @param parent the instance parent in data-repository
   */
  CFLStatistics( const string & name, dataRepository::Group * const parent );
};

/**
 * @brief Reponsible of computing physical statistics over the grid, registering the result in the
 *        data repository, but not storing / outputing it by itself. It does not have mutable state
 *        except the encountered issues.
 */
class StatsAggregator : public StatsAggregatorBase< StatsAggregator >
{
public:

  using Base = StatsAggregatorBase< StatsAggregator >;

  /**
   * @brief the associated view keys
   */
  struct ViewKeys
  {
    /// String for the cfl statistics group
    constexpr static char const * cflStatisticsString() { return "cflStatistics"; }
  };

  /**
   * @brief Construct a new Stats Aggregator object
   * @param ownerName the unique name of the entity requesting the statistics.
   *                  An error is thrown if not unique in this context.
   */
  StatsAggregator( dataRepository::DataContext const & ownerDataContext );

  /**
   * @brief Enable the computation of any statistics, initialize data structure to collect them.
   *        Register the resulting data wrappers so they will be targeted by TimeHistory output
   * @param solver flow solver object to retrieve:
                   - the simulated regions,
                   - fields for statistics computation.
   * @param meshBodies The Group containing the MeshBody objects
   */
  void initStatisticsAggregation( dataRepository::Group & meshBodies,
                                  CompositionalMultiphaseBase & solver );

  /**
   * @brief Enable the computation of region statistics, initialize data structure to collect them.
   *        Register the resulting data wrappers so they will be targeted by TimeHistory output
   * @note Must be called in or after the "registerDataOnMesh" initialization phase
   * @param meshBodies The Group containing the MeshBody objects
   */
  void enableRegionStatisticsAggregation();

  /**
   * @brief Register the results structs & wrappers so they will be targeted by TimeHistory output
   * @note Must be called in or after the "registerDataOnMesh" initialization phase
   * @param meshBodies The Group containing the MeshBody objects
   */
  void enableCFLStatistics();

  /**
   * @brief set the statistics as dirty, ensuring isComputed() will be false until the next computation.
   */
  void setDirty();

  /**
   * @brief Compute CFL numbers
   * @param[in] time current time
   * @param[in] dt the time step size
   * @param[in] domain the domain partition
   * @return false if there was a problem that prevented the statistics to be computed correctly.
   */
  bool computeCFLNumbers( real64 const time,
                          real64 const dt,
                          DomainPartition & domain );

  integer getNumPhases() const
  { return m_numPhases; }

  integer getNumComponents() const
  { return m_numComponents; }

  CFLStatistics * getCFLStatistics( DomainPartition & domain ) const;

  CFLStatistics & getCflStatistics( MeshLevel & mesh ) const;

  // template implementations
  /// @cond DO_NOT_DOCUMENT

  void initStats( RegionStatistics & stats, real64 time ) const;
  void computeSubRegionRankStats( CellElementSubRegion & subRegion, RegionStatistics & subRegionStats ) const;
  void aggregateStats( RegionStatistics & stats, RegionStatistics const & other ) const;
  void mpiAggregateStats( RegionStatistics & stats ) const;
  void postAggregateStats( RegionStatistics & stats );

  /// @endcond

private:

  SolverType * m_solver = nullptr;

  AggregatorParameters m_params;

  StatsState m_cflStatsState;

  integer m_numPhases;

  integer m_numComponents;

};

} /* namespace compositionalMultiphaseStatistics */

} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASESTATISTICSAGGREGATOR_HPP_ */
