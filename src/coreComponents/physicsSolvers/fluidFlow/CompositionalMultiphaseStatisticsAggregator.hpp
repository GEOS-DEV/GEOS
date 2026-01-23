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
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASESTATISTICSAGGREGATOR_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASESTATISTICSAGGREGATOR_HPP_

#include "common/DataTypes.hpp"
#include "common/StdContainerWrappers.hpp"
#include "dataRepository/DataContext.hpp"
#include "dataRepository/Group.hpp"
#include "mesh/CellElementRegion.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "mesh/MeshLevel.hpp"

namespace geos
{

class CompositionalMultiphaseBase;

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
 *        Attributes are public since the class is a POD. Can it be replaced by a wrapped-struct?
 */
class RegionStatistics : public dataRepository::Group
{
public:

  /// Time of statistics computation
  real64 m_time;

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
  /// total region uncompacted pore volume
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

  /**
   * @return the name of the data-repository object that is targeted by the statistics
   *         (mesh level / region / sub-region).
   */
  string_view getTargetName() const
  { return getName(); }

};

/**
 * @brief Output data group to contain the result of a given stat aggregator on the dataRepository.
 *        Attributes are public since the class is a POD. Can it be replaced by a wrapped-struct?
 */
class CFLStatistics : public dataRepository::Group
{
public:
  /// Time of statistics computation
  real64 m_time;

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
 * @todo repair 1D HDF5 outputs to enable stats HDF5 outputs
 */
class StatsAggregator
{
public:

  /**
   * @brief the associated view keys
   */
  struct ViewKeys
  {
    /// String for the discretization statistics group
    constexpr static char const * statisticsString() { return "statistics"; }
    /// String for the region statistics group
    constexpr static char const * regionsStatisticsString() { return "regionsStatistics"; }
    /// String for the cfl statistics group
    constexpr static char const * cflStatisticsString() { return "cflStatistics"; }
  };

  /**
   * @brief Standard function signature for any functor that applies on RegionStatistics instances
   *        param 0: TODO document
   *        param 1: TODO document
   * @tparam OWNER_T
   */
  template< typename OWNER_T >
  using RegionStatisticsFunctor = std::function< void ( OWNER_T &, RegionStatistics & ) >;

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

  void forRegionStatistics( RegionStatisticsFunctor< MeshLevel > const & functor ) const;

  void forRegionStatistics( MeshLevel & mesh,
                            RegionStatistics & meshRegionsStatistics,
                            RegionStatisticsFunctor< CellElementRegion > const & functor ) const;

  void forRegionStatistics( CellElementRegion & region,
                            RegionStatistics & regionStatistics,
                            RegionStatisticsFunctor< CellElementSubRegion > const & functor ) const;

  /**
   * @param[in] timeRequest The time for which we want to know if the statistics are computed.
   * @param[in] stats the statistics data structure we want to know if it has been computed
   * @return true if the statistics have been computed.
   */
  bool isComputed( real64 const timeRequest, RegionStatistics const & stats );

  /**
   * @brief Compute statistics on the mesh discretizations (average field pressure, etc)
   *        Results are reduced on rank 0, and broadcasted over all ranks.
   * @param[in] timeRequest The time for which we want to compute the statistics.
   * @return false if there was a problem that prevented the statistics to be computed correctly.
   */
  bool computeRegionsStatistics( real64 const timeRequest );

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

  CFLStatistics * getCFLStatistics( DomainPartition & domain ) const
  {return domain.getGroupPointer< CFLStatistics >( ViewKeys::cflStatisticsString() ); }

  CompositionalMultiphaseBase const * getSolver() const
  { return m_solver; }

  /**
   * @return The encountered issues during the last computing method call.
   */
  stdVector< string > const & getWarnings() const
  { return m_warnings; }

  /**
   * @return the name of the entity that needs the statistics.
   */
  string const & getOwnerName() const
  { return m_ownerDataContext.getTargetName(); }

  integer getNumPhases() const
  { return m_numPhases; }

  integer getNumComponents() const
  { return m_numComponents; }

  dataRepository::Group & getInstanceStatisticsGroup( MeshLevel & mesh ) const;

  RegionStatistics & getMeshRegionsStatistics( MeshLevel & mesh ) const;

  /**
   * @brief TODO
   * @throw InputError if no statistics data is found for the given region name.
   * @param mesh TODO
   * @param regionNname TODO
   * @return TODO
   */
  RegionStatistics & getRegionStatistics( MeshLevel & mesh, string_view regionName ) const;

  CFLStatistics & getCflStatisticsGroup( MeshLevel & mesh ) const;

private:

  struct StatsState {
    bool m_isEnabled = false;
    bool m_isDirty = false;
  };

  /// @see getOwnerName()
  dataRepository::DataContext const & m_ownerDataContext;

  CompositionalMultiphaseBase * m_solver = nullptr;

  dataRepository::Group * m_meshBodies = nullptr;

  AggregatorParameters m_params;

  stdVector< string > m_warnings;

  StatsState m_regionStatsState;

  StatsState m_cflStatsState;

  integer m_numPhases;

  integer m_numComponents;

  /**
   * @brief Initialize all statistics values to aggregable default values,
   *        before any computation / reduction for the current timestep.
   * @param stats the statistics instance
   * @param time start time of the current timestep (s)
   */
  void initStats( RegionStatistics & stats, real64 time ) const;

  void computeSubRegionRankStats( CellElementSubRegion & subRegion, RegionStatistics & subRegionStats ) const;

  /**
   * @brief Aggregate all instance statistics with those of another instance on the current rank.
   * @param stats the statistics instance
   * @param other the other instance to aggregate with.
   */
  void aggregateStats( RegionStatistics & stats, RegionStatistics const & other ) const;

  /**
   * @brief Aggregate all instance statistics with those of other ranks.
   * @param stats the statistics instance
   */
  void mpiAggregateStats( RegionStatistics & stats ) const;

  /**
   * @brief Do the final computations for the statistics. Must be called after computations & aggregations.
   * @param stats the statistics instance
   */
  void postAggregateStats( RegionStatistics & stats );

};

} /* namespace compositionalMultiphaseStatistics */

} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASESTATISTICSAGGREGATOR_HPP_ */
