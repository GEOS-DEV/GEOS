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

#include "common/StdContainerWrappers.hpp"
#include "dataRepository/Group.hpp"
#include "mesh/DomainPartition.hpp"
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

  /// average region pressure
  real64 m_averagePressure;
  /// minimum region pressure
  real64 m_minPressure;
  /// maximum region pressure
  real64 m_maxPressure;

  /// minimum region delta pressure
  real64 m_minDeltaPressure;
  /// maximum region delta pressure
  real64 m_maxDeltaPressure;

  /// average region temperature
  real64 m_averageTemperature;
  /// minimum region temperature
  real64 m_minTemperature;
  /// maximum region temperature
  real64 m_maxTemperature;

  /// total region pore volume
  real64 m_totalPoreVolume;
  /// total region uncompacted pore volume
  real64 m_totalUncompactedPoreVolume;
  /// phase region phase pore volume
  array1d< real64 > m_phasePoreVolume;

  /// region phase mass (trapped and non-trapped, immobile and mobile)
  array1d< real64 > m_phaseMass;
  /// trapped region phase mass
  array1d< real64 > m_trappedPhaseMass;
  /// non-trapped region phase mass
  array1d< real64 > m_nonTrappedPhaseMass;
  /// immobile region phase mass
  array1d< real64 > m_immobilePhaseMass;
  /// mobile region phase mass
  array1d< real64 > m_mobilePhaseMass;
  /// region component mass
  array2d< real64 > m_componentMass;

  // TODO: -> split to struct PressureStats...MassStats:
  // - VKS for struct name ("pressureStats"..."massStats")
  // - current RegionStatistics struct bits

  /**
   * @brief Construct a new Region Statistics object
   * @param name instance name in data-repository
   * @param parent the instance parent in data-repository
   */
  RegionStatistics( const string & name, dataRepository::Group * const parent );
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
 */
class StatsAggregator
{
public:

  /**
   * @brief the associated view keys
   */
  struct VKStruct
  {
    /// String for the region statistics group
    constexpr static char const * regionStatisticsString() { return "regionStatistics"; }
    /// String for the cfl statistics group
    constexpr static char const * cflStatisticsString() { return "cflStatistics"; }
  };

  using RegionFunctor = std::function< void (string_view, RegionStatistics const &) >;

  StatsAggregator();

  /**
   * @brief Set the reference flow solver object to retrieve:
            - the simulated regions,
            - fields for statistics computation.
   * @param solver The reference flow solver
   */
  void setFlowSolver( CompositionalMultiphaseBase & solver )
  { m_solver = &solver; }

  // void forDiscretizations( DomainPartition const &,
  //                          std::function< void(MeshLevel const &,
  //                                              string_array const & regionNames) > functor ) const;

  void forRegionStatistics( MeshLevel const & mesh,
                            RegionFunctor functor ) const;

  /**
   * @brief Enable the computation of region statistics, initialize data structure to collect them.
   *        Register the resulting data wrappers so they will be targeted by TimeHistory output
   * @note Must be called in "registerDataOnMesh" initialization phase
   * @param meshBodies The Group containing the MeshBody objects
   */
  void enableRegionStatistics( dataRepository::Group & meshBodies );

  /**
   * @brief Register the results structs & wrappers so they will be targeted by TimeHistory output
   * @note Must be called in "registerDataOnMesh" initialization phase
   * @param meshBodies The Group containing the MeshBody objects
   */
  void enableCFLStatistics( dataRepository::Group & meshBodies );

  /**
   * @brief Compute some statistics on a given mesh discretization (average field pressure, etc)
   *        Results are reduced on rank 0, and broadcasted over all ranks.
   * @param[in] time current time
   * @param[in] mesh the mesh level object
   * @param[in] regionNames the array of target region names
   * @return false if there was a problem that prevented the statistics to be computed correctly.
   */
  bool computeRegionsStatistics( real64 const time,
                                 MeshLevel & mesh,
                                 string_array const & regionNames );

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
  {return domain.getGroupPointer< CFLStatistics >( VKStruct::cflStatisticsString() ); }

  CompositionalMultiphaseBase const * getSolver() const
  { return m_solver; }

  /**
   * @return The encountered issues during the last computing method call.
   */
  stdVector< string > const & getIssues() const
  { return m_issues; }

private:

  CompositionalMultiphaseBase * m_solver;

  AggregatorParameters m_params;

  stdVector< string > m_issues;

  bool m_isRegionStatsEnabled;

  bool m_isCFLNumberEnabled;

};

} /* namespace compositionalMultiphaseStatistics */

} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASESTATISTICSAGGREGATOR_HPP_ */
