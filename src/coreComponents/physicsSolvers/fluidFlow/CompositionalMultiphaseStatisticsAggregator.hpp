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
 * @file CompositionalMultiphaseStatistics.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASESTATISTICSAGGREGATOR_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASESTATISTICSAGGREGATOR_HPP_

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

struct RegionStatistics
{
  /// Time of statistics computation
  real64 time;

  /// average region pressure
  real64 averagePressure;
  /// minimum region pressure
  real64 minPressure;
  /// maximum region pressure
  real64 maxPressure;

  /// minimum region delta pressure
  real64 minDeltaPressure;
  /// maximum region delta pressure
  real64 maxDeltaPressure;

  /// average region temperature
  real64 averageTemperature;
  /// minimum region temperature
  real64 minTemperature;
  /// maximum region temperature
  real64 maxTemperature;

  /// total region pore volume
  real64 totalPoreVolume;
  /// total region uncompacted pore volume
  real64 totalUncompactedPoreVolume;
  /// phase region phase pore volume
  array1d< real64 > phasePoreVolume;

  /// region phase mass (trapped and non-trapped, immobile and mobile)
  array1d< real64 > phaseMass;
  /// trapped region phase mass
  array1d< real64 > trappedPhaseMass;
  /// immobile region phase mass
  array1d< real64 > immobilePhaseMass;
  /// region component mass
  array2d< real64 > componentMass;

  // TODO: -> split to struct PressureStats...MassStats:
  // - VKS for struct name ("pressureStats"..."massStats")
  // - current RegionStatistics struct bits
};

struct CFLStatistics
{
  /// Time of statistics computation
  real64 time;

  /// Maximum Courant Friedrichs Lewy number in the grid for each phase
  real64 maxPhaseCFL;

  /// Maximum Courant-Friedrichs-Lewy number in the grid for each component
  real64 maxCompCFL;
};

class StatsAggregator
{
public:

  struct viewKeyStruct
  {
    /// String for the region statistics
    constexpr static char const * regionStatisticsString() { return "regionStatistics"; }
  };

  using RegionFunctor = std::function< void (string, RegionStatistics const &) >;

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

  // void forRegionStatistics( DomainPartition const &,
  //                           MeshLevel const & mesh,
  //                           RegionFunctor functor ) const;

  /**
   * @brief Register the results structs & wrappers so they will be targeted by TimeHistory output
   * @note Must be called in "registerDataOnMesh" initialization phase
   * @param solver The flow solver
   * @param meshBodies The Group containing the MeshBody objects
   */
  void enableRegionStatistics( dataRepository::Group & meshBodies );

  /**
   * @brief Register the results structs & wrappers so they will be targeted by TimeHistory output
   * @note Must be called in "registerDataOnMesh" initialization phase
   * @param solver The flow solver
   * @param meshBodies The Group containing the MeshBody objects
   */
  void enableCFLStatistics( dataRepository::Group & meshBodies );

  /**
   * @brief Compute some statistics on a given mesh discretization (average field pressure, etc)
   * @param[in] time current time
   * @param[in] mesh the mesh level object
   * @param[in] regionNames the array of target region names
   */
  void computeDiscretizationStatistics( real64 const time,
                                        MeshLevel & mesh,
                                        string_array const & regionNames ) const;

  /**
   * @brief Compute CFL numbers
   * @param[in] time current time
   * @param[in] dt the time step size
   * @param[in] domain the domain partition
   */
  void computeCFLNumbers( real64 const time,
                          real64 const dt,
                          DomainPartition & domain ) const;

  CFLStatistics const & getCFLStatistics() const
  { return m_cflStats; }

private:

  CompositionalMultiphaseBase * m_solver;

  AggregatorParameters m_params;

  CFLStatistics m_cflStats;

  bool m_isRegionStatsEnabled;

  bool m_isCFLNumberEnabled;

};

} /* namespace compositionalMultiphaseStatistics */

} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASESTATISTICSAGGREGATOR_HPP_ */
