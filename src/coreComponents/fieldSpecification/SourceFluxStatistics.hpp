/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhaseStatistics.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASESTATISTICS_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASESTATISTICS_HPP_

#include "physicsSolvers/FieldStatisticsBase.hpp"

class SourceFluxBoundaryCondition;

namespace geos
{

/**
 * @class SourceFluxStatistics
 *
 * Task class allowing for the computation of aggregate statistics of SourceFluxBoundaryCondition
 */
class SourceFluxStatistics : public FieldStatisticsBase< FlowSolverBase >
{
public:

  /**
   * @brief Statistics flux data in the whole mesh or in a precise region
   */
  struct Stats
  {
    /// fluid mass produced by the flux(es) (kg). Negative if injecting.
    real64 producedMass;
    /// flux(es) production rate (kg/s). Negative if injecting.
    real64 productionRate;
  };

  /**
   * @brief Constructor for the statistics class
   * @param[in] name the name of the task coming from the xml
   * @param[in] parent the parent group of the task
   */
  SourceFluxStatistics( const string & name,
                        Group * const parent );

  /// Accessor for the catalog name
  static string catalogName() { return "SourceFluxStatistics"; }

  /**
   * @defgroup Tasks Interface Functions
   *
   * This function implements the interface defined by the abstract TaskBase class
   */
  /**@{*/

  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  /**@}*/

  /**
   * @brief View keys
   */
  struct viewKeyStruct
  {
    /// @return The key for setName
    constexpr static char const * setNamesString() { return "setNames"; }
  };

protected:

  void postProcessInput() override;

private:
  using Base = FieldStatisticsBase< SinglePhaseBase >;

  /// the names of the SourceFlux(s) for which we want the statistics
  string_array m_fluxNames;
  /// Internal array of the SourceFlux(s) for which we want the statistics
  std::vector< SourceFluxBoundaryCondition * > m_fluxes;

  /**
   * @return a string used to name the wrapper that is added to each region that is simulated by
   * the solver. The string is based on the name of this Group so it is unique within the region.
   */
  string getRegionStatsName( string_view fluxName )
  {
    return GEOS_FMT( "{}_region_{}_Stats",
                     getName(), fluxName );
  }

  /**
   * @copydoc Group::registerDataOnMesh(Group &)
   */
  void registerDataOnMesh( Group & meshBodies ) override;

  /**
   * @param Stats the statistics that must be output in the log
   */
  static void logStats( Stats const & stats ) override;

};


/******************************** SourceFluxStatisticsKernel ********************************/
struct SourceFluxStatisticsKernel
{};


} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASESTATISTICS_HPP_ */
