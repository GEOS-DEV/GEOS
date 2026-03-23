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
 * @file CompositionalMultiphaseStatisticsTask.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASESTATISTICSTASK_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASESTATISTICSTASK_HPP_

#include "common/DataTypes.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "physicsSolvers/FieldStatisticsBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseStatisticsAggregator.hpp"
#include <memory>

namespace geos
{

class CompositionalMultiphaseBase;

namespace compositionalMultiphaseStatistics
{

/**
 * @class compositionalMultiphaseStatistics::StatsStats
 *
 * Task class allowing for the computation of aggregate statistics in compositional multiphase simulations
 */
class StatsTask : public FieldStatisticsBase< CompositionalMultiphaseBase >
{
public:

  /**
   * @brief Constructor for the statistics class
   * @param[in] name the name of the task coming from the xml
   * @param[in] parent the parent group of the task
   */
  StatsTask( const string & name, dataRepository::Group * const parent );

  /// Accessor for the catalog name
  static string catalogName() { return "CompositionalMultiphaseStatistics"; }

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

  StatsAggregator & getStatisticsAggregator()
  { return *m_aggregator; }

  StatsAggregator const & getStatisticsAggregator() const
  { return *m_aggregator; }

private:

  using Base = FieldStatisticsBase< CompositionalMultiphaseBase >;

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct
  {
    /// String for optionnal targeted element set(s)
    constexpr static char const * setNamesString() { return "setNames"; }
    /// String for the flag deciding the computation of the CFL numbers
    constexpr static char const * computeCFLNumbersString() { return "computeCFLNumbers"; }
    /// String for the flag deciding the computation of the region statistics
    constexpr static char const * computeRegionStatisticsString() { return "computeRegionStatistics"; }
    /// String for the relperm threshold
    constexpr static char const * relpermThresholdString() { return "relpermThreshold"; }
  };

  void postInputInitialization() override;

  void registerDataOnMesh( Group & meshBodies ) override;

  void prepareFluidMetaData();

  void prepareLogTableLayouts( string_view tableName );

  void prepareCsvTableLayouts( string_view tableName );

  string getCsvFileName( string_view meshName ) const;

  void outputLogStats( real64 statsTime,
                       MeshLevel & mesh,
                       RegionStatistics & meshSetsStatistics );

  void outputCsvStats( real64 statsTime,
                       MeshLevel & mesh,
                       RegionStatistics & meshSetsStatistics );

  /// For each discretization (MeshLevel name), table formatter for log output.
  stdMap< string, std::unique_ptr< TableTextFormatter > > m_logFormatters;

  /// For each discretization (MeshLevel name), table formatter for csv output.
  stdMap< string, std::unique_ptr< TableCSVFormatter > > m_csvFormatters;

  // mesh statistics aggregator
  std::unique_ptr< StatsAggregator > m_aggregator;

  /// Optionnal targeted element set(s)
  string_array m_setNames;

  /// Flag to decide whether CFL numbers are computed or not
  integer m_computeCFLNumbers;

  /// Flag to decide whether region statistics are computed or not
  integer m_computeRegionStatistics;

  /// Threshold to decide whether a phase is considered "mobile" or not
  real64 m_relpermThreshold;

  struct FluidMetaData
  {
    integer m_numPhases;
    integer m_numComps;
    stdVector< string > m_phaseNames;
    stdVector< string > m_compNames;
    stdVector< string > m_phaseCompNames;

    /**
     * @param phaseId index of the phase for which we want the name
     * @param compId index of the component for which we want the name
     * @return string& reference of the name string
     */
    string & phaseCompName( integer phaseId, integer compId )
    { return m_phaseCompNames[phaseId * m_numComps + compId]; }

  } m_fluid;

};

} /* namespace compositionalMultiphaseStatistics */

} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASESTATISTICSTASK_HPP_ */
