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
 * @file SinglePhaseStatisticsTask.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASESTATISTICS_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASESTATISTICS_HPP_

#include "physicsSolvers/FieldStatisticsBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseStatisticsAggregator.hpp"

namespace geos
{

class SinglePhaseBase;

namespace singlePhaseStatistics
{

/**
 * Task class allowing for the computation of aggregate statistics in single phase simulations
 */
class StatsTask : public FieldStatisticsBase< SinglePhaseBase >
{
public:

  /**
   * @brief Constructor for the statistics class
   * @param[in] name the name of the task coming from the xml
   * @param[in] parent the parent group of the task
   */
  StatsTask( const string & name, dataRepository::Group * const parent );

  /// Accessor for the catalog name
  static string catalogName() { return "SinglePhaseStatistics"; }

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

  using Base = FieldStatisticsBase< SinglePhaseBase >;

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct
  {
    /// String for optionnal targeted element set(s)
    constexpr static char const * setNamesString() { return "setNames"; }
  };

  void postInputInitialization() override;

  void registerDataOnMesh( Group & meshBodies ) override;

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

};

} /* namespace singlePhaseStatistics */

} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASESTATISTICS_HPP_ */
