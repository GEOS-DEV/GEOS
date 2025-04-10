
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
 * @file MemoryStatsOutput.cpp
 */

#include "MemoryStatsOutput.hpp"

#include "fileIO/LogLevelsInfo.hpp"
#include "common/MemoryInfos.hpp"

namespace geos
{

using namespace dataRepository;

namespace logInfo
{
struct MemoryStatsOutputTimer : public OutputTimerBase
{
  std::string_view getDescription() const override { return "MemoryStats output timing"; }
};
}

logInfo::OutputTimerBase const & MemoryStatsOutput::getTimerCategory() const
{
  static logInfo::MemoryStatsOutputTimer timer;
  return timer;
}

MemoryStatsOutput::MemoryStatsOutput( string const & name,
                                      Group * const parent ):
  OutputBase( name, parent ),
  m_writeCSV( 0 )
{
  bool const umpireStatsDefault = MemoryLogging::getInstance().isUmpireStatsLogOutputEnabled();
  addLogLevel< logInfo::UmpireStatistics >();
  getWrapper< integer >( string( Group::viewKeyStruct::logLevelString() ) ).
    setApplyDefaultValue( umpireStatsDefault ? 1 : 0 );

  bool const csvOutputDefault = MemoryLogging::getInstance().isUmpireStatsLogOutputEnabled();
  this->registerWrapper( viewKeysStruct::writeCSV, &m_writeCSV ).
    setApplyDefaultValue( csvOutputDefault ? 1 : 0 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "When set to 1, write the same statistics as the 'logLevel' allows to output in a CSV file" );
}

bool MemoryStatsOutput::execute( real64,
                                 real64,
                                 integer,
                                 integer,
                                 real64,
                                 DomainPartition & )
{
  auto & memLogging = MemoryLogging::getInstance();
  memLogging.enableUmpireStatsLogReport( isLogLevelActive< logInfo::UmpireStatistics >( getLogLevel() ) );
  memLogging.enableUmpireStatsCsvReport( m_writeCSV );
  memLogging.setUmpireStatsCsvReportFilename( GEOS_FMT( "{}_umpireStats.csv", getName() ) );

  memLogging.memoryStatsReport();

  return false;
}

REGISTER_CATALOG_ENTRY( OutputBase, MemoryStatsOutput, string const &, Group * const )

} /* namespace geos */
