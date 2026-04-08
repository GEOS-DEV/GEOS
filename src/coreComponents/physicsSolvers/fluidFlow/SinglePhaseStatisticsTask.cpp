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
 * @file SinglePhaseStatisticsTask.cpp
 */

#include "SinglePhaseStatisticsTask.hpp"

#include "physicsSolvers/LogLevelsInfo.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"

namespace geos
{

using namespace dataRepository;

namespace singlePhaseStatistics
{

StatsTask::StatsTask( const string & name, Group * const parent ):
  Base( name, parent )
{
  registerWrapper( viewKeyStruct::setNamesString(), &m_setNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Optional targeted mesh element set(s) for which the statistics will be restricted. "
                    "A set can be be defined by a 'Geometry' component, or correspond to an imported index set "
                    "of mesh elements in case of an external mesh. "
                    "If empty, all mesh regions will be processed. "
                    "Be aware that only the regions that are processed by the solver will be taken into account." );

  addLogLevel< logInfo::Statistics >();
}

void StatsTask::postInputInitialization()
{
  Base::postInputInitialization();

  GEOS_THROW_IF_EQ_MSG( m_solver, nullptr,
                        "To identify simulated regions, a solver must be provided.",
                        InputError, getWrapperDataContext( getSolverWrapperKey() ) );

  if( !dynamicCast< SinglePhaseBase * >( m_solver ) )
  {
    GEOS_THROW( "Incompatible solver selected, a single-phase solver is expected",
                InputError, getDataContext() );
  }

  m_aggregator = std::make_unique< StatsAggregator >( getDataContext(), true );
}

void StatsTask::registerDataOnMesh( Group & meshBodies )
{
  // for now, this guard is needed to avoid breaking the xml schema generation
  if( m_solver == nullptr || m_aggregator == nullptr )
    return;

  if( m_writeCSV || isLogLevelActive< logInfo::Statistics >( this->getLogLevel()) )
  {
    // expected to work as check is done in postInputInitialization()
    SinglePhaseBase * castedSolver = dynamicCast< SinglePhaseBase * >( m_solver );
    GEOS_ERROR_IF_EQ_MSG( castedSolver, nullptr,
                          GEOS_FMT( "{} {}: Unexpected error (solver pointer changed?)", catalogName(), getDataContext() ) );
    m_aggregator->initStatisticsAggregation( meshBodies, *castedSolver );
  }
  else
  {
    GEOS_WARNING( GEOS_FMT( "{} {}: No computing option enabled, no output is scheduled.",
                            catalogName(), getDataContext() ) );
  }

  m_aggregator->enableRegionStatisticsAggregation( m_setNames );

  m_aggregator->forRegionStatistics( [&] ( MeshLevel & mesh, RegionStatistics & )
  {
    prepareLogTableLayouts( mesh.getName() );
    prepareCsvTableLayouts( mesh.getName() );
  } );
}

void StatsTask::prepareLogTableLayouts( string_view meshName )
{
  // only output from rank 0
  if( MpiWrapper::commRank() != 0 )
    return;

  TableLayout const tableLayout = TableLayout()
                                    .setTitle( GEOS_FMT( "Statistics: {} / Discretization: {}",
                                                         getName(), meshName ) );

  m_logFormatters.emplace( meshName, std::make_unique< TableTextFormatter >( tableLayout ) );
}

void StatsTask::prepareCsvTableLayouts( string_view meshName )
{
  // only output from rank 0
  if( MpiWrapper::commRank() != 0 || !m_writeCSV )
    return;

  string_view massUnit = units::getSymbol( m_solver->getMassUnit() );

  TableLayout tableLayout( {
        TableLayout::Column( GEOS_FMT( "Time [{}]", units::getSymbol( units::Unit::Time ))),
        TableLayout::Column( "Set" ),
        TableLayout::Column( "Region" ),
        TableLayout::Column( GEOS_FMT( "Min pressure [{}]", units::getSymbol( units::Unit::Pressure ))),
        TableLayout::Column( GEOS_FMT( "Average pressure [{}]", units::getSymbol( units::Unit::Pressure )) ),
        TableLayout::Column( GEOS_FMT( "Max pressure [{}]", units::getSymbol( units::Unit::Pressure ) ) ),
        TableLayout::Column( GEOS_FMT( "Min delta pressure [{}]", units::getSymbol( units::Unit::Pressure ))),
        TableLayout::Column( GEOS_FMT( "Max delta pressure [{}]", units::getSymbol( units::Unit::Pressure ))),
        TableLayout::Column( GEOS_FMT( "Min temperature [{}]", units::getSymbol( units::Unit::Temperature ) )),
        TableLayout::Column( GEOS_FMT( "Average temperature [{}]", units::getSymbol( units::Unit::Temperature ) )),
        TableLayout::Column( GEOS_FMT( "Max temperature [{}]", units::getSymbol( units::Unit::Temperature ) )),
        TableLayout::Column( GEOS_FMT( "Total dynamic pore volume [{}]", units::getSymbol( units::Unit::ReservoirVolume ) )),
        TableLayout::Column( GEOS_FMT( "Total fluid mass [{}]", massUnit )),
      } );

  auto & csvFormatter = m_csvFormatters.get_inserted( string( meshName ) );
  csvFormatter = std::make_unique< TableCSVFormatter >( tableLayout );

  // output CSV header
  std::ofstream outputFile( getCsvFileName( meshName ) );
  outputFile << csvFormatter->headerToString();
}

string StatsTask::getCsvFileName( string_view meshName ) const
{ return GEOS_FMT( "{}/{}.csv", m_outputDir, meshName ); }

bool StatsTask::execute( real64 const time_n,
                         real64 const dt,
                         integer const GEOS_UNUSED_PARAM( cycleNumber ),
                         integer const GEOS_UNUSED_PARAM( eventCounter ),
                         real64 const GEOS_UNUSED_PARAM( eventProgress ),
                         DomainPartition & )
{
  // current statistics time is after solver resolution: time_n (timestep start) + dt
  real64 statsTime = time_n + dt;

  GEOS_ERROR_IF( !m_aggregator,
                 "No statistics aggregator initialized!", getDataContext() );

  m_aggregator->computeRegionsStatistics( statsTime );

  m_aggregator->forRegionStatistics( [&] ( MeshLevel & mesh, RegionStatistics & meshLevelStats )
  {
    outputLogStats( statsTime, mesh, meshLevelStats );
    outputCsvStats( statsTime, mesh, meshLevelStats );
  } );

  return false;
}

void StatsTask::outputLogStats( real64 const statsTime,
                                MeshLevel & mesh,
                                RegionStatistics & meshLevelStats )
{
  if( MpiWrapper::commRank() > 0 || !isLogLevelActive< logInfo::Statistics >( this->getLogLevel() ) )
    return;

  auto const formatterIter = m_logFormatters.find( mesh.getName() );
  if( formatterIter==m_logFormatters.end())
    return;

  TableTextFormatter const & formatter = *formatterIter->second;
  TableData tableData;
  static constexpr auto merge = CellType::MergeNext;

  string_view massUnit = units::getSymbol( m_solver->getMassUnit() );
  string_view pressureUnit = units::getSymbol( units::Pressure );
  string_view tempUnit = units::getSymbol( units::Temperature );
  string_view resVolUnit = units::getSymbol( units::ReservoirVolume );

  tableData.getErrorsList().appendErrors( m_aggregator->getWarnings() );

  tableData.addRow( "Statistics time", merge, merge, statsTime );

  // lamda to apply for each region statistics
  auto const outputRegionStats = [&] ( string_view title,
                                       RegionStatistics & stats )
  {
    tableData.addSeparator();
    tableData.addRow( merge, merge, merge, "" );
    tableData.addRow( merge, merge, merge, title );
    tableData.addSeparator();

    tableData.addRow( "statistics", "min", "average", "max" );
    tableData.addSeparator();

    tableData.addRow( GEOS_FMT( "Pressure [{}]", pressureUnit ),
                      stats.m_minPressure, stats.m_averagePressure, stats.m_maxPressure );
    tableData.addRow( GEOS_FMT( "Delta pressure [{}]", pressureUnit ),
                      stats.m_minDeltaPressure, "/", stats.m_maxDeltaPressure );
    tableData.addRow( GEOS_FMT( "Temperature [{}]", tempUnit ),
                      stats.m_minTemperature, stats.m_averageTemperature, stats.m_maxTemperature );

    tableData.addSeparator();

    tableData.addRow( GEOS_FMT( "Total dynamic pore volume [{}]", resVolUnit ),
                      "all", CellType::MergeNext, stats.m_totalDynamicPoreVolume );

    tableData.addRow( GEOS_FMT( "Total fluid mass [{}]", massUnit ),
                      "all", CellType::MergeNext, stats.m_totalMass );
  };

  // apply the output lambda for the mesh sets and regions
  string const meshTitle = GEOS_FMT( "Discretization: '{}'", mesh.getName() );
  bool const logPerRegion = isLogLevelActive< logInfo::StatisticsPerRegion >( this->getLogLevel() );

  if( m_aggregator->isTargetingMultipleSets() )
  {
    string const allSetsTitle = m_aggregator->isRestrictedToSets() ? "Selected sets" : "All elements";
    outputRegionStats( GEOS_FMT( "{} / {} ({} elements)",
                                 meshTitle, allSetsTitle, meshLevelStats.m_elemCount ),
                       meshLevelStats );
  }

  m_aggregator->forRegionStatistics( mesh, meshLevelStats, false,
                                     [&] ( StatsAggregator::MeshLevelSet meshSet,
                                           RegionStatistics & meshSetStats )
  {
    string const setTitle = m_aggregator->isRestrictedToSets() ?
                            GEOS_FMT( "Element set: '{}'", meshSet.setName ) :
                            "All elements";

    outputRegionStats( GEOS_FMT( "{} / {} ({} elements)",
                                 meshTitle, setTitle, meshSetStats.m_elemCount ),
                       meshSetStats );

    if( logPerRegion )
    {
      m_aggregator->forRegionStatistics( meshSet, meshSetStats,
                                         [&] ( CellElementRegion & region,
                                               RegionStatistics & setRegionStats )
      {
        outputRegionStats( GEOS_FMT( "{} / {} / Region: '{}' ({} elements)",
                                     meshTitle, setTitle, region.getName(), setRegionStats.m_elemCount ),
                           setRegionStats );
      } );
    }
  } );

  // output to log
  GEOS_LOG_RANK_0( formatter.toString( tableData ) );
}

void StatsTask::outputCsvStats( real64 statsTime,
                                MeshLevel & mesh,
                                RegionStatistics & meshLevelStats )
{
  if( MpiWrapper::commRank() > 0 || m_writeCSV == 0 )
    return;

  auto const formatterIter = m_csvFormatters.find( mesh.getName() );
  if( formatterIter==m_csvFormatters.end())
    return;

  TableCSVFormatter const & formatter = *formatterIter->second;
  TableData tableData;

  stdVector< string > row;
  row.reserve( formatter.getLayout().getTotalLowermostColumnCount() );

  // lamda to apply for each region statistics
  auto const outputRegionStats = [&] ( string_view targetName,
                                       string_view setName,
                                       RegionStatistics & stats )
  {
    row.clear();
    row.insert( row.begin(),
                { std::to_string( statsTime ),
                  string( setName ),
                  string( targetName ),
                  std::to_string( stats.m_minPressure ),
                  std::to_string( stats.m_averagePressure ),
                  std::to_string( stats.m_maxPressure ),
                  std::to_string( stats.m_minDeltaPressure ),
                  std::to_string( stats.m_maxDeltaPressure ),
                  std::to_string( stats.m_minTemperature ),
                  std::to_string( stats.m_averageTemperature ),
                  std::to_string( stats.m_maxTemperature ),
                  std::to_string( stats.m_totalDynamicPoreVolume ),
                  std::to_string( stats.m_totalMass ),
                } );

    tableData.addRow( row );
  };

  // apply the output lambda for the mesh sets and regions
  string_view allSetsTitle = "Selected sets";

  outputRegionStats( mesh.getName(), allSetsTitle, meshLevelStats );

  m_aggregator->forRegionStatistics( mesh, meshLevelStats, false,
                                     [&] ( StatsAggregator::MeshLevelSet meshSet,
                                           RegionStatistics & meshSetStats )
  {
    outputRegionStats( meshSet.mesh.getName(), meshSet.setName, meshSetStats );

    m_aggregator->forRegionStatistics( meshSet, meshSetStats,
                                       [&] ( CellElementRegion & region,
                                             RegionStatistics & setRegionStats )
    {
      outputRegionStats( region.getName(), setRegionStats.getSetName(), setRegionStats );
    } );
  } );

  // append to csv file
  std::ofstream outputFile( getCsvFileName( mesh.getName() ), std::ios_base::app );
  outputFile << formatter.dataToString( tableData );
  outputFile.close();
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        StatsTask,
                        string const &, dataRepository::Group * const )

} /* namespace singlePhaseStatistics */

} /* namespace geos */
