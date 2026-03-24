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
 * @file CompositionalMultiphaseStatistics.cpp
 */

#include "CompositionalMultiphaseStatisticsTask.hpp"

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "physicsSolvers/fluidFlow/LogLevelsInfo.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseHybridFVM.hpp"


namespace geos
{

using namespace dataRepository;

namespace compositionalMultiphaseStatistics
{

StatsTask::StatsTask( string const & name, Group * const parent ):
  Base( name, parent ),
  m_aggregator(),
  m_computeCFLNumbers( 0 ),
  m_computeRegionStatistics( 1 )
{
  registerWrapper( viewKeyStruct::setNamesString(), &m_setNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Optional targeted mesh element set(s) for which the statistics will be restricted. "
                    "A set can be be defined by a 'Geometry' component, or correspond to imported sets in case of an external mesh. "
                    "If empty, all mesh regions will be processed. "
                    "Be aware that only the regions that are computed by the solver will be taken into account." );

  registerWrapper( viewKeyStruct::computeCFLNumbersString(), &m_computeCFLNumbers ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to decide whether CFL numbers are computed or not" );

  registerWrapper( viewKeyStruct::computeRegionStatisticsString(), &m_computeRegionStatistics ).
    setApplyDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to decide whether region statistics are computed or not" );

  registerWrapper( viewKeyStruct::relpermThresholdString(), &m_relpermThreshold ).
    setApplyDefaultValue( 1e-6 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to decide whether a phase is considered mobile (when the relperm is above the threshold) or immobile (when the relperm is below the threshold) in metric 2" );

  addLogLevel< logInfo::CFL >();
  addLogLevel< logInfo::Statistics >();
}

void StatsTask::postInputInitialization()
{
  Base::postInputInitialization();

  GEOS_THROW_IF_EQ_MSG( m_solver, nullptr,
                        "To identify simulated regions, a solver must be provided.",
                        InputError, getWrapperDataContext( getSolverWrapperKey() ) );

  if( !dynamicCast< CompositionalMultiphaseBase * >( m_solver ) )
  {
    GEOS_THROW( "Incompatible solver selected, a compositional multiphase solver is expected",
                InputError, getDataContext() );
  }
  else if( dynamicCast< CompositionalMultiphaseHybridFVM * >( m_solver ) && m_computeCFLNumbers != 0 )
  {
    GEOS_THROW( "The option to compute CFL numbers is incompatible with CompositionalMultiphaseHybridFVM",
                InputError, getDataContext() );
  }

  m_aggregator = std::make_unique< StatsAggregator >( getDataContext(), true );
}

void StatsTask::registerDataOnMesh( Group & meshBodies )
{
  // for now, this guard is needed to avoid breaking the xml schema generation
  if( m_solver == nullptr || m_aggregator == nullptr )
    return;

  prepareFluidMetaData();

  if( m_computeRegionStatistics || m_computeCFLNumbers )
  {
    // expected to work as check is done in postInputInitialization()
    CompositionalMultiphaseBase * castedSolver = dynamicCast< CompositionalMultiphaseBase * >( m_solver );
    GEOS_ERROR_IF_EQ_MSG( castedSolver, nullptr,
                          GEOS_FMT( "{} {}: Unexpected error (solver pointer changed?)", catalogName(), getDataContext() ) );
    m_aggregator->initStatisticsAggregation( meshBodies, *castedSolver );
  }
  else
  {
    GEOS_WARNING( GEOS_FMT( "{} {}: No computing option enabled, no output is scheduled.",
                            catalogName(), getDataContext() ) );
  }

  if( m_computeRegionStatistics )
    m_aggregator->enableRegionStatisticsAggregation( m_setNames );

  // if we have to compute CFL numbers later, we need to register additional variables
  if( m_computeCFLNumbers )
    m_aggregator->enableCFLStatistics();

  m_aggregator->forRegionStatistics( [&] ( MeshLevel & mesh, RegionStatistics & )
  {
    prepareLogTableLayouts( mesh.getName() );
    prepareCsvTableLayouts( mesh.getName() );
  } );
}

void StatsTask::prepareFluidMetaData()
{
  using namespace constitutive;

  ConstitutiveManager const & constitutiveManager = this->getGroupByPath< ConstitutiveManager >( "/Problem/domain/Constitutive" );
  MultiFluidBase const & fluid = constitutiveManager.getGroup< MultiFluidBase >( m_solver->referenceFluidModelName() );

  m_fluid.m_numPhases = fluid.numFluidPhases();
  m_fluid.m_numComps = fluid.numFluidComponents();

  m_fluid.m_phaseNames = fluid.phaseNames();
  m_fluid.m_compNames = fluid.componentNames();

  m_fluid.m_phaseCompNames.resize( m_fluid.m_numPhases * m_fluid.m_numComps, string() );

  for( int ip = 0; ip < m_fluid.m_numPhases; ++ip )
    for( int ic = 0; ic < m_fluid.m_numComps; ++ic )
      m_fluid.phaseCompName( ip, ic ) = GEOS_FMT( "{} / {}", m_fluid.m_phaseNames[ip], m_fluid.m_compNames[ic] );
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

  integer const numPhases = m_solver->numFluidPhases();
  integer const numComps = m_solver->numFluidComponents();

  auto addPhaseColumns = []( TableLayout & ptableLayout, string const & description, string_view punit,
                             integer pnumPhases )
  {
    for( int ip = 0; ip < pnumPhases; ++ip )
      ptableLayout.addColumn( GEOS_FMT( "{} (phase {}) [{}]",
                                        description, ip, punit ) );
  };
  auto addPhaseCompColumns = []( TableLayout & ptableLayout, string const & description, string_view punit,
                                 integer pnumPhases, integer pnumComps )
  {
    for( int ip = 0; ip < pnumPhases; ++ip )
      for( int ic = 0; ic < pnumComps; ++ic )
        ptableLayout.addColumn( GEOS_FMT( "{} (component {} / phase {}) [{}]",
                                          description, ic, ip, punit ) );
  };

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
      } );
  addPhaseColumns( tableLayout, "Phase dynamic pore volume", units::getSymbol( units::Unit::ReservoirVolume ), numPhases );
  addPhaseColumns( tableLayout, "Phase mass", massUnit, numPhases );
  addPhaseColumns( tableLayout, "Trapped phase mass (metric 1)", massUnit, numPhases );
  addPhaseColumns( tableLayout, "Non-trapped phase mass (metric 1)", massUnit, numPhases );
  addPhaseColumns( tableLayout, "Immobile phase mass (metric 2)", massUnit, numPhases );
  addPhaseColumns( tableLayout, "Mobile phase mass (metric 2)", massUnit, numPhases );
  addPhaseCompColumns( tableLayout, "Component mass", massUnit, numPhases, numComps );

  auto & csvFormatter = m_csvFormatters.get_inserted( string( meshName ) );
  csvFormatter = std::make_unique< TableCSVFormatter >( tableLayout );

  // output CSV header
  std::ofstream outputFile( getCsvFileName( meshName ) );
  outputFile << csvFormatter->headerToString();
  GEOS_LOG( GEOS_FMT( "table {} : {}", meshName, csvFormatter->headerToString() ) );     // TODO : remove this log
}

string StatsTask::getCsvFileName( string_view meshName ) const
{ return GEOS_FMT( "{}/{}.csv", m_outputDir, meshName ); }

bool StatsTask::execute( real64 const time_n,
                         real64 const dt,
                         integer const GEOS_UNUSED_PARAM( cycleNumber ),
                         integer const GEOS_UNUSED_PARAM( eventCounter ),
                         real64 const GEOS_UNUSED_PARAM( eventProgress ),
                         DomainPartition & domain )
{
  // current time is time_n + dt. TODO: verify implication of events ordering in 'time_n+dt' validity
  real64 statsTime = time_n + dt;

  GEOS_ERROR_IF( !m_aggregator,
                 "No statistics aggregator initialized!", getDataContext() );

  m_aggregator->computeRegionsStatistics( statsTime );

  m_aggregator->forRegionStatistics( [&] ( MeshLevel & mesh, RegionStatistics & meshLevelStats )
  {
    if( m_computeRegionStatistics )
    {
      outputLogStats( statsTime, mesh, meshLevelStats );
      outputCsvStats( statsTime, mesh, meshLevelStats );
    }
  } );

  if( m_computeCFLNumbers )
    m_aggregator->computeCFLNumbers( statsTime, dt, domain );

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
                      "all",
                      CellType::MergeNext,
                      stats.m_totalPoreVolume );
    tableData.addRow( GEOS_FMT( "Phase dynamic pore volume [{}]", resVolUnit ),
                      stringutilities::joinLambda( m_fluid.m_phaseNames, "\n", []( auto data ) { return data[0]; } ),
                      CellType::MergeNext,
                      stringutilities::joinLambda( stats.m_phaseDynamicPoreVolume, "\n", []( auto data ) { return data[0]; } ) );

    tableData.addSeparator();

    tableData.addRow( GEOS_FMT( "Phase mass [{}]", massUnit ),
                      stringutilities::joinLambda( m_fluid.m_phaseNames, "\n", []( auto data ) { return data[0]; } ),
                      CellType::MergeNext,
                      stringutilities::joinLambda( stats.m_phaseMass, "\n", []( auto data ) { return data[0]; } ) );

    tableData.addSeparator();

    tableData.addRow( GEOS_FMT( "Trapped phase mass (metric 1) [{}]", massUnit ),
                      stringutilities::joinLambda( m_fluid.m_phaseNames, "\n", []( auto value ) { return value[0]; } ),
                      CellType::MergeNext,
                      stringutilities::joinLambda( stats.m_trappedPhaseMass, "\n", []( auto value ) { return value[0]; } ) );
    tableData.addRow( GEOS_FMT( "Non-trapped phase mass (metric 1) [{}]", massUnit ),
                      stringutilities::joinLambda( m_fluid.m_phaseNames, "\n", []( auto value ) { return value[0]; } ),
                      CellType::MergeNext,
                      stringutilities::joinLambda( stats.m_nonTrappedPhaseMass, "\n", []( auto value ) { return value[0]; } ) );

    tableData.addSeparator();

    tableData.addRow( GEOS_FMT( "Immobile phase mass (metric 2) [{}]", massUnit ),
                      stringutilities::joinLambda( m_fluid.m_phaseNames, "\n", []( auto value ) { return value[0]; } ),
                      CellType::MergeNext,
                      stringutilities::joinLambda( stats.m_immobilePhaseMass, "\n", []( auto value ) { return value[0]; } )  );
    tableData.addRow( GEOS_FMT( "Mobile phase mass (metric 2) [{}]", massUnit ),
                      stringutilities::joinLambda( m_fluid.m_phaseNames, "\n", []( auto value ) { return value[0]; } ),
                      CellType::MergeNext,
                      stringutilities::joinLambda( stats.m_mobilePhaseMass, "\n", []( auto value ) { return value[0]; } ) );

    tableData.addSeparator();

    tableData.addRow( GEOS_FMT( "Component mass [{}]", massUnit ),
                      stringutilities::join( m_fluid.m_phaseCompNames, '\n' ),
                      CellType::MergeNext,
                      stringutilities::join( stats.m_componentMass, '\n' ) );
  };

  // apply the output lambda for the mesh sets and regions
  string const meshTitle = GEOS_FMT( "Discretization: '{}'", mesh.getName() );
  bool const logPerRegion = isLogLevelActive< logInfo::StatisticsPerRegion >( this->getLogLevel() );

  if( m_aggregator->isTargetingMultipleSets() )
  {
    string const allSetsTitle = m_aggregator->isRestrictedToSets() ? "Selected sets" : "All elements";
    outputRegionStats( GEOS_FMT( "{} / {}", meshTitle, allSetsTitle ),
                       meshLevelStats );
  }

  m_aggregator->forRegionStatistics( mesh, meshLevelStats, false,
                                     [&] ( StatsAggregator::MeshLevelSet meshSet,
                                           RegionStatistics & meshSetStats )
  {
    string const setTitle = m_aggregator->isRestrictedToSets() ?
                            GEOS_FMT( "Element set: '{}'", meshSet.setName ) :
                            "All elements";

    outputRegionStats( GEOS_FMT( "{} / {}", meshTitle, setTitle ),
                       meshLevelStats );

    if( logPerRegion )
    {
      m_aggregator->forRegionStatistics( meshSet, meshSetStats,
                                         [&] ( CellElementRegion & region,
                                               RegionStatistics & setRegionStats )
      {
        outputRegionStats( GEOS_FMT( "{} / {} / Region: '{}'", meshTitle, setTitle, region.getName() ),
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

    auto addPhaseValues = []( auto & list, auto const & values )
    {
      for( auto value : values )
        list.emplace_back( std::to_string( value ) );
    };

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
                  std::to_string( stats.m_totalPoreVolume ),
                } );
    addPhaseValues( row, stats.m_phaseDynamicPoreVolume );
    addPhaseValues( row, stats.m_phaseMass );
    addPhaseValues( row, stats.m_trappedPhaseMass );
    addPhaseValues( row, stats.m_nonTrappedPhaseMass );
    addPhaseValues( row, stats.m_immobilePhaseMass );
    addPhaseValues( row, stats.m_mobilePhaseMass );
    addPhaseValues( row, stats.m_componentMass ); // TODO verify phase / comp ordering

    tableData.addRow( row );
  };

  // apply the lambda for each region and, finally, the mesh summary
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

}   /* namespace compositionalMultiphaseStatistics */

} /* namespace geos */
