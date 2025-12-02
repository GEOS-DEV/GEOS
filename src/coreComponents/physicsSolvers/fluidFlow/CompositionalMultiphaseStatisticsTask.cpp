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

#include "common/DataTypes.hpp"
#include "common/StdContainerWrappers.hpp"
#include "common/format/Format.hpp"
#include "mesh/DomainPartition.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "physicsSolvers/fluidFlow/LogLevelsInfo.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseHybridFVM.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "common/format/table/TableData.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "common/format/table/TableLayout.hpp"
#include <memory>


namespace geos
{

using namespace constitutive;
using namespace fields;
using namespace dataRepository;

namespace compositionalMultiphaseStatistics
{

StatsTask::StatsTask( const string & name, Group * const parent ):
  Base( name, parent ),
  m_aggregator(),
  m_computeCFLNumbers( 0 ),
  m_computeRegionStatistics( 1 )
{
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

  m_aggregator.setFlowSolver( *m_solver );

  if( dynamicCast< CompositionalMultiphaseHybridFVM * >( m_solver ) && m_computeCFLNumbers != 0 )
  {
    GEOS_THROW( GEOS_FMT( "{} {}: the option to compute CFL numbers is incompatible with CompositionalMultiphaseHybridFVM",
                          catalogName(), getDataContext() ),
                InputError );
  }
}

void StatsTask::registerDataOnMesh( Group & meshBodies )
{
  // for now, this guard is needed to avoid breaking the xml schema generation
  if( m_solver == nullptr )
  {
    return;
  }

  prepareFluidMetaData();

  if( m_computeRegionStatistics )
    m_aggregator.enableRegionStatistics( meshBodies );

  // if we have to compute CFL numbers later, we need to register additional variables
  if( m_computeCFLNumbers )
    m_solver->registerDataForCFL( meshBodies );

  m_solver->forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                              MeshLevel & mesh,
                                                              string_array const & regionNames )
  {
    prepareLogTableLayouts( mesh.getName() );
    prepareCsvTableLayouts( mesh.getName() );
  } );
}

void StatsTask::prepareFluidMetaData()
{
  ConstitutiveManager const & constitutiveManager = this->getGroupByPath< ConstitutiveManager >( "/Problem/domain/Constitutive" );
  MultiFluidBase const & fluid = constitutiveManager.getGroup< MultiFluidBase >( m_solver->referenceFluidModelName() );

  m_fluid.m_numPhases = fluid.numFluidPhases();
  m_fluid.m_numComps = fluid.numFluidComponents();

  m_fluid.m_phaseNames = fluid.phaseNames();
  m_fluid.m_compNames = fluid.componentNames();

  m_fluid.m_phaseCompNames.resize(
    m_fluid.m_numPhases,
    stdVector< string >( m_fluid.m_numComps, {} ) );

  for( int ip = 0; ip < m_fluid.m_numPhases; ++ip )
    for( int ic = 0; ic < m_fluid.m_numComps; ++ic )
      m_fluid.m_phaseCompNames[ip][ic] = GEOS_FMT( "{}/{}", m_fluid.m_phaseNames[ip], m_fluid.m_compNames[ic] );
}

void StatsTask::prepareLogTableLayouts( string_view meshName )
{
  // only output from rank 0
  if( MpiWrapper::commRank() != 0 )
    return;

  TableLayout const tableLayout = TableLayout()
                                    .setTitle( GEOS_FMT( "{}: mesh {}", getName(), meshName ) );

  m_csvFormatters.emplace( meshName, std::make_unique< TableTextFormatter >( tableLayout ) );
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
        TableLayout::Column().setName( "Time [s]" ),
        TableLayout::Column().setName( "Region" ),   // TODO : mention this change in PR description
        TableLayout::Column().setName( "Min pressure [Pa]" ),
        TableLayout::Column().setName( "Average pressure [Pa]" ),
        TableLayout::Column().setName( "Max pressure [Pa]" ),
        TableLayout::Column().setName( "Min delta pressure [Pa]" ),
        TableLayout::Column().setName( "Max delta pressure [Pa]" ),
        TableLayout::Column().setName( "Min temperature [Pa]" ),
        TableLayout::Column().setName( "Average temperature [Pa]" ),
        TableLayout::Column().setName( "Max temperature [Pa]" ),
        TableLayout::Column().setName( "Total dynamic pore volume [rm^3]" ),
      } );
  addPhaseColumns( tableLayout, "Phase dynamic pore volume", "rm^3", numPhases );
  addPhaseColumns( tableLayout, "Phase mass", massUnit, numPhases );
  addPhaseColumns( tableLayout, "Trapped phase mass (metric 1)", massUnit, numPhases );
  addPhaseColumns( tableLayout, "Non-trapped phase mass (metric 1)", massUnit, numPhases );
  addPhaseColumns( tableLayout, "Immobile phase mass (metric 2)", massUnit, numPhases );
  addPhaseColumns( tableLayout, "Mobile phase mass (metric 2)", massUnit, numPhases );
  addPhaseCompColumns( tableLayout, "Component mass", massUnit, numPhases, numComps );

  auto csvFormatter = std::make_unique< TableCSVFormatter >( tableLayout );
  m_csvFormatters.emplace( meshName, std::move( csvFormatter ) );

  // output CSV header
  std::ofstream outputFile( GEOS_FMT( "{}/{}.csv", m_outputDir, meshName ));
  outputFile << csvFormatter->headerToString();
  GEOS_LOG( GEOS_FMT( "table {} : {}", meshName, csvFormatter->headerToString() ) );     // TODO : remove this log
}

bool StatsTask::execute( real64 const time_n,
                         real64 const dt,
                         integer const GEOS_UNUSED_PARAM( cycleNumber ),
                         integer const GEOS_UNUSED_PARAM( eventCounter ),
                         real64 const GEOS_UNUSED_PARAM( eventProgress ),
                         DomainPartition & domain )
{
  // current time is time_n + dt. TODO: verify implication of events ordering in 'time_n+dt' validity
  real64 statsTime = time_n + dt;

  m_solver->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                          MeshLevel & mesh,
                                                                          string_array const & regionNames )
  {
    if( m_computeRegionStatistics )
    {
      m_aggregator.computeRegionsStatistics( statsTime, mesh, regionNames );
      outputLogStats( statsTime, mesh, regionNames );
      outputCsvStats( statsTime, mesh, regionNames );
    }
  } );

  if( m_computeCFLNumbers )
    m_aggregator.computeCFLNumbers( statsTime, dt, domain );

  return false;
}

void StatsTask::outputLogStats( real64 const statsTime,
                                MeshLevel & mesh,
                                string_array const & regionNames )
{
  auto const formatterIter = m_logFormatters.find( mesh.getName() );
  if( formatterIter==m_logFormatters.end())
    return;

  TableFormatter const & formatter = *formatterIter->second;
  TableData tableData;
  static constexpr auto merge = CellType::MergeNext;

  string_view massUnit = units::getSymbol( m_aggregator.getSolver()->getMassUnit() );
  string_view pressureUnit = units::getSymbol( units::Pressure );
  string_view tempUnit = units::getSymbol( units::Temperature );
  string_view resVolUnit = units::getSymbol( units::ReservoirVolume );

  tableData.addRow( "Statistics time", merge, merge, statsTime );
  m_aggregator.forRegionStatistics( MeshLevel & mesh,
                                    [&]( string_view regionName, RegionStatistics const & stats )
  {

    tableData.addRow( merge, merge, merge, "" );

    tableData.addSeparator();
    tableData.addRow( merge, merge, merge, GEOS_FMT( "Region '{}'", regionName ) );
    tableData.addRow( "statistics", "min", "average", "max" );
    tableData.addSeparator();

    tableData.addRow( GEOS_FMT( "Pressure [{}]", pressureUnit ),
                      stats.minPressure, stats.averagePressure, stats.maxPressure );
    tableData.addRow( GEOS_FMT( "Delta pressure [{}]", pressureUnit ),
                      stats.minDeltaPressure, "/", stats.maxDeltaPressure );
    tableData.addRow( GEOS_FMT( "Temperature [{}]", tempUnit ),
                      stats.minTemperature, stats.averageTemperature, stats.maxTemperature );

    tableData.addSeparator();

    tableData.addRow( GEOS_FMT( "Total dynamic pore volume [{}]", resVolUnit ), CellType::MergeNext, CellType::MergeNext, stats.totalPoreVolume );
    tableData.addRow( GEOS_FMT( "Phase dynamic pore volume [{}]", resVolUnit ),
                      stringutilities::joinLambda( m_fluid.m_phaseNames, "\n", []( auto data ) { return data[0]; } ),
                      CellType::MergeNext,
                      stringutilities::joinLambda( stats.phasePoreVolume, "\n", []( auto data ) { return data[0]; } ) );

    tableData.addSeparator();

    tableData.addRow( GEOS_FMT( "Phase mass [{}]", massUnit ),
                      stringutilities::joinLambda( m_fluid.m_phaseNames, "\n", []( auto data ) { return data[0]; } ),
                      CellType::MergeNext,
                      stringutilities::joinLambda( stats.phaseMass, "\n", []( auto data ) { return data[0]; } ) );

    tableData.addSeparator();

    tableData.addRow( GEOS_FMT( "Trapped phase mass (metric 1) [{}]", massUnit ),
                      stringutilities::joinLambda( m_fluid.m_phaseNames, "\n", []( auto value ) { return value[0]; } ),
                      CellType::MergeNext,
                      stringutilities::joinLambda( stats.trappedPhaseMass, "\n", []( auto value ) { return value[0]; } ) );
    tableData.addRow( GEOS_FMT( "Non-trapped phase mass (metric 1) [{}]", massUnit ),
                      stringutilities::joinLambda( m_fluid.m_phaseNames, "\n", []( auto value ) { return value[0]; } ),
                      CellType::MergeNext,
                      stringutilities::joinLambda( stats.nonTrappedPhaseMass, "\n", []( auto value ) { return value[0]; } ) );

    tableData.addSeparator();

    tableData.addRow( GEOS_FMT( "Immobile phase mass (metric 2) [{}]", massUnit ),
                      stringutilities::joinLambda( m_fluid.m_phaseNames, "\n", []( auto value ) { return value[0]; } ),
                      CellType::MergeNext,
                      stringutilities::joinLambda( stats.immobilePhaseMass, "\n", []( auto value ) { return value[0]; } )  );
    tableData.addRow( GEOS_FMT( "Mobile phase mass (metric 2) [{}]", massUnit ),
                      stringutilities::joinLambda( m_fluid.m_phaseNames, "\n", []( auto value ) { return value[0]; } ),
                      CellType::MergeNext,
                      stringutilities::joinLambda( stats.mobilePhaseMass, "\n", []( auto value ) { return value[0]; } ) );

    tableData.addSeparator();

    tableData.addRow( GEOS_FMT( "Component mass [{}]", massUnit ),
                      stringutilities::join( m_fluid.m_phaseCompNames, '\n' ),
                      CellType::MergeNext,
                      stringutilities::join( stats.componentMass, '\n' ) );

    tableData.addSeparator();
  } );
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        StatsTask,
                        string const &, dataRepository::Group * const )

}   /* namespace compositionalMultiphaseStatistics */

} /* namespace geos */
