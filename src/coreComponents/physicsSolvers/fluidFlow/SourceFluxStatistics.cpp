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
 * @file SourceFluxStatistics.cpp
 */

#include "SourceFluxStatistics.hpp"

#include "fieldSpecification/SourceFluxBoundaryCondition.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "physicsSolvers/fluidFlow/LogLevelsInfo.hpp"

namespace geos
{
using namespace dataRepository;

SourceFluxStatsAggregator::SourceFluxStatsAggregator( const string & name,
                                                      Group * const parent ):
  Base( name, parent )
{
  registerWrapper( viewKeyStruct::fluxNamesString().data(), &m_fluxNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDefaultValue( "*" ).
    setDescription( GEOS_FMT( "Name(s) array of the {0}(s) for which we want the statistics. "
                              "Use \"*\" to target all {0}.",
                              SourceFluxBoundaryCondition::catalogName() ) );

  addLogLevel< logInfo::DetailedSourceFluxStats >();
  addLogLevel< logInfo::DetailedRegionsSourceFluxStats >();
  addLogLevel< logInfo::AggregatedSourceFluxStats >();
}

void SourceFluxStatsAggregator::postInputInitialization()
{
  Base::postInputInitialization();

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  if( m_fluxNames.size() == 1 && m_fluxNames[0] == "*" )
  {
    m_fluxNames.clear();
    fsManager.forSubGroups< SourceFluxBoundaryCondition >( [&]( SourceFluxBoundaryCondition & sourceFlux )
    {
      m_fluxNames.emplace_back( string( sourceFlux.getName() ) );
    } );
    GEOS_WARNING_IF( m_fluxNames.empty(),
                     GEOS_FMT( "No {} was found in {}.",
                               SourceFluxBoundaryCondition::catalogName(), fsManager.getName() ),
                     getDataContext(), fsManager.getDataContext(), getDataContext() );
  }
  else
  {
    for( string const & fluxName : m_fluxNames )
    {
      GEOS_ERROR_IF( !fsManager.hasGroup< SourceFluxBoundaryCondition >( fluxName ),
                     GEOS_FMT( "No {} named {} was found in {}.",
                               SourceFluxBoundaryCondition::catalogName(),
                               fluxName, fsManager.getName() ),
                     getDataContext(), fsManager.getDataContext(), getDataContext()  );
    }
  }

  m_csvFilename = GEOS_FMT( "{}/{}_fluxes_stats.csv",
                            m_outputDir, getName() );
}

Wrapper< SourceFluxStatsAggregator::WrappedStats > &
SourceFluxStatsAggregator::registerWrappedStats( Group & group,
                                                 string_view fluxName )
{
  string const wrapperName = getStatWrapperName( fluxName );
  Wrapper< WrappedStats > & statsWrapper = group.registerWrapper< WrappedStats >( wrapperName );
  statsWrapper.setRestartFlags( RestartFlags::NO_WRITE );

  WrappedStats & stats = statsWrapper.reference();
  stats.setTarget( getName(), fluxName );

  { //tableLayout initialisation
    string_view massUnit = units::getSymbol( m_solver->getMassUnit() );
    string const massColumn = GEOS_FMT( "Prod. mass [{}]", massUnit );
    string const rateColumn = GEOS_FMT( "Prod. rate [{}/s]", massUnit );

    // log output preparation
    TableLayout statsLogLayout( "", { "Flux(es)", "Region", "Element Count", massColumn, rateColumn } );
    m_logLayout = statsLogLayout;

    // CSV output preparation
    TableLayout const statsCSVLayout( "", {"Time [s]", "Flux(es)", "Region", "Element Count", massColumn, rateColumn} );
    m_csvLayout = statsCSVLayout;
  }

  return statsWrapper;
}
void SourceFluxStatsAggregator::registerDataOnMesh( Group & meshBodies )
{
  if( m_solver == nullptr )
  {
    return;
  }

  m_solver->forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                              MeshLevel & mesh,
                                                              string_array const & )
  {
    registerWrappedStats( mesh,
                          viewKeyStruct::fluxSetWrapperString() );

    for( string const & fluxName : m_fluxNames )
    {
      registerWrappedStats( mesh,
                            fluxName );

      mesh.getElemManager().forElementRegions( [&]( ElementRegionBase & region )
      {
        Wrapper< WrappedStats > & regionStatsWrapper =
          registerWrappedStats( region, fluxName );
        region.excludeWrappersFromPacking( { regionStatsWrapper.getName() } );

        region.forElementSubRegions( [&]( ElementSubRegionBase & subRegion )
        {
          Wrapper< WrappedStats > & subRegionStatsWrapper =
            registerWrappedStats( subRegion, fluxName );
          subRegion.excludeWrappersFromPacking( { subRegionStatsWrapper.getName() } );
        } );
      } );
    }
  } );

  if( m_writeCSV > 0 && MpiWrapper::commRank() == 0 )
  {
    std::ofstream outputFile( m_csvFilename );
    TableCSVFormatter const tableStatFormatter( m_csvLayout );
    outputFile << tableStatFormatter.headerToString();
    outputFile.close();
  }
}

void SourceFluxStatsAggregator::gatherStatsForLog( bool logLevelActive,
                                                   string_view fluxName,
                                                   string_view elementSetName,
                                                   TableData & tableData,
                                                   WrappedStats const & wrappedStats )
{
  if( logLevelActive && MpiWrapper::commRank()  == 0 )
  {
    if( wrappedStats.stats().m_producedMass.size() == 1 )
    {
      tableData.addRow( fluxName,
                        elementSetName,
                        GEOS_FMT( "{}", wrappedStats.stats().m_elementCount ),
                        GEOS_FMT( "{}", wrappedStats.stats().m_producedMass[0] ),
                        GEOS_FMT( "{}", wrappedStats.stats().m_productionRate[0] ) );
    }
    else
    {
      tableData.addRow( fluxName,
                        elementSetName,
                        GEOS_FMT( "{}", wrappedStats.stats().m_elementCount ),
                        GEOS_FMT( "{}", wrappedStats.stats().m_producedMass ),
                        GEOS_FMT( "{}", wrappedStats.stats().m_productionRate ) );
    }
  }
}

void SourceFluxStatsAggregator::gatherStatsForCSV( string_view fluxName,
                                                   string_view elementSetName,
                                                   TableData & tableData,
                                                   WrappedStats const & wrappedStats )
{
  if( m_writeCSV > 0 && MpiWrapper::commRank() == 0 )
  {
    if( wrappedStats.stats().m_producedMass.size() == 1 )
    {
      tableData.addRow( wrappedStats.getStatsPeriodStart(),
                        fluxName,
                        elementSetName,
                        wrappedStats.stats().m_elementCount,
                        wrappedStats.stats().m_producedMass[0],
                        wrappedStats.stats().m_productionRate[0] );
    }
    else
    {
      tableData.addRow( wrappedStats.getStatsPeriodStart(),
                        fluxName,
                        elementSetName,
                        wrappedStats.stats().m_elementCount,
                        wrappedStats.stats().m_producedMass,
                        wrappedStats.stats().m_productionRate );
    }
  }
}

void SourceFluxStatsAggregator::outputStatsToLog( bool logLevelActive,
                                                  string_view fluxesStr,
                                                  TableData const & statsData )
{
  GEOS_UNUSED_VAR( statsData );
  if( logLevelActive && MpiWrapper::commRank() == 0 )
  {
    string const title = GEOS_FMT( "{}, flux statistics for: {}",
                                   getName(), fluxesStr );
    string const multilineTitle = wrapTextToMaxLength( title, 80 );
    m_logLayout.setTitle( multilineTitle );
    TableTextFormatter const tableStatFormatter( m_logLayout );
    GEOS_LOG_RANK( tableStatFormatter.toString( statsData ) );
  }
}
void SourceFluxStatsAggregator::outputStatsToCSV( TableData & csvData )
{
  if( m_writeCSV > 0 && MpiWrapper::commRank() == 0 )
  {
    std::ofstream outputFile( m_csvFilename, std::ios::app );
    TableCSVFormatter const tableStatFormatter( m_csvLayout );
    outputFile << tableStatFormatter.dataToString( csvData );
    outputFile.close();
    csvData.clear();
  }
}

bool SourceFluxStatsAggregator::execute( real64 const GEOS_UNUSED_PARAM( time_n ),
                                         real64 const GEOS_UNUSED_PARAM( dt ),
                                         integer const GEOS_UNUSED_PARAM( cycleNumber ),
                                         integer const GEOS_UNUSED_PARAM( eventCounter ),
                                         real64 const GEOS_UNUSED_PARAM( eventProgress ),
                                         DomainPartition & domain )
{
  bool const fluxMeshesStatsOn = isLogLevelActive< logInfo::AggregatedSourceFluxStats >( getLogLevel() );
  bool const fluxesStatsOn = isLogLevelActive< logInfo::DetailedSourceFluxStats >( getLogLevel() );
  bool const regionsStatsOn = isLogLevelActive< logInfo::DetailedRegionsSourceFluxStats >( getLogLevel() );
  string_view allRegionsStr = viewKeyStruct::allRegionWrapperString();
  string_view fluxesSetStr = viewKeyStruct::fluxSetWrapperString();

  forMeshLevelStatsWrapper( domain,
                            [&] ( MeshLevel & meshLevel, WrappedStats & meshLevelStats )
  {
    TableData logData;
    TableData csvData;
    meshLevelStats.stats().reset();
    forAllFluxStatsWrappers( meshLevel,
                             [&] ( MeshLevel &, WrappedStats & fluxStats )
    {
      fluxStats.stats().reset();
      forAllRegionStatsWrappers( meshLevel, fluxStats.getFluxName(),
                                 [&] ( ElementRegionBase & region, WrappedStats & regionStats )
      {
        regionStats.stats().reset();

        forAllSubRegionStatsWrappers( region, regionStats.getFluxName(),
                                      [&] ( ElementSubRegionBase &, WrappedStats & subRegionStats )
        {
          subRegionStats.stats().reset();
          subRegionStats.finalizePeriod();
          regionStats.combine( subRegionStats );
        } );
        fluxStats.combine( regionStats );

        gatherStatsForLog( regionsStatsOn,
                           fluxStats.getFluxName(), region.getName(), logData, regionStats );
        gatherStatsForCSV( fluxStats.getFluxName(), region.getName(), csvData, regionStats );
      } );

      meshLevelStats.combine( fluxStats );

      gatherStatsForLog( fluxesStatsOn,
                         fluxStats.getFluxName(), allRegionsStr, logData, fluxStats );
      gatherStatsForCSV( fluxStats.getFluxName(), allRegionsStr, csvData, fluxStats );
    } );
    gatherStatsForLog( fluxMeshesStatsOn,
                       fluxesSetStr, allRegionsStr, logData, meshLevelStats );
    gatherStatsForCSV( fluxesSetStr, allRegionsStr, csvData, meshLevelStats );

    outputStatsToLog( fluxMeshesStatsOn,
                      stringutilities::join( m_fluxNames, ", " ), logData );
    outputStatsToCSV( csvData );
  } );

  return false;
}



void SourceFluxStatsAggregator::StatData::allocate( integer phaseCount )
{
  if( m_producedMass.size() < phaseCount )
  {
    m_producedMass.resize( phaseCount );
    m_productionRate.resize( phaseCount );
  }
}
void SourceFluxStatsAggregator::StatData::reset()
{
  for( int ip = 0; ip < getPhaseCount(); ++ip )
  {
    m_producedMass[ip] = 0.0;
    m_productionRate[ip] = 0.0;
  }
  m_elementCount = 0;
}
void SourceFluxStatsAggregator::StatData::combine( StatData const & other )
{
  allocate( other.getPhaseCount() );

  for( int ip = 0; ip < other.getPhaseCount(); ++ip )
  {
    m_producedMass[ip] += other.m_producedMass[ip];
    m_productionRate[ip] += other.m_productionRate[ip];
  }
  m_elementCount += other.m_elementCount;
}
void SourceFluxStatsAggregator::StatData::mpiReduce()
{
  allocate( MpiWrapper::max( m_producedMass.size() ) );
  for( int ip = 0; ip < getPhaseCount(); ++ip )
  {
    m_producedMass[ip] = MpiWrapper::sum( m_producedMass[ip] );
    m_productionRate[ip] = MpiWrapper::sum( m_productionRate[ip] );
  }
  m_elementCount = MpiWrapper::sum( m_elementCount );
}

void SourceFluxStatsAggregator::WrappedStats::setTarget( string_view aggregatorName,
                                                         string_view fluxName )
{
  m_aggregatorName = aggregatorName;
  m_fluxName = fluxName;
}
void SourceFluxStatsAggregator::WrappedStats::gatherTimeStepStats( real64 const currentTime, real64 const dt,
                                                                   arrayView1d< real64 const > const & producedMass,
                                                                   integer const elementCount )
{
  m_periodStats.allocate( producedMass.size() );

  if( !m_periodStats.m_isGathering )
  {
    // if beginning a new period, we must initialize constant values over the period
    m_periodStats.m_periodStart = currentTime;
    m_periodStats.m_elementCount = elementCount;
    m_periodStats.m_isGathering = true;
  }
  else
  {
    GEOS_WARNING_IF( currentTime< m_periodStats.m_timeStepStart, GEOS_FMT( "{}: Time seems to have rollback, stats will be wrong.", m_aggregatorName ) );
    if( currentTime > m_periodStats.m_timeStepStart )
    {
      // if beginning a new timestep, we must accumulate the stats from previous timesteps (mass & dt) before collecting the new ones
      for( int ip = 0; ip < m_periodStats.getPhaseCount(); ++ip )
      {
        m_periodStats.m_periodPendingMass[ip] += m_periodStats.m_timeStepMass[ip];
      }
      m_periodStats.m_periodPendingDeltaTime += m_periodStats.m_timeStepDeltaTime;
    }
  }
  // current timestep stats to take into account (overriding if not begining a new timestep)
  m_periodStats.m_timeStepStart = currentTime;
  m_periodStats.m_timeStepDeltaTime = dt;
  for( int ip = 0; ip < m_periodStats.getPhaseCount(); ++ip )
  {
    m_periodStats.m_timeStepMass = producedMass;
  }
}
void SourceFluxStatsAggregator::WrappedStats::finalizePeriod()
{
  // init phase data memory allocation if needed
  m_stats.allocate( m_periodStats.getPhaseCount() );

  // produce the period stats of this rank
  m_stats.m_elementCount = m_periodStats.m_elementCount;
  // MPI ranks without cells in the flux region may bypass gatherTimeStepStats() entirely,
  // leaving their m_periodStart behind. Take the max so all ranks agree on the most advanced time.
  m_statsPeriodStart = MpiWrapper::max( m_periodStats.m_periodStart );
  m_statsPeriodDT = m_periodStats.m_timeStepDeltaTime + m_periodStats.m_periodPendingDeltaTime;

  real64 const timeDivisor = m_statsPeriodDT > 0.0 ? 1.0 / m_statsPeriodDT : 0.0;
  for( int ip = 0; ip < m_periodStats.getPhaseCount(); ++ip )
  {
    real64 periodMass = m_periodStats.m_timeStepMass[ip] + m_periodStats.m_periodPendingMass[ip];
    m_stats.m_producedMass[ip] = periodMass;
    m_stats.m_productionRate[ip] = periodMass * timeDivisor;
  }

  // combine period results from all MPI ranks
  m_stats.mpiReduce();

  // start a new timestep
  m_periodStats.reset();
}
void SourceFluxStatsAggregator::WrappedStats::combine( WrappedStats const & other )
{
  stats().combine( other.stats() );
  // Some fluxes may lag behind (their ranks may have skipped gatherTimeStepStats()),
  // so take the most advanced period start time across all combined stats.
  m_statsPeriodStart = LvArray::math::max( m_statsPeriodStart, other.m_statsPeriodStart );
}
void SourceFluxStatsAggregator::WrappedStats::PeriodStats::allocate( integer phaseCount )
{
  if( m_timeStepMass.size() < phaseCount )
  {
    m_timeStepMass.resize( phaseCount );
    m_periodPendingMass.resize( phaseCount );
  }
}
void SourceFluxStatsAggregator::WrappedStats::PeriodStats::reset()
{
  for( int ip = 0; ip < getPhaseCount(); ++ip )
  {
    m_timeStepMass[ip] = 0.0;
    m_periodPendingMass[ip] = 0.0;
  }
  m_periodPendingDeltaTime = 0.0;
  m_elementCount = 0;
  m_timeStepStart = 0.0;
  m_timeStepDeltaTime = 0.0;
  m_isGathering = false;
}


REGISTER_CATALOG_ENTRY( TaskBase,
                        SourceFluxStatsAggregator,
                        string const &,
                        dataRepository::Group * const )

} /* namespace geos */
