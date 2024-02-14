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
 * @file SourceFluxStatistics.cpp
 */

#include "SourceFluxStatistics.hpp"

#include "SourceFluxBoundaryCondition.hpp"
#include "FieldSpecificationManager.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{
using namespace dataRepository;

SourceFluxStatsAggregator::SourceFluxStatsAggregator( const string & name,
                                                      Group * const parent ):
  Base( name, parent )
{
  getWrapper< integer >( Group::viewKeyStruct::logLevelString() ).
    appendDescription( GEOS_FMT( "\n- Log Level 1 outputs the sum of all {0}(s) produced rate & mass,\n"
                                 "- Log Level 2 details values for each {0},\n"
                                 "- Log Level 3 details values for each region,\n"
                                 "- Log Level 4 details values for each sub-region.",
                                 SourceFluxBoundaryCondition::catalogName() ) );

  registerWrapper( viewKeyStruct::fluxNamesString().data(), &m_fluxNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( GEOS_FMT( "Name(s) array of the {0}(s) for which we want the statistics. "
                              "Use \"all\" to target all {0}.",
                              SourceFluxBoundaryCondition::catalogName() ) );
}

void SourceFluxStatsAggregator::postProcessInput()
{
  Base::postProcessInput();

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  if( m_fluxNames.size() == 1 && m_fluxNames[0] == "all" )
  {
    m_fluxNames.clear();
    fsManager.forSubGroups< SourceFluxBoundaryCondition >( [&]( SourceFluxBoundaryCondition & sourceFlux )
    {
      m_fluxNames.emplace_back( string( sourceFlux.getName() ) );
    } );
    GEOS_WARNING_IF( m_fluxNames.empty(),
                     GEOS_FMT( "{}: No {} was found in {}.",
                               getDataContext(), SourceFluxBoundaryCondition::catalogName(),
                               fsManager.getDataContext() ) );
  }
  else
  {
    for( string const & fluxName : m_fluxNames )
    {
      GEOS_ERROR_IF( !fsManager.hasGroup< SourceFluxBoundaryCondition >( fluxName ),
                     GEOS_FMT( "{}: No {} named {} was found in {}.",
                               getDataContext(), SourceFluxBoundaryCondition::catalogName(),
                               fluxName, fsManager.getDataContext() ) );
    }
  }
}

Wrapper< SourceFluxStatsAggregator::WrappedStats > &
SourceFluxStatsAggregator::registerWrappedStats( Group & group, string_view fluxName )
{
  string const wrapperName = getStatWrapperName( fluxName );
  Wrapper< WrappedStats > & statsWrapper = group.registerWrapper< WrappedStats >( wrapperName );
  statsWrapper.setRestartFlags( RestartFlags::NO_WRITE );
  statsWrapper.reference().setTarget( getName(), fluxName );
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
                                                              arrayView1d< string const > const & )
  {
    registerWrappedStats( mesh, viewKeyStruct::fluxSetWrapperString() );
    for( string const & fluxName : m_fluxNames )
    {
      registerWrappedStats( mesh, fluxName );

      mesh.getElemManager().forElementRegions( [&]( ElementRegionBase & region )
      {
        Wrapper< WrappedStats > & regionStatsWrapper = registerWrappedStats( region, fluxName );
        region.excludeWrappersFromPacking( { regionStatsWrapper.getName() } );

        region.forElementSubRegions( [&]( ElementSubRegionBase & subRegion )
        {
          Wrapper< WrappedStats > & subRegionStatsWrapper = registerWrappedStats( subRegion, fluxName );
          subRegion.excludeWrappersFromPacking( { subRegionStatsWrapper.getName() } );
        } );
      } );
    }
  } );
}

void SourceFluxStatsAggregator::writeStatData( integer minLogLevel,
                                               string_view elementSetName,
                                               WrappedStats const & wrappedStats )
{
  if( getLogLevel() >= minLogLevel && logger::internal::rank == 0 )
  {
    GEOS_LOG_RANK( GEOS_FMT( "{} {} (of {}, in {}): Producing on {} elements",
                             catalogName(), getName(), wrappedStats.getFluxName(), elementSetName,
                             wrappedStats.stats().m_elementCount ) );

    // we want to format differently if we have got multiple phases or not
    if( wrappedStats.stats().m_producedMass.size() == 1 )
    {
      GEOS_LOG_RANK( GEOS_FMT( "{} {} (of {}, in {}): Produced mass = {} kg",
                               catalogName(), getName(), wrappedStats.getFluxName(), elementSetName,
                               wrappedStats.stats().m_producedMass[0] ) );
      GEOS_LOG_RANK( GEOS_FMT( "{} {} (of {}, in {}): Production rate = {} kg/s",
                               catalogName(), getName(), wrappedStats.getFluxName(), elementSetName,
                               wrappedStats.stats().m_productionRate[0] ) );
    }
    else
    {
      GEOS_LOG_RANK( GEOS_FMT( "{} {} (of {}, in {}): Produced mass = {} kg",
                               catalogName(), getName(), wrappedStats.getFluxName(), elementSetName,
                               wrappedStats.stats().m_producedMass ) );
      GEOS_LOG_RANK( GEOS_FMT( "{} {} (of {}, in {}): Production rate = {} kg/s",
                               catalogName(), getName(), wrappedStats.getFluxName(), elementSetName,
                               wrappedStats.stats().m_productionRate ) );
    }
  }
}

bool SourceFluxStatsAggregator::execute( real64 const GEOS_UNUSED_PARAM( time_n ),
                                         real64 const GEOS_UNUSED_PARAM( dt ),
                                         integer const GEOS_UNUSED_PARAM( cycleNumber ),
                                         integer const GEOS_UNUSED_PARAM( eventCounter ),
                                         real64 const GEOS_UNUSED_PARAM( eventProgress ),
                                         DomainPartition & domain )
{
  forMeshLevelStatsWrapper( domain,
                            [&] ( MeshLevel & meshLevel, WrappedStats & meshLevelStats )
  {
    meshLevelStats.stats() = StatData();

    forAllFluxStatsWrappers( meshLevel,
                             [&] ( MeshLevel &, WrappedStats & fluxStats )
    {
      fluxStats.stats() = StatData();

      forAllRegionStatsWrappers( meshLevel, fluxStats.getFluxName(),
                                 [&] ( ElementRegionBase & region, WrappedStats & regionStats )
      {
        regionStats.stats() = StatData();

        forAllSubRegionStatsWrappers( region, regionStats.getFluxName(),
                                      [&] ( ElementSubRegionBase & subRegion, WrappedStats & subRegionStats )
        {
          subRegionStats.finalizePeriod();

          regionStats.stats().combine( subRegionStats.stats() );
          writeStatData( 4, subRegion.getName(), subRegionStats );
        } );

        fluxStats.stats().combine( regionStats.stats() );
        writeStatData( 3, region.getName(), regionStats );
      } );

      meshLevelStats.stats().combine( fluxStats.stats() );
      writeStatData( 2, viewKeyStruct::allRegionWrapperString(), fluxStats );
    } );

    writeStatData( 1, viewKeyStruct::allRegionWrapperString(), meshLevelStats );
  } );

  return false;
}



void SourceFluxStatsAggregator::StatData::allocate( integer phaseCount )
{
  if( m_producedMass.size() != phaseCount )
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

  // if beginning a new timestep, we must aggregate the stats from previous timesteps (mass & dt) before collecting the new ones
  bool isBeginingNewTS = currentTime >= ( m_periodStats.m_timeStepStart + m_periodStats.m_timeStepDeltaTime );
  if( isBeginingNewTS )
  {
    for( int ip = 0; ip < m_periodStats.getPhaseCount(); ++ip )
    {
      m_periodStats.m_periodPendingMass[ip] += m_periodStats.m_timeStepMass[ip];
    }
    m_periodStats.m_elementCount = elementCount;
    m_periodStats.m_periodPendingDeltaTime += m_periodStats.m_timeStepDeltaTime;
  }

  // new timestep stats to take into account
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

  // produce timestep stats of this ranks
  m_stats.m_elementCount = m_periodStats.m_elementCount;

  real64 const dt = m_periodStats.m_timeStepDeltaTime + m_periodStats.m_periodPendingDeltaTime;
  real64 const timeDivisor = dt > 0.0 ? 1.0 / dt : 0.0;
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
void SourceFluxStatsAggregator::WrappedStats::PeriodStats::allocate( integer phaseCount )
{
  if( m_timeStepMass.size() != phaseCount )
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
}


REGISTER_CATALOG_ENTRY( TaskBase,
                        SourceFluxStatsAggregator,
                        string const &,
                        dataRepository::Group * const )

} /* namespace geos */
