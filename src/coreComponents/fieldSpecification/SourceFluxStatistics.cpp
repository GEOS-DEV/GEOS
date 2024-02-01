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
  // the fields have to be registered in "registerDataOnMesh" (and not later)
  // otherwise they cannot be targeted by TimeHistory

  // for now, this guard is needed to avoid breaking the xml schema generation
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

// SourceFluxStatsAggregator::WrappedStats &
// SourceFluxStatsAggregator::getFluxStatData( Group & container,
//                                             string_view fluxName )
// {
//   WrappedStats * r = nullptr;
//   container.forWrappers< WrappedStats >( [&]( dataRepository::Wrapper< WrappedStats > & statsWrapper )
//   {
//     WrappedStats & statsWrapperView = statsWrapper.referenceAsView();
//     if( statsWrapperView.getFluxName() == fluxName && statsWrapperView.getAggregatorName() == getName() )
//     {
//       r = &statsWrapperView;
//     }
//   } );
//   // Error if SourceFluxStatsAggregator::registerDataOnMesh() did not work as expected
//   GEOS_ERROR_IF( r == nullptr, GEOS_FMT( "{}: {} data wrongly registered on mesh (no flux stats wrapper was found for {} named {}).",
//                                          getName(), catalogName(),
//                                          SourceFluxBoundaryCondition::catalogName(), fluxName ) );
//   return *r;
// }

void SourceFluxStatsAggregator::writeStatData( integer minLogLevel,
                                               string_view elementSetName,
                                               WrappedStats const & wrappedStats )
{
  if( getLogLevel() >= minLogLevel && logger::internal::rank == 0 )
  {
    GEOS_LOG_RANK( GEOS_FMT( "{} {} (of {}, in {}): Producting on {} elements",
                             catalogName(), getName(), wrappedStats.getFluxName(), elementSetName,
                             wrappedStats.stats().m_elementCount ) );
    GEOS_LOG_RANK( GEOS_FMT( "{} {} (of {}, in {}): Produced mass = {} kg",
                             catalogName(), getName(), wrappedStats.getFluxName(), elementSetName,
                             wrappedStats.stats().m_producedMass ) );
    GEOS_LOG_RANK( GEOS_FMT( "{} {} (of {}, in {}): Production rate = {} kg/s",
                             catalogName(), getName(), wrappedStats.getFluxName(), elementSetName,
                             wrappedStats.stats().m_productionRate ) );
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



void SourceFluxStatsAggregator::StatData::combine( StatData const & other )
{
  m_producedMass += other.m_producedMass;
  m_productionRate += other.m_productionRate;
  m_elementCount += other.m_elementCount;
}
void SourceFluxStatsAggregator::StatData::mpiReduce()
{
  m_producedMass = MpiWrapper::sum( m_producedMass );
  m_productionRate = MpiWrapper::sum( m_productionRate );
  m_elementCount = MpiWrapper::sum( m_elementCount );
}

void SourceFluxStatsAggregator::WrappedStats::setTarget( string_view aggregatorName,
                                                         string_view fluxName )
{
  m_aggregatorName = aggregatorName;
  m_fluxName = fluxName;
}
void SourceFluxStatsAggregator::WrappedStats::gatherTimeStepStats( real64 dt,
                                                                   real64 productedMass,
                                                                   integer elementCount,
                                                                   bool overwriteTimeStepStats )
{
  // we are beginning a new timestep, so we must aggregate the pending stats (mass & dt) before collecting the current stats data
  if( !overwriteTimeStepStats )
  {
    m_periodStats.m_periodPendingMass += m_periodStats.m_timeStepMass;
    m_periodStats.m_timeStepMass = 0;
    m_periodStats.m_periodDeltaTime += dt;
    m_periodStats.m_elementCount = elementCount;
  }

  m_periodStats.m_timeStepMass = productedMass;
}
void SourceFluxStatsAggregator::WrappedStats::finalizePeriod()
{
  // produce timestep stats of this ranks
  real64 periodMass = m_periodStats.m_timeStepMass + m_periodStats.m_periodPendingMass;
  m_stats.m_producedMass = periodMass;
  m_stats.m_productionRate = m_periodStats.m_periodDeltaTime > 0.0 ?
                             periodMass / m_periodStats.m_periodDeltaTime :
                             0.0;
  m_stats.m_elementCount = m_periodStats.m_elementCount;

  // combine period results from all MPI ranks
  m_stats.mpiReduce();

  // start a new timestep
  m_periodStats = WrappedStats::PeriodStats();

}


REGISTER_CATALOG_ENTRY( TaskBase,
                        SourceFluxStatsAggregator,
                        string const &,
                        dataRepository::Group * const )

} /* namespace geos */
