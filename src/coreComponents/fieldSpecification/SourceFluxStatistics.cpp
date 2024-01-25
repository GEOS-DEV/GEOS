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

//!\\ TODO add a test for this component in SinglePhase and MultiPhase

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

  registerWrapper( viewKeyStruct::fluxNamesString(), &m_fluxNames ).
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
    // Adding, on each sub-region, a wrapper to hold the stats of each wrapper
    mesh.getElemManager().forElementSubRegions( [&]( ElementSubRegionBase & subRegion )
    {
      for( string const & fluxName : m_fluxNames )
      {
        string const wrapperName = getRegionStatDataName( fluxName );
        Wrapper< WrappedStats > & statsWrapper = subRegion.registerWrapper< WrappedStats >( wrapperName );
        statsWrapper.setRestartFlags( RestartFlags::NO_WRITE );
        statsWrapper.reference().setTarget( getName(), fluxName );
        subRegion.excludeWrappersFromPacking( { wrapperName } );
      }
    } );
    // Do we need to add a similar wrapper for each statistics ? (region-level, mesh-level... to hold as many stats the log levels offer)
  } );
}

SourceFluxStatsAggregator::WrappedStats &
SourceFluxStatsAggregator::getFluxStatData( Group & container,
                                            string_view fluxName )
{
  WrappedStats * r = nullptr;
  container.forWrappers< WrappedStats >( [&]( dataRepository::Wrapper< WrappedStats > & statsWrapper )
  {
    WrappedStats & statsWrapperView = statsWrapper.referenceAsView();
    if( statsWrapperView.getFluxName() == fluxName && statsWrapperView.getAggregatorName() == getName() )
    {
      r = &statsWrapperView;
    }
  } );
  // Error if SourceFluxStatsAggregator::registerDataOnMesh() did not work as expected
  GEOS_ERROR_IF( r == nullptr, GEOS_FMT( "{}: {} data wrongly registered on mesh (no flux stats wrapper was found for {} named {}).",
                                         getName(), catalogName(),
                                         SourceFluxBoundaryCondition::catalogName(), fluxName ) );
  return *r;
}

void SourceFluxStatsAggregator::writeStatData( integer minLogLevel, string_view subSetName,
                                               string_view fluxName, StatData const & stats )
{
  if( getLogLevel() >= minLogLevel && logger::internal::rank == 0 )
  {
    GEOS_LOG_RANK( GEOS_FMT( "{}, {}, {}: applied on {} elements",
                             getName(), fluxName, subSetName, stats.m_elementCount ) );
    GEOS_LOG_RANK( GEOS_FMT( "{}, {}, {}: Produced mass = {} kg",
                             getName(), fluxName, subSetName, stats.m_producedMass ) );
    GEOS_LOG_RANK( GEOS_FMT( "{}, {}, {}: Production rate = {} kg/s",
                             getName(), fluxName, subSetName, stats.m_productionRate ) );
  }
}

bool SourceFluxStatsAggregator::execute( real64 const GEOS_UNUSED_PARAM( time_n ),
                                         real64 const GEOS_UNUSED_PARAM( dt ),
                                         integer const GEOS_UNUSED_PARAM( cycleNumber ),
                                         integer const GEOS_UNUSED_PARAM( eventCounter ),
                                         real64 const GEOS_UNUSED_PARAM( eventProgress ),
                                         DomainPartition & domain )
{
  if( getLogLevel() >= 1 )
  {
    StatData allStats;
    m_solver->forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                              [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & )
    {
      for( string const & fluxName : m_fluxNames )
      {
        StatData fluxStats;

        mesh.getElemManager().forElementRegions( [&]( ElementRegionBase & region )
        {
          StatData regionStats;
          region.forElementSubRegions( [&]( ElementSubRegionBase & subRegion )
          {
            StatData subRegionStats = getFluxStatData( subRegion, fluxName ).finalizePeriod();
            subRegionStats.mpiReduce();

            writeStatData( 4, subRegion.getName(), fluxName, subRegionStats );
            regionStats.combine( subRegionStats );
          } );

          writeStatData( 3, region.getName(), fluxName, regionStats );
          fluxStats.combine( regionStats );
        } );

        writeStatData( 2, "Whole mesh", fluxName, fluxStats );
        allStats.combine( fluxStats );
      }

      writeStatData( 1, "Whole mesh", "Fluxes sum", allStats );
    } );
  }
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
  m_producedMass += MpiWrapper::sum( m_producedMass );
  m_productionRate += MpiWrapper::sum( m_productionRate );
  m_elementCount += MpiWrapper::sum( m_elementCount );
}

void SourceFluxStatsAggregator::WrappedStats::setTarget( string_view aggregatorName,
                                                         string_view fluxName )
{
  m_aggregatorName = aggregatorName;
  m_fluxName = fluxName;
}
void SourceFluxStatsAggregator::WrappedStats::setTimeStepStats( real64 dt,
                                                                real64 productedMass,
                                                                integer elementCount,
                                                                bool overwriteTimeStepStats )
{
  // we are beginning a new timestep, so we must aggregate the pending stats before overwriting the current stats data
  if( !overwriteTimeStepStats )
  {
    m_pendingPeriodMass += m_currentTimeStepMass;
    m_periodDeltaTime += dt;
    m_elementCount = elementCount;
  }

  m_currentTimeStepMass = productedMass;
}
SourceFluxStatsAggregator::StatData SourceFluxStatsAggregator::WrappedStats::finalizePeriod()
{
  real64 periodMass = m_currentTimeStepMass + m_pendingPeriodMass;
  StatData periodStats;
  periodStats.m_producedMass = periodMass;
  periodStats.m_productionRate = periodMass / m_periodDeltaTime;
  periodStats.m_elementCount = m_elementCount;

  m_currentTimeStepMass = 0.0;
  m_pendingPeriodMass = 0.0;
  m_periodDeltaTime = 0.0;
  m_elementCount = 0;

  return periodStats;
}


REGISTER_CATALOG_ENTRY( TaskBase,
                        SourceFluxStatsAggregator,
                        string const &,
                        dataRepository::Group * const )

} /* namespace geos */
