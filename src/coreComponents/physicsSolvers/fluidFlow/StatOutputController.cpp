/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file StatOutputController.cpp
 */

#include "StatOutputController.hpp"
#include "events/EventManager.hpp"
#include "events/tasks/TasksManager.hpp"
#include "events/EventBase.hpp"
#include "fileIO/Outputs/TimeHistoryOutput.hpp"

#include "physicsSolvers/solidMechanics/SolidMechanicsStatistics.hpp"
//#include "physicsSolvers/SolverStatistics.hpp"


namespace geos
{

using namespace constitutive;
using namespace dataRepository;

PackCollection * generatePackCollection( TasksManager & taskManager,
                                         string const key,
                                         string_view path,
                                         string_view fieldName )
{
  PackCollection * packCollection = &taskManager.registerGroup< PackCollection >( key );
  string & pcObjectPath = packCollection->getReference< string >( PackCollection::viewKeysStruct::objectPathString());
  pcObjectPath = path;
  string & pcName = packCollection->getReference< string >( PackCollection::viewKeysStruct::fieldNameString());
  pcName = fieldName;
  return packCollection;
}

TimeHistoryOutput * generateTimeHistory( OutputManager & outputManager,
                                         string const key,
                                         string_array sourceTasks,
                                         string_view filename )
{
  TimeHistoryOutput * timeHistory = &outputManager.registerGroup< TimeHistoryOutput >( key );
  string_array & collectorPaths =  timeHistory->getReference< string_array >( TimeHistoryOutput::viewKeys::timeHistoryOutputTargetString() );
  collectorPaths = sourceTasks;
  string & outputFile =  timeHistory->getReference< string >( TimeHistoryOutput::viewKeys::timeHistoryOutputFilenameString() );
  outputFile = filename;
  return timeHistory;
}

StatOutputController::StatOutputController( const string & name,
                                            Group * const parent ):
  TaskBase( name, parent ),
  m_statistics( nullptr )
{}

void StatOutputController::initializePreSubGroups()
{
  Group & problemManager = this->getGroupByPath( "/Problem" );
  DomainPartition & domain = problemManager.getGroup< DomainPartition >( "domain" );
  Group & meshBodies = domain.getMeshBodies();

  TasksManager & taskManager = this->getGroupByPath< TasksManager >( "/Tasks" );
  OutputManager & outputManager = this->getGroupByPath< OutputManager >( "/Outputs" );
  
  std::vector< string > const groupNames = this->getSubGroupsNames();
  m_statistics = &this->getGroup< TaskBase >( groupNames[0] );

  forSubStats( [&]( auto & statistics ) {
    using STATSTYPE = typename TYPEOFREF( statistics );
    statistics.getSolver()->forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                                              MeshLevel & mesh,
                                                                              arrayView1d< string const > const & regionNames )
    {
      ElementRegionManager & elemManager = mesh.getElemManager();
      for( string const & regionName : regionNames )
      {
        ElementRegionBase & region = elemManager.getRegion( regionName );
        string const regionStatPath = GEOS_FMT( "{}/regionStatistics", region.getPath() );
        typename STATSTYPE::RegionStatistics & regionStats = this->getGroupByPath< typename STATSTYPE::RegionStatistics >( regionStatPath );
        string_array sourceTasks;

        regionStats.forWrappers( [&]( WrapperBase const & wrapper )
        { //PackCollection generation
          string const taskManagerKey = GEOS_FMT( "packCollection{}{}", regionName, wrapper.getName());
          PackCollection * packCollection = generatePackCollection( taskManager, taskManagerKey, regionStatPath, wrapper.getName());
          m_packCollections.push_back( packCollection );
          sourceTasks.emplace_back( GEOS_FMT( "{}/", packCollection->getPath()) );
        } );

        { //TimeHistory generation
          string const outputManagerKey = GEOS_FMT( "compFlowHistory{}", regionName );
          string const filename = GEOS_FMT( "generatedHDFStat{}", regionName );
          TimeHistoryOutput * timeHistory = generateTimeHistory( outputManager, outputManagerKey, sourceTasks, filename );
          m_timeHistories.push_back( timeHistory );
        }
      }
    } );
  } );

}

Group * StatOutputController::createChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0( "Adding Statistics: " << childKey << ", " << childName );
  std::unique_ptr< TaskBase > task = TaskBase::CatalogInterface::factory( childKey, childName, this );
  return &this->registerGroup< TaskBase >( childName, std::move( task ) );
}

template< typename LAMBDA >
void StatOutputController::forSubStats( LAMBDA lambda )
{
  forSubGroups< SinglePhaseStatistics,
                CompositionalMultiphaseStatistics >( lambda );
}

void StatOutputController::expandObjectCatalogs()
{
  createChild( SinglePhaseStatistics::catalogName(), SinglePhaseStatistics::catalogName() );
  createChild( CompositionalMultiphaseStatistics::catalogName(), CompositionalMultiphaseStatistics::catalogName() );
}

bool StatOutputController::execute( real64 const time_n,
                                    real64 const dt,
                                    integer const cycleNumber,
                                    integer const eventCounter,
                                    real64 const eventProgress,
                                    DomainPartition & domain )
{
  m_statistics->execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );
  for( PackCollection * packCollection : m_packCollections )
  {
    packCollection->execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );
  }
  for( TimeHistoryOutput * timeHistory : m_timeHistories )
  {
    timeHistory->execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );
  }
  return false;
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        StatOutputController,
                        string const &, dataRepository::Group * const )

}
