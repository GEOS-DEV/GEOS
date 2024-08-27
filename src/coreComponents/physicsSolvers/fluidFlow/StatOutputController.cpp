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

#include "CompositionalMultiphaseStatistics.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsStatistics.hpp"
#include "SinglePhaseStatistics.hpp"
//#include "physicsSolvers/SolverStatistics.hpp"


namespace geos
{

using namespace constitutive;
using namespace dataRepository;


StatOutputController::StatOutputController( const string & name,
                                            Group * const parent ):
  TaskBase( name, parent )
{}

void StatOutputController::postInputInitialization()
{}
void StatOutputController::initializePreSubGroups()
{
  Group & problemManager = this->getGroupByPath( "/Problem" );
  compMultiphaseStatistics = &this->getGroupByPath< CompositionalMultiphaseStatistics >( "/Tasks/testStats/compflowStatistics" );
  OutputManager & outputManager = this->getGroupByPath< OutputManager >( "/Outputs" );
  TasksManager & taskManager = this->getGroupByPath< TasksManager >( "/Tasks" );
  DomainPartition & domain = problemManager.getGroup< DomainPartition >( "domain" );
  Group & meshBodies = domain.getMeshBodies();

  compMultiphaseStatistics->getSolver()->forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                                                           MeshLevel & mesh,
                                                                                           arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    for( string const & regionName : regionNames )
    {
      ElementRegionBase & region = elemManager.getRegion( regionName );
      string const regionStatPath = GEOS_FMT( "{}/regionStatistics", region.getPath() );
      CompositionalMultiphaseStatistics::RegionStatistics & regionStats =
        this->getGroupByPath< CompositionalMultiphaseStatistics::RegionStatistics >( regionStatPath );

      string_array sourceTasks;
      regionStats.forWrappers( [&]( WrapperBase const & wrapper )
      {
        //PackCollection generation
        string const taskManagerKey = GEOS_FMT( "packCollection{}{}", regionName, wrapper.getName());
        PackCollection * packCollection = &taskManager.registerGroup< PackCollection >( taskManagerKey );
        string & pcObjectPath = packCollection->getReference< string >( PackCollection::viewKeysStruct::objectPathString());
        pcObjectPath=  regionStatPath;
        string & pcName = packCollection->getReference< string >( PackCollection::viewKeysStruct::fieldNameString());
        pcName = wrapper.getName();
        m_packCollections.push_back( packCollection );

        sourceTasks.emplace_back( GEOS_FMT( "{}/", packCollection->getPath()) );
      } );

      { //TimeHistory generation
        string const outputManagerKey = GEOS_FMT( "compFlowHistory{}", regionName );
        TimeHistoryOutput * timeHistory = &outputManager.registerGroup< TimeHistoryOutput >( outputManagerKey );
        string_array & collectorPaths =  timeHistory->getReference< string_array >( TimeHistoryOutput::viewKeys::timeHistoryOutputTargetString() );
        collectorPaths = sourceTasks;
        string & outputFile =  timeHistory->getReference< string >( TimeHistoryOutput::viewKeys::timeHistoryOutputFilenameString() );
        outputFile =  GEOS_FMT( "generatedHDFStat{}", regionName );
        m_timeHistories.push_back( timeHistory );
        // std::cout << "TEST PATH -- " << joinPath( FunctionBase::getOutputDirectory(), regionName + ".hdf5" ) << std::endl;
      }
    }
  } );
}

Group * StatOutputController::createChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0( "Adding Statistics: " << childKey << ", " << childName );
  std::unique_ptr< TaskBase > task = TaskBase::CatalogInterface::factory( childKey, childName, this );
  return &this->registerGroup< TaskBase >( childName, std::move( task ) );
}


void StatOutputController::expandObjectCatalogs()
{
  createChild( CompositionalMultiphaseStatistics::catalogName(), CompositionalMultiphaseStatistics::catalogName() );
  // createChild( SolidMechanicsStatistics::catalogName(), SolidMechanicsStatistics::catalogName() );
  // createChild( SinglePhaseStatistics::catalogName(), SinglePhaseStatistics::catalogName() );
  //createChild( SolverBase::groupKeyStruct::solverStatisticsString(), SolverBase::groupKeyStruct::solverStatisticsString() );
}

bool StatOutputController::execute( real64 const time_n,
                                    real64 const dt,
                                    integer const cycleNumber,
                                    integer const eventCounter,
                                    real64 const eventProgress,
                                    DomainPartition & domain )
{
  compMultiphaseStatistics->execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );
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
