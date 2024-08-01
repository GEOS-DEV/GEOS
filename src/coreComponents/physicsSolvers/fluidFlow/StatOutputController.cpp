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
#include "mainInterface/ProblemManager.hpp"
#include "events/EventManager.hpp"
#include "events/tasks/TasksManager.hpp"
#include "events/EventBase.hpp"

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
  TaskBase( name, parent ),
  m_timeHistory( nullptr )
{}

void StatOutputController::postInputInitialization()
{}

void StatOutputController::registerDataOnMesh( Group & meshBodies )
{

  CompositionalMultiphaseStatistics & compMultiStats =
    this->getGroupByPath< CompositionalMultiphaseStatistics >( "/Tasks/testStats/compflowStatistics" );

  // for now, this guard is needed to avoid breaking the xml schema generation
  if( compMultiStats.getSolver() == nullptr )
  {
    return;
  }

  OutputManager & outputManager = this->getGroupByPath< OutputManager >( "/Outputs" );
  //EventManager & eventManager = this->getGroupByPath< EventManager >( "/Events" );

  compMultiStats.getSolver()->forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                                                MeshLevel & GEOS_UNUSED_PARAM( mesh ),
                                                                                arrayView1d< string const > const & regionNames )
  {
    TasksManager & taskManager = this->getGroupByPath< TasksManager >( "/Tasks" );
    std::vector< string > timeHistoryStatsSource;
    for( string const & regionName : regionNames )
    {
      CompositionalMultiphaseStatistics::RegionStatistics & regionStats =
        this->getGroupByPath< CompositionalMultiphaseStatistics::RegionStatistics >( regionName );
      regionStats.forWrappers( [&]( WrapperBase const & wrapper )
      {
        string keyGroup = GEOS_FMT( "packCollection{}", wrapper.getName());

        //PackCollection generation
        PackCollection * packCollection = &taskManager.registerGroup< PackCollection >( keyGroup );
        packCollection->getWrapper< string >( PackCollection::viewKeysStruct::objectPathString()).
          setApplyDefaultValue( GEOS_FMT( "{}{}", regionStats.getPath(), "/regionStatistics" ) );
        packCollection->getWrapper< string >( PackCollection::viewKeysStruct::fieldNameString()).
          setApplyDefaultValue( wrapper.getName());
        m_packCollections.push_back( packCollection );

        string packFullPath = GEOS_FMT( "{}{}", packCollection->getPath(), packCollection->getName());
        timeHistoryStatsSource.push_back( packFullPath );
      }
                               );
    }
  } );

  //TimeHistory Generation
  m_timeHistory =  &outputManager.registerGroup< TimeHistoryOutput >( "generatedStatsHistoryOutput" );
  m_timeHistory->getWrapper< string >( TimeHistoryOutput::viewKeys::timeHistoryOutputTargetString()).
    setApplyDefaultValue( "4.32e6" );
  m_timeHistory->getWrapper< string >( TimeHistoryOutput::viewKeys::timeHistoryOutputFilenameString()).
    setApplyDefaultValue( "generatedHDFStat" );

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
  std::cout << " Execute " << std::endl;

  for( PackCollection * packCollection : m_packCollections )
  {
    packCollection->execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );
  }
  m_timeHistory->execute( time_n, dt, cycleNumber, eventCounter, eventProgress, domain );
  return false;
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        StatOutputController,
                        string const &, dataRepository::Group * const )

}
