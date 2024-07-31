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
#include "fileIO/timeHistory/PackCollection.hpp"

#include "CompositionalMultiphaseStatistics.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsStatistics.hpp"
#include "SinglePhaseStatistics.hpp"
//#include "physicsSolvers/SolverStatistics.hpp"
//#include "physicsSolvers/SolverBase.hpp"

namespace geos
{
StatOutputController::StatOutputController( const string & name,
                                            Group * const parent ):
  TaskBase( name, parent ),
  m_timeHistory( nullptr )
{}

void StatOutputController::postInputInitialization()
{
  std::cout << " postInputInitialization " << std::endl;

  //TODO refactor hard coded outputs string (defined in problem manager)
  // TODO unicité
  OutputManager & outputManager = this->getGroupByPath< OutputManager >( "/Outputs" );
  m_timeHistory = &outputManager.registerGroup< TimeHistoryOutput >( "generatedCompFlowHistoryOutput" );

  EventManager & eventManager = this->getGroupByPath< EventManager >( "/Events" );
  m_periodicEvents{&eventManager.registerGroup< PeriodicEvent >( "generatedCompflowStats" ),
                   &eventManager.registerGroup< PeriodicEvent >( "generatedCompflowOutput" )};

  PeriodicEvent * periodicFlowStats = m_periodicEvents[0];
  periodicFlowStats->getWrapper< string >( viewKeyStruct::timeFrequencyString()).
    setApplyDefaultValue( "1e6" );
  periodicFlowStats->getWrapper< string >( viewKeyStruct::eventTargetString()).
    setApplyDefaultValue( "/Tasks/compflowStatistics" );
  PeriodicEvent * periodicFlowOutput = m_periodicEvents[1];
  periodicFlowStats->getWrapper< string >( viewKeyStruct::forceDtString()).
    setApplyDefaultValue( "4.32e6" );
  periodicFlowStats->getWrapper< string >( viewKeyStruct::eventTargetString()).
    setApplyDefaultValue( "/Outputs/compflowHistory" );

  TasksManager & taskManager = this->getGroupByPath< TasksManager >( "/Tasks" );
  RegionStatistics & regionStats = this->getGroupByPath( "/domain/MeshBodies/mesh/meshLevels/Level0/ElementRegions/elementRegionsGroup/aquiferTop/regionStatistics" );
  regionStats.forWrapper( [&]( WrapperBase const & wrapper )
  {
    string keyGroup = GEOS_FMT( "packCollection{}", wrapper.getName());

    PackCollection * packCollection = &taskManager.registerGroup< PackCollection >( keyGroup );
    packCollection->getWrapper< string >( viewKeyStruct::objectPathString()).
      setApplyDefaultValue( "/domain/MeshBodies/mesh/meshLevels/Level0/ElementRegions/elementRegionsGroup/aquiferTop/regionStatistics" );
    packCollection->getWrapper< string >( viewKeyStruct::fieldNameString()).
      setApplyDefaultValue( wrapper.getName());

    m_packCollections.push_back( packCollection );
  }
                          );
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
  createChild( SolidMechanicsStatistics::catalogName(), SolidMechanicsStatistics::catalogName() );
  createChild( SinglePhaseStatistics::catalogName(), SinglePhaseStatistics::catalogName() );
  //createChild( SolverBase::groupKeyStruct::solverStatisticsString(), SolverBase::groupKeyStruct::solverStatisticsString() );
}

bool StatOutputController::execute( real64 const time_n,
                                    real64 const dt,
                                    integer const GEOS_UNUSED_PARAM( cycleNumber ),
                                    integer const GEOS_UNUSED_PARAM( eventCounter ),
                                    real64 const GEOS_UNUSED_PARAM( eventProgress ),
                                    DomainPartition & domain )
{
  std::cout << " Execute " << std::endl;
  std::cout << domain.getName()<< time_n << dt <<std::endl;

  return false;
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        StatOutputController,
                        string const &, dataRepository::Group * const )

}
