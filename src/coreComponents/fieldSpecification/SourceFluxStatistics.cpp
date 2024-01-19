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

//TODO add a test for this component in SinglePhase and MultiPhase

#include "SourceFluxStatistics.hpp"

#include "SourceFluxBoundaryCondition.hpp"
#include "FieldSpecificationManager.hpp"

namespace geos
{
using namespace dataRepository;

SourceFluxStatistics::SourceFluxStatistics( const string & name,
                                            Group * const parent ):
  Base( name, parent )
{
  getWrapper< integer >( Group::viewKeyStruct::logLevelString ).
    appendDescription( GEOS_FMT( "\n- Log Level 1 outputs the sum of all {0}(s) produced rate & mass,\n"
                                 "- Log Level 2 outputs detailed values for each {0}.",
                                 SourceFluxBoundaryCondition::catalogName() ) );

  registerWrapper( viewKeyStruct::setNamesString(), &m_fluxNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( GEOS_FMT( "Name(s) array of the {0}(s) for which we want the statistics."
                              "Use \"{ all }\" to target all {0}.",
                              SourceFluxBoundaryCondition::catalogName() ) );
}

void SourceFluxStatistics::postProcessInput()
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  if( m_fluxNames.size() == 1 && m_fluxNames[0] == "all" )
  {
    m_fluxNames.clear();
    fsManager.forSubGroups< SourceFluxBoundaryCondition >( [&]( SourceFluxBoundaryCondition & sourceFlux )
    {
      m_set.push_back( &sourceFlux );
    } );
  }
  else
  {
    for( string const & fluxName : m_fluxNames )
    {
      //TODO : test this error
      GEOS_ERROR_IF( !fsManager.hasGroup< SourceFluxBoundaryCondition >( fluxName ),
                     GEOS_FMT( "{}: No {} named {} was found in {}.",
                               getDataContext(), SourceFluxBoundaryCondition::catalogName(),
                               fluxName, fs.getDataContext() ) );
    }
  }
  //TODO : remove that
  for( string const & flux : m_set )
  {
    GEOS_LOG( flux.getName());
  }
}

void SourceFluxStatistics::registerDataOnMesh( Group & meshBodies )
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
    mesh.getElemManager().forElementRegions( [&]( ElementRegionBase const & region )
    {

      WrapperBase & regionStats = region.registerWrapper< Stats >( getRegionStatsName() );
      regionStats.setRestartFlags( RestartFlags::NO_WRITE );

      region.excludeWrappersFromPacking( { getRegionStatsName() } );
    } );
  } );
}

void SourceFluxStatistics::writeStats( SourceFluxStatistics::Stats const & stats )

bool SourceFluxStatistics::execute( real64 const GEOS_UNUSED_PARAM( time_n ),
                                    real64 const GEOS_UNUSED_PARAM( dt ),
                                    integer const GEOS_UNUSED_PARAM( cycleNumber ),
                                    integer const GEOS_UNUSED_PARAM( eventCounter ),
                                    real64 const GEOS_UNUSED_PARAM( eventProgress ),
                                    DomainPartition & domain )
{
  if( getLogLevel() >= 1 )
  {
    // m_solver->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
    //                                                                         MeshLevel & mesh,
    //                                                                         arrayView1d< string const > const & regionNames )
    {
      //TODO : get back all regions info on sourceflux region stat wrappers
      GEOS_LOG( GEOS_FMT( "{} {}: SourceFluxStatistics::execute::allFluxStats not yet implemented",
                          catalogName(), getName() ) );
      if( getLogLevel() >= 2 )
      {
        GEOS_LOG( GEOS_FMT( "{} {}: SourceFluxStatistics::execute::perFluxStats not yet implemented",
                            catalogName(), getName() ) );
      }
    } //);
  }
  return false;
}

} /* namespace geos */
