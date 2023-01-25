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
 * @file WaveSolverBase.cpp
 */

#include "WaveSolverBase.hpp"

#include "dataRepository/KeyNames.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/PerfectlyMatchedLayer.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

namespace geosx
{

using namespace dataRepository;

WaveSolverBase::WaveSolverBase( const std::string & name,
                                Group * const parent ):
  SolverBase( name,
              parent )
{

  registerWrapper( viewKeyStruct::sourceCoordinatesString(), &m_sourceCoordinates ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Coordinates (x,y,z) of the sources" );

  registerWrapper( viewKeyStruct::sourceValueString(), &m_sourceValue ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setSizedFromParent( 0 ).
    setDescription( "Source Value of the sources" );

  registerWrapper( viewKeyStruct::timeSourceFrequencyString(), &m_timeSourceFrequency ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Central frequency for the time source" );

  registerWrapper( viewKeyStruct::receiverCoordinatesString(), &m_receiverCoordinates ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Coordinates (x,y,z) of the receivers" );

  registerWrapper( viewKeyStruct::rickerOrderString(), &m_rickerOrder ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 2 ).
    setDescription( "Flag that indicates the order of the Ricker to be used o, 1 or 2. Order 2 by default" );

  registerWrapper( viewKeyStruct::outputSeismoTraceString(), &m_outputSeismoTrace ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag that indicates if we write the seismo trace in a file .txt, 0 no output, 1 otherwise" );

  registerWrapper( viewKeyStruct::dtSeismoTraceString(), &m_dtSeismoTrace ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Time step for output pressure at receivers" );

  registerWrapper( viewKeyStruct::indexSeismoTraceString(), &m_indexSeismoTrace ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( 0 ).
    setDescription( "Count for output pressure at receivers" );

  registerWrapper( viewKeyStruct::forwardString(), &m_forward ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1 ).
    setDescription( "Set to 1 to compute forward propagation" );

  registerWrapper( viewKeyStruct::saveFieldsString(), &m_saveFields ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Set to 1 to save fields during forward and restore them during backward" );


  registerWrapper( viewKeyStruct::shotIndexString(), &m_shotIndex ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Set the current shot for temporary files" );

  registerWrapper( viewKeyStruct::usePMLString(), &m_usePML ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag to apply PML" );

  registerWrapper( viewKeyStruct::useDASString(), &m_useDAS ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag to indicate if DAS type of data will be modeled" );

  registerWrapper( viewKeyStruct::linearDASGeometryString(), &m_linearDASGeometry ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Geometry parameters for a linear DAS fiber (dip, azimuth, gauge length)" );

}

WaveSolverBase::~WaveSolverBase()
{
  // TODO Auto-generated destructor stub
}

void WaveSolverBase::reinit()
{
  initializePreSubGroups();
  postProcessInput();
  initializePostInitialConditionsPreSubGroups();
}

void WaveSolverBase::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();
}

void WaveSolverBase::postProcessInput()
{
  SolverBase::postProcessInput();

  /// set flag PML to one if a PML field is specified in the xml
  /// if counter>1, an error will be thrown as one single PML field is allowed
  integer counter = 0;
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  fsManager.forSubGroups< PerfectlyMatchedLayer >( [&] ( PerfectlyMatchedLayer const & )
  {
    counter++;
  } );
  GEOSX_THROW_IF( counter > 1,
                  "One single PML field specification is allowed",
                  InputError );

  m_usePML = counter;

  if( m_linearDASGeometry.size( 1 ) > 0 )
  {
    m_useDAS = 1;
  }

  if( m_useDAS )
  {
    GEOSX_LOG_LEVEL_RANK_0( 1, "Modeling linear DAS data is activated" );

    GEOSX_ERROR_IF( m_linearDASGeometry.size( 1 ) != 3,
                    "Invalid number of geometry parameters for the linear DAS fiber. Three parameters are required: dip, azimuth, gauge length" );

    GEOSX_ERROR_IF( m_linearDASGeometry.size( 0 ) != m_receiverCoordinates.size( 0 ),
                    "Invalid number of geometry parameters instances for the linear DAS fiber. It should match the number of receivers." );

    /// initialize DAS geometry
    initializeDAS();
  }
}

void WaveSolverBase::initializeDAS()
{
  /// double the number of receivers and modify their coordinates
  /// so to have two receivers on each side of each DAS channel
  localIndex const numReceiversGlobal = m_receiverCoordinates.size( 0 );
  m_receiverCoordinates.resize( 2*numReceiversGlobal, 3 );

  arrayView2d< real64 > const receiverCoordinates = m_receiverCoordinates.toView();
  arrayView2d< real64 const > const linearDASGeometry = m_linearDASGeometry.toViewConst();

  for( localIndex ircv = 0; ircv < numReceiversGlobal; ++ircv )
  {
    /// updated xyz of receivers on the far end of a DAS channel
    receiverCoordinates[numReceiversGlobal+ircv][0] = receiverCoordinates[ircv][0]
                                                      + cos( linearDASGeometry[ircv][0] ) * cos( linearDASGeometry[ircv][1] ) * linearDASGeometry[ircv][2] / 2.0;
    receiverCoordinates[numReceiversGlobal+ircv][1] = receiverCoordinates[ircv][1]
                                                      + cos( linearDASGeometry[ircv][0] ) * sin( linearDASGeometry[ircv][1] ) * linearDASGeometry[ircv][2] / 2.0;
    receiverCoordinates[numReceiversGlobal+ircv][2] = receiverCoordinates[ircv][2]
                                                      + sin( linearDASGeometry[ircv][0] ) * linearDASGeometry[ircv][2] / 2.0;

    /// updated xyz of receivers on the near end of a DAS channel
    receiverCoordinates[ircv][0] = receiverCoordinates[ircv][0]
                                   - cos( linearDASGeometry[ircv][0] ) * cos( linearDASGeometry[ircv][1] ) * linearDASGeometry[ircv][2] / 2.0;
    receiverCoordinates[ircv][1] = receiverCoordinates[ircv][1]
                                   - cos( linearDASGeometry[ircv][0] ) * sin( linearDASGeometry[ircv][1] ) * linearDASGeometry[ircv][2] / 2.0;
    receiverCoordinates[ircv][2] = receiverCoordinates[ircv][2]
                                   - sin( linearDASGeometry[ircv][0] ) * linearDASGeometry[ircv][2] / 2.0;
  }
}

real64 WaveSolverBase::solverStep( real64 const & time_n,
                                   real64 const & dt,
                                   integer const cycleNumber,
                                   DomainPartition & domain )
{
  return explicitStep( time_n, dt, cycleNumber, domain );
}

real64 WaveSolverBase::explicitStep( real64 const & time_n,
                                     real64 const & dt,
                                     integer const cycleNumber,
                                     DomainPartition & domain )
{
  if( m_forward )
  {
    return explicitStepForward( time_n, dt, cycleNumber, domain, m_saveFields );
  }
  else
  {
    return explicitStepBackward( time_n, dt, cycleNumber, domain, m_saveFields );
  }
}
} /* namespace geosx */
