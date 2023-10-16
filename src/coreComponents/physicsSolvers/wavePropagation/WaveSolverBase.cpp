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
#include "WaveSolverUtils.hpp"

#include <limits>

namespace geos
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

  registerWrapper( viewKeyStruct::receiverCoordinatesString(), &m_receiverCoordinates ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Coordinates (x,y,z) of the receivers" );

  registerWrapper( viewKeyStruct::sourceValueString(), &m_sourceValue ).
    setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setSizedFromParent( 0 ).
    setDescription( "Source Value of the sources" );

  registerWrapper( viewKeyStruct::timeSourceFrequencyString(), &m_timeSourceFrequency ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Central frequency for the time source" );

  registerWrapper( viewKeyStruct::timeSourceDelayString(), &m_timeSourceDelay ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( -1 ).
    setDescription( "Source time delay (1 / f0 by default)" );

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

  registerWrapper( viewKeyStruct::enableLifoString(), &m_enableLifo ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Set to 1 to enable LIFO storage feature" );

  registerWrapper( viewKeyStruct::lifoSizeString(), &m_lifoSize ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( std::numeric_limits< int >::max() ).
    setDescription( "Set the capacity of the lifo storage (should be the total number of buffers to store in the LIFO)" );

  registerWrapper( viewKeyStruct::lifoOnDeviceString(), &m_lifoOnDevice ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( -80 ).
    setDescription( "Set the capacity of the lifo device storage (if negative, opposite of percentage of remaining memory)" );

  registerWrapper( viewKeyStruct::lifoOnHostString(), &m_lifoOnHost ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( -80 ).
    setDescription( "Set the capacity of the lifo host storage (if negative, opposite of percentage of remaining memory)" );

  registerWrapper( viewKeyStruct::usePMLString(), &m_usePML ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag to apply PML" );

  registerWrapper( viewKeyStruct::useDASString(), &m_useDAS ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag to indicate if DAS data will be modeled, and which DAS type to use: 1 for strain integration, 2 for displacement difference" );

  registerWrapper( viewKeyStruct::linearDASSamplesString(), &m_linearDASSamples ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 5 ).
    setDescription( "Number of sample points to be used for strain integration when integrating the strain for the DAS signal" );

  registerWrapper( viewKeyStruct::linearDASGeometryString(), &m_linearDASGeometry ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Geometry parameters for a linear DAS fiber (dip, azimuth, gauge length)" );

  registerWrapper( viewKeyStruct::sourceNodeIdsString(), &m_sourceNodeIds ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Indices of the nodes (in the right order) for each source point" );

  registerWrapper( viewKeyStruct::sourceConstantsString(), &m_sourceConstants ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Constant part of the source for the nodes listed in m_sourceNodeIds" );

  registerWrapper( viewKeyStruct::sourceIsAccessibleString(), &m_sourceIsAccessible ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Flag that indicates whether the source is local to this MPI rank" );

  registerWrapper( viewKeyStruct::receiverNodeIdsString(), &m_receiverNodeIds ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Indices of the nodes (in the right order) for each receiver point" );

  registerWrapper( viewKeyStruct::sourceConstantsString(), &m_sourceConstants ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Constant part of the receiver for the nodes listed in m_receiverNodeIds" );

  registerWrapper( viewKeyStruct::receiverIsLocalString(), &m_receiverIsLocal ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Flag that indicates whether the receiver is local to this MPI rank" );


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

void WaveSolverBase::registerDataOnMesh( Group & meshBodies )
{
  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    nodeManager.registerField< fields::referencePosition32 >( this->getName() );
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

    nodeManager.getField< fields::referencePosition32 >().resizeDimension< 1 >( X.size( 1 ) );
    arrayView2d< wsCoordType, nodes::REFERENCE_POSITION_USD > const X32  = nodeManager.getField< fields::referencePosition32 >();
    for( int i = 0; i < X.size( 0 ); i++ )
    {
      for( int j = 0; j < X.size( 1 ); j++ )
      {
        X32[i][j] = X[i][j];
      }
    }
  } );
}

void WaveSolverBase::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  localIndex const numNodesPerElem = WaveSolverBase::getNumNodesPerElem();

  localIndex const numSourcesGlobal = m_sourceCoordinates.size( 0 );
  m_sourceNodeIds.resize( numSourcesGlobal, numNodesPerElem );
  m_sourceConstants.resize( numSourcesGlobal, numNodesPerElem );
  m_sourceIsAccessible.resize( numSourcesGlobal );

  localIndex const numReceiversGlobal = m_receiverCoordinates.size( 0 );
  m_receiverNodeIds.resize( numReceiversGlobal, numNodesPerElem );
  m_receiverConstants.resize( numReceiversGlobal, numNodesPerElem );
  m_receiverIsLocal.resize( numReceiversGlobal );

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
  GEOS_THROW_IF( counter > 1,
                 getDataContext() << ": One single PML field specification is allowed",
                 InputError );

  m_usePML = counter;

  if( m_useDAS == 0 && m_linearDASGeometry.size( 0 ) > 0 )
  {
    m_useDAS = 1;
    m_linearDASSamples = 5;
  }

  if( m_useDAS == 2 )
  {
    m_linearDASSamples = 2;
  }

  if( m_useDAS > 0 )
  {
    GEOS_LOG_LEVEL_RANK_0( 1, "Modeling linear DAS data is activated" );
    GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "Linear DAS formulation: {}", m_useDAS == 1 ? "strain integration" : "displacement difference" ) );

    GEOS_ERROR_IF( m_linearDASGeometry.size( 1 ) != 3,
                   getWrapperDataContext( viewKeyStruct::linearDASGeometryString() ) <<
                   ": Invalid number of geometry parameters for the linear DAS fiber. Three parameters are required: dip, azimuth, gauge length" );

    GEOS_ERROR_IF( m_linearDASGeometry.size( 0 ) != m_receiverCoordinates.size( 0 ),
                   getWrapperDataContext( viewKeyStruct::linearDASGeometryString() ) <<
                   ": Invalid number of geometry parameters instances for the linear DAS fiber. It should match the number of receivers." );

    /// initialize DAS geometry by adding receiver locations along the DAS, as needed
    initializeDAS();
  }

  GEOS_THROW_IF( m_sourceCoordinates.size( 1 ) != 3,
                 getWrapperDataContext( viewKeyStruct::sourceCoordinatesString() ) <<
                 ": Invalid number of physical coordinates for the sources",
                 InputError );

  GEOS_THROW_IF( m_receiverCoordinates.size( 1 ) != 3,
                 getWrapperDataContext( viewKeyStruct::receiverCoordinatesString() ) <<
                 ": Invalid number of physical coordinates for the receivers",
                 InputError );

  EventManager const & event = this->getGroupByPath< EventManager >( "/Problem/Events" );
  real64 const & maxTime = event.getReference< real64 >( EventManager::viewKeyStruct::maxTimeString() );
  real64 const & minTime = event.getReference< real64 >( EventManager::viewKeyStruct::minTimeString() );
  real64 dt = 0;
  for( localIndex numSubEvent = 0; numSubEvent < event.numSubGroups(); ++numSubEvent )
  {
    EventBase const * subEvent = static_cast< EventBase const * >( event.getSubGroups()[numSubEvent] );
    if( subEvent->getEventName() == "/Solvers/" + this->getName() )
    {
      dt = subEvent->getReference< real64 >( EventBase::viewKeyStruct::forceDtString() );
    }
  }

  GEOS_THROW_IF( dt < epsilonLoc*maxTime, getDataContext() << ": Value for dt: " << dt <<" is smaller than local threshold: " << epsilonLoc, std::runtime_error );

  if( m_dtSeismoTrace > 0 )
  {
    m_nsamplesSeismoTrace = int( maxTime / m_dtSeismoTrace) + 1;
  }
  else
  {
    m_nsamplesSeismoTrace = 0;
  }
  localIndex const nsamples = int( (maxTime-minTime) /dt) + 1;

  localIndex const numSourcesGlobal = m_sourceCoordinates.size( 0 );
  m_sourceValue.resize( nsamples, numSourcesGlobal );

}

void WaveSolverBase::initializeDAS()
{
  /// increase the number of receivers and modify their coordinates
  /// so to have two correct number of receivers along each DAS channel
  localIndex const numReceiversGlobal = m_receiverCoordinates.size( 0 );
  m_receiverCoordinates.resize( m_linearDASSamples*numReceiversGlobal, 3 );

  arrayView2d< real64 > const receiverCoordinates = m_receiverCoordinates.toView();
  arrayView2d< real64 const > const linearDASGeometry = m_linearDASGeometry.toViewConst();

  array1d< real64 > const samplePointLocationsA( m_linearDASSamples );
  arrayView1d< real64 > const samplePointLocations = samplePointLocationsA.toView();
  if( m_linearDASSamples == 1 )
  {
    samplePointLocations[ 0 ] = 0;
  }
  else
  {
    for( integer i = 0; i < m_linearDASSamples; ++i )
    {
      samplePointLocations[ i ] = -0.5 + (real64) i / ( m_linearDASSamples - 1 );
    }
  }

  for( localIndex ircv = 0; ircv < numReceiversGlobal; ++ircv )
  {
    R1Tensor receiverCenter = { receiverCoordinates[ ircv ][ 0 ], receiverCoordinates[ ircv ][ 1 ], receiverCoordinates[ ircv ][ 2 ] };
    for( integer i = 0; i < m_linearDASSamples; ++i )
    {
      /// updated xyz of sample receiver along DAS
      R1Tensor dasVector = WaveSolverUtils::computeDASVector( linearDASGeometry[ ircv ][ 0 ], linearDASGeometry[ ircv ][ 1 ] );
      receiverCoordinates[ i*numReceiversGlobal + ircv ][ 0 ] = receiverCenter[ 0 ] + dasVector[ 0 ] * linearDASGeometry[ ircv ][ 2 ] * samplePointLocations[ i ];
      receiverCoordinates[ i*numReceiversGlobal + ircv ][ 1 ] = receiverCenter[ 1 ] + dasVector[ 1 ] * linearDASGeometry[ ircv ][ 2 ] * samplePointLocations[ i ];
      receiverCoordinates[ i*numReceiversGlobal + ircv ][ 2 ] = receiverCenter[ 2 ] + dasVector[ 2 ] * linearDASGeometry[ ircv ][ 2 ] * samplePointLocations[ i ];
    }
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

localIndex WaveSolverBase::getNumNodesPerElem()
{
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();

  FiniteElementDiscretizationManager const &
  feDiscretizationManager = numericalMethodManager.getFiniteElementDiscretizationManager();

  FiniteElementDiscretization const * const
  feDiscretization = feDiscretizationManager.getGroupPointer< FiniteElementDiscretization >( m_discretizationName );
  GEOS_THROW_IF( feDiscretization == nullptr,
                 getDataContext() << ": FE discretization not found: " << m_discretizationName,
                 InputError );

  localIndex numNodesPerElem = 0;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                  [&]( string const &,
                                       MeshLevel const & mesh,
                                       arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();
    elemManager.forElementRegions( regionNames,
                                   [&] ( localIndex const,
                                         ElementRegionBase const & elemRegion )
    {
      elemRegion.forElementSubRegions( [&]( ElementSubRegionBase const & elementSubRegion )
      {
        finiteElement::FiniteElementBase const &
        fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
        localIndex const numSupportPoints = fe.getNumSupportPoints();
        if( numSupportPoints > numNodesPerElem )
        {
          numNodesPerElem = numSupportPoints;
        }
      } );
    } );


  } );
  return numNodesPerElem;

}

bool WaveSolverBase::directoryExists( std::string const & directoryName )
{
  struct stat buffer;
  return stat( directoryName.c_str(), &buffer ) == 0;
}


} /* namespace geos */
