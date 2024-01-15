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
 * @file ElasticWaveEquationSEM.cpp
 */

#include "ElasticWaveEquationSEM.hpp"
#include "ElasticWaveEquationSEMKernel.hpp"

#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/ElementType.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "WaveSolverUtils.hpp"

namespace geos
{

using namespace dataRepository;

ElasticWaveEquationSEM::ElasticWaveEquationSEM( const std::string & name,
                                                Group * const parent ):

  WaveSolverBase( name,
                  parent )
{
  registerWrapper( viewKeyStruct::sourceConstantsString(), &m_sourceConstantsx ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Constant part of the source for the nodes listed in m_sourceNodeIds in x-direction" );

  registerWrapper( viewKeyStruct::sourceConstantsString(), &m_sourceConstantsy ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Constant part of the source for the nodes listed in m_sourceNodeIds in y-direction" );

  registerWrapper( viewKeyStruct::sourceConstantsString(), &m_sourceConstantsz ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Constant part of the source for the nodes listed in m_sourceNodeIds in z-direction" );

  registerWrapper( viewKeyStruct::displacementXNp1AtReceiversString(), &m_displacementXNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Displacement value at each receiver for each timestep (x-component)" );

  registerWrapper( viewKeyStruct::displacementYNp1AtReceiversString(), &m_displacementYNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Displacement value at each receiver for each timestep (y-component)" );

  registerWrapper( viewKeyStruct::displacementZNp1AtReceiversString(), &m_displacementZNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Displacement value at each receiver for each timestep (z-component)" );

  registerWrapper( viewKeyStruct::sourceForceString(), &m_sourceForce ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setApplyDefaultValue( { 0.0, 0.0, 0.0 } ).
    setDescription( "Force of the source: 3 real values for a vector source, and 6 real values for a tensor source (in Voigt notation)."
                    "The default value is { 0, 0, 0 } (no net force)." );

  registerWrapper( viewKeyStruct::sourceMomentString(), &m_sourceMoment ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setApplyDefaultValue( { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 } ).
    setDescription( "Moment of the source: 6 real values describing a symmetric tensor in Voigt notation."
                    "The default value is { 1, 1, 1, 0, 0, 0 } (diagonal moment, corresponding to a pure explosion)." );
}

ElasticWaveEquationSEM::~ElasticWaveEquationSEM()
{
  // TODO Auto-generated destructor stub
}

void ElasticWaveEquationSEM::initializePreSubGroups()
{

  WaveSolverBase::initializePreSubGroups();

  localIndex const numNodesPerElem = getNumNodesPerElem();

  localIndex const numSourcesGlobal = m_sourceCoordinates.size( 0 );
  m_sourceConstantsx.resize( numSourcesGlobal, numNodesPerElem );
  m_sourceConstantsy.resize( numSourcesGlobal, numNodesPerElem );
  m_sourceConstantsz.resize( numSourcesGlobal, numNodesPerElem );

  localIndex const numReceiversGlobal = m_receiverCoordinates.size( 0 );
  m_receiverConstants.resize( numReceiversGlobal, numNodesPerElem );

}


void ElasticWaveEquationSEM::registerDataOnMesh( Group & meshBodies )
{
  WaveSolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    nodeManager.registerField< fields::Displacementx_nm1,
                               fields::Displacementy_nm1,
                               fields::Displacementz_nm1,
                               fields::Displacementx_n,
                               fields::Displacementy_n,
                               fields::Displacementz_n,
                               fields::Displacementx_np1,
                               fields::Displacementy_np1,
                               fields::Displacementz_np1,
                               fields::ForcingRHSx,
                               fields::ForcingRHSy,
                               fields::ForcingRHSz,
                               fields::ElasticMassVector,
                               fields::DampingVectorx,
                               fields::DampingVectory,
                               fields::DampingVectorz,
                               fields::StiffnessVectorx,
                               fields::StiffnessVectory,
                               fields::StiffnessVectorz,
                               fields::ElasticFreeSurfaceNodeIndicator >( getName() );

    FaceManager & faceManager = mesh.getFaceManager();
    faceManager.registerField< fields::ElasticFreeSurfaceFaceIndicator >( getName() );

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerField< fields::ElasticVelocityVp >( getName() );
      subRegion.registerField< fields::ElasticVelocityVs >( getName() );
      subRegion.registerField< fields::ElasticDensity >( getName() );
    } );

  } );
}



void ElasticWaveEquationSEM::postProcessInput()
{
  WaveSolverBase::postProcessInput();

  EventManager const & event = getGroupByPath< EventManager >( "/Problem/Events" );
  real64 const & maxTime = event.getReference< real64 >( EventManager::viewKeyStruct::maxTimeString() );
  real64 dt = 0;
  for( localIndex numSubEvent = 0; numSubEvent < event.numSubGroups(); ++numSubEvent )
  {
    EventBase const * subEvent = static_cast< EventBase const * >( event.getSubGroups()[numSubEvent] );
    if( subEvent->getEventName() == "/Solvers/" + getName() )
    {
      dt = subEvent->getReference< real64 >( EventBase::viewKeyStruct::forceDtString() );
    }
  }

  GEOS_THROW_IF( dt < epsilonLoc*maxTime, getDataContext() << ": Value for dt: " << dt <<" is smaller than local threshold: " << epsilonLoc, std::runtime_error );

  if( m_dtSeismoTrace > 0 )
  {
    m_nsamplesSeismoTrace = int( maxTime / m_dtSeismoTrace ) + 1;
  }
  else
  {
    m_nsamplesSeismoTrace = 0;
  }
  localIndex const nsamples = int( maxTime / dt ) + 1;

  localIndex const numReceiversGlobal = m_receiverCoordinates.size( 0 );
  localIndex const numSourcesGlobal = m_sourceCoordinates.size( 0 );
  m_sourceIsAccessible.resize( numSourcesGlobal );

  m_displacementXNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal + 1 );
  m_displacementYNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal + 1 );
  m_displacementZNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal + 1 );
  m_sourceValue.resize( nsamples, numSourcesGlobal );

  /// The receivers are initialized to zero.
  /// This is essential for DAS modeling as MPI_Allreduce is called when computing DAS data
  m_displacementXNp1AtReceivers.zero();
  m_displacementYNp1AtReceivers.zero();
  m_displacementZNp1AtReceivers.zero();
}


void ElasticWaveEquationSEM::precomputeSourceAndReceiverTerm( MeshLevel & mesh, arrayView1d< string const > const & regionNames )
{
  NodeManager const & nodeManager = mesh.getNodeManager();
  FaceManager const & faceManager = mesh.getFaceManager();

  arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const X =
    nodeManager.getField< fields::referencePosition32 >().toViewConst();
  arrayView2d< real64 const > const faceNormal  = faceManager.faceNormal();
  arrayView2d< real64 const > const faceCenter  = faceManager.faceCenter();

  arrayView2d< real64 const > const sourceCoordinates = m_sourceCoordinates.toViewConst();
  arrayView2d< localIndex > const sourceNodeIds = m_sourceNodeIds.toView();
  arrayView2d< real64 > const sourceConstantsx = m_sourceConstantsx.toView();
  arrayView2d< real64 > const sourceConstantsy = m_sourceConstantsy.toView();
  arrayView2d< real64 > const sourceConstantsz = m_sourceConstantsz.toView();
  arrayView1d< localIndex > const sourceIsAccessible = m_sourceIsAccessible.toView();

  sourceNodeIds.setValues< EXEC_POLICY >( -1 );
  sourceConstantsx.setValues< EXEC_POLICY >( 0 );
  sourceConstantsy.setValues< EXEC_POLICY >( 0 );
  sourceConstantsz.setValues< EXEC_POLICY >( 0 );
  sourceIsAccessible.zero();

  arrayView2d< real64 const > const receiverCoordinates = m_receiverCoordinates.toViewConst();
  arrayView2d< localIndex > const receiverNodeIds = m_receiverNodeIds.toView();
  arrayView2d< real64 > const receiverConstants = m_receiverConstants.toView();
  arrayView1d< localIndex > const receiverIsLocal = m_receiverIsLocal.toView();
  receiverNodeIds.setValues< EXEC_POLICY >( -1 );
  receiverConstants.setValues< EXEC_POLICY >( 0 );
  receiverIsLocal.zero();

  arrayView2d< real32 > const sourceValue = m_sourceValue.toView();

  real64 dt = 0;
  EventManager const & event = getGroupByPath< EventManager >( "/Problem/Events" );
  for( localIndex numSubEvent = 0; numSubEvent < event.numSubGroups(); ++numSubEvent )
  {
    EventBase const * subEvent = static_cast< EventBase const * >( event.getSubGroups()[numSubEvent] );
    if( subEvent->getEventName() == "/Solvers/" + getName() )
    {
      dt = subEvent->getReference< real64 >( EventBase::viewKeyStruct::forceDtString() );
    }
  }

  mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                        CellElementSubRegion & elementSubRegion )
  {

    GEOS_THROW_IF( elementSubRegion.getElementType() != ElementType::Hexahedron,
                   getDataContext() << ": Invalid type of element, the elastic solver is designed for hexahedral meshes only (C3D8) ",
                   InputError );

    arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();
    arrayView2d< real64 const > const elemCenter = elementSubRegion.getElementCenter();
    arrayView1d< integer const > const elemGhostRank = elementSubRegion.ghostRank();

    finiteElement::FiniteElementBase const &
    fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
    finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
    {
      using FE_TYPE = TYPEOFREF( finiteElement );

      localIndex const numFacesPerElem = elementSubRegion.numFacesPerElement();

      elasticWaveEquationSEMKernels::
        PrecomputeSourceAndReceiverKernel::
        launch< EXEC_POLICY, FE_TYPE >
        ( elementSubRegion.size(),
        numFacesPerElem,
        X,
        elemGhostRank,
        elemsToNodes,
        elemsToFaces,
        elemCenter,
        faceNormal,
        faceCenter,
        sourceCoordinates,
        sourceIsAccessible,
        sourceNodeIds,
        sourceConstantsx,
        sourceConstantsy,
        sourceConstantsz,
        receiverCoordinates,
        receiverIsLocal,
        receiverNodeIds,
        receiverConstants,
        sourceValue,
        dt,
        m_timeSourceFrequency,
        m_timeSourceDelay,
        m_rickerOrder,
        m_useDAS,
        m_linearDASSamples,
        m_linearDASGeometry.toViewConst(),
        m_sourceForce,
        m_sourceMoment );
    } );
  } );
}

void ElasticWaveEquationSEM::computeDAS( arrayView2d< real32 > const xCompRcv,
                                         arrayView2d< real32 > const yCompRcv,
                                         arrayView2d< real32 > const zCompRcv )
{
  arrayView2d< real64 const > const linearDASGeometry = m_linearDASGeometry.toViewConst();
  arrayView1d< localIndex const > const receiverIsLocal = m_receiverIsLocal.toViewConst();

  localIndex const numReceiversGlobal = m_receiverIsLocal.size();
  localIndex const numDAS = linearDASGeometry.size( 0 );
  localIndex const nsamplesSeismoTrace = m_nsamplesSeismoTrace;
  localIndex const nlinearDASSamples = m_linearDASSamples;
  localIndex const useDAS = m_useDAS;
   

  array1d< real64 > const samplePointIntegrationConstantsA( m_linearDASSamples );
  /// for displacement difference DAS (m_useDAS==2), take the difference of the pair of geophones
  if( m_useDAS == 2 )
  {
    samplePointIntegrationConstantsA[ 0 ] = -1.0;
    samplePointIntegrationConstantsA[ 1 ] = 1.0;
  }
  /// for strain integration DAS (m_useDAS==1), take the average of strains to average strain data
  else if( m_linearDASSamples == 1 )
  {
    samplePointIntegrationConstantsA[ 0 ] = 1.0;
  }
  else
  {
    for( integer i = 0; i < m_linearDASSamples; i++ )
    {
      samplePointIntegrationConstantsA[ i ] = 1.0 / m_linearDASSamples;
    }
  }
  arrayView1d< real64 const > const samplePointIntegrationConstants = samplePointIntegrationConstantsA.toViewConst();


  if( m_nsamplesSeismoTrace > 0 )
  {
    // set time column to zero on two components, to avoid summing them in the following MPI reduce
    if( MpiWrapper::commRank( MPI_COMM_GEOSX ) != 0)
    {
      for( localIndex iTimeSample = 0; iTimeSample < nsamplesSeismoTrace; ++iTimeSample )
      {
        xCompRcv( iTimeSample, numReceiversGlobal ) = 0;
        yCompRcv( iTimeSample, numReceiversGlobal ) = 0;
        zCompRcv( iTimeSample, numReceiversGlobal ) = 0;
      }
    }
    /// synchronize receivers across MPI ranks
    MpiWrapper::allReduce( xCompRcv.data(),
                           xCompRcv.data(),
                           (numReceiversGlobal + 1)*nsamplesSeismoTrace,
                           MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
                           MPI_COMM_GEOSX );

    MpiWrapper::allReduce( yCompRcv.data(),
                           yCompRcv.data(),
                           (numReceiversGlobal + 1)*nsamplesSeismoTrace,
                           MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
                           MPI_COMM_GEOSX );

    MpiWrapper::allReduce( zCompRcv.data(),
                           zCompRcv.data(),
                           (numReceiversGlobal + 1)*nsamplesSeismoTrace,
                           MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
                           MPI_COMM_GEOSX );
    localIndex const myrank = MpiWrapper::commRank( MPI_COMM_GEOSX );
    forAll< EXEC_POLICY >( numDAS, [=] GEOS_HOST_DEVICE ( localIndex const ircv )
    {
      if( receiverIsLocal[ ircv ] == 1 )
      {
        R1Tensor dasVector = WaveSolverUtils::computeDASVector( linearDASGeometry[ ircv ][ 0 ], linearDASGeometry[ ircv ][ 1 ] );

        for( localIndex iTimeSample = 0; iTimeSample < nsamplesSeismoTrace; ++iTimeSample )
        {
          real64 das = 0;
          /// convert point data to DAS data
          for( localIndex iSample = 0; iSample < nlinearDASSamples; ++iSample )
          {
            // store strain data in the z-component of the receiver (copied to x after resize)
            das += ( dasVector[ 0 ] * xCompRcv( iTimeSample, iSample*numDAS+ircv )
                     + dasVector[ 1 ] * yCompRcv( iTimeSample, iSample*numDAS+ircv )
                     + dasVector[ 2 ] * zCompRcv( iTimeSample, iSample*numDAS+ircv ) ) * samplePointIntegrationConstants[ iSample ];
          }
          if( useDAS == 2 )
          {
            das /= linearDASGeometry[ ircv ][ 2 ];
          }
          zCompRcv( iTimeSample, ircv ) = das;
          // also copy the last receiver, which contains the seismogram time
          if( ircv == 0 )
          {
            zCompRcv( iTimeSample, numDAS ) = zCompRcv( iTimeSample, numReceiversGlobal );
          }
        }
      }
    } );
  }

  /// resize the receiver arrays by dropping the extra data to avoid confusion
  /// the remaining x-component contains DAS data, the other components are set to zero
  m_displacementXNp1AtReceivers.resize( m_nsamplesSeismoTrace, numDAS + 1 );
  arrayView2d< real32 > const dasReceiver = m_displacementXNp1AtReceivers.toView();
  forAll< EXEC_POLICY >( numDAS, [=] GEOS_HOST_DEVICE ( localIndex const ircv )
  {
    if( receiverIsLocal[ ircv ] == 1 )
    {
      /// convert dipole data (pairs of geophones) to average strain data and
      for( localIndex iTimeSample = 0; iTimeSample < nsamplesSeismoTrace; ++iTimeSample )
      {
        // store strain data in the z-component of the receiver (copied to x after resize)
        dasReceiver( iTimeSample, ircv ) = zCompRcv( iTimeSample, ircv );
        if( ircv == 0 )
        {
          dasReceiver( iTimeSample, numDAS ) = zCompRcv( iTimeSample, numDAS );
        }
      }
    }
  } );
  /// set the y and z components to zero to avoid any confusion
  m_displacementYNp1AtReceivers.resize( m_nsamplesSeismoTrace, numDAS + 1 );
  m_displacementYNp1AtReceivers.zero();
  m_displacementZNp1AtReceivers.resize( m_nsamplesSeismoTrace, numDAS + 1 );
  m_displacementZNp1AtReceivers.zero();
}

void ElasticWaveEquationSEM::addSourceToRightHandSide( integer const & cycleNumber,
                                                       arrayView1d< real32 > const rhsx,
                                                       arrayView1d< real32 > const rhsy,
                                                       arrayView1d< real32 > const rhsz )
{
  arrayView2d< localIndex const > const sourceNodeIds = m_sourceNodeIds.toViewConst();
  arrayView2d< real64 const > const sourceConstantsx   = m_sourceConstantsx.toViewConst();
  arrayView2d< real64 const > const sourceConstantsy   = m_sourceConstantsy.toViewConst();
  arrayView2d< real64 const > const sourceConstantsz   = m_sourceConstantsz.toViewConst();

  arrayView1d< localIndex const > const sourceIsAccessible = m_sourceIsAccessible.toViewConst();
  arrayView2d< real32 const > const sourceValue   = m_sourceValue.toViewConst();

  GEOS_THROW_IF( cycleNumber > sourceValue.size( 0 ), getDataContext() << ": Too many steps compared to array size", std::runtime_error );
  forAll< EXEC_POLICY >( m_sourceConstantsx.size( 0 ), [=] GEOS_HOST_DEVICE ( localIndex const isrc )
  {
    if( sourceIsAccessible[isrc] == 1 )
    {
      for( localIndex inode = 0; inode < sourceConstantsx.size( 1 ); ++inode )
      {
        real32 const localIncrementx = sourceConstantsx[isrc][inode] * sourceValue[cycleNumber][isrc];
        RAJA::atomicAdd< ATOMIC_POLICY >( &rhsx[sourceNodeIds[isrc][inode]], localIncrementx );
        real32 const localIncrementy = sourceConstantsy[isrc][inode] * sourceValue[cycleNumber][isrc];
        RAJA::atomicAdd< ATOMIC_POLICY >( &rhsy[sourceNodeIds[isrc][inode]], localIncrementy );
        real32 const localIncrementz = sourceConstantsz[isrc][inode] * sourceValue[cycleNumber][isrc];
        RAJA::atomicAdd< ATOMIC_POLICY >( &rhsz[sourceNodeIds[isrc][inode]], localIncrementz );
      }
    }
  } );
}

void ElasticWaveEquationSEM::initializePostInitialConditionsPreSubGroups()
{

  WaveSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = getGroupByPath< DomainPartition >( "/Problem/domain" );

  applyFreeSurfaceBC( 0.0, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    precomputeSourceAndReceiverTerm( mesh, regionNames );

    NodeManager & nodeManager = mesh.getNodeManager();
    FaceManager & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords = nodeManager.getField< fields::referencePosition32 >().toViewConst();

    // mass matrix to be computed in this function
    arrayView1d< real32 > const mass = nodeManager.getField< fields::ElasticMassVector >();
    mass.zero();
    /// damping matrix to be computed for each dof in the boundary of the mesh
    arrayView1d< real32 > const dampingx = nodeManager.getField< fields::DampingVectorx >();
    arrayView1d< real32 > const dampingy = nodeManager.getField< fields::DampingVectory >();
    arrayView1d< real32 > const dampingz = nodeManager.getField< fields::DampingVectorz >();
    dampingx.zero();
    dampingy.zero();
    dampingz.zero();

    /// get array of indicators: 1 if face is on the free surface; 0 otherwise
    arrayView1d< localIndex const > const freeSurfaceFaceIndicator    = faceManager.getField< fields::ElasticFreeSurfaceFaceIndicator >();
    arrayView1d< integer const > const & facesDomainBoundaryIndicator = faceManager.getDomainBoundaryIndicator();
    ArrayOfArraysView< localIndex const > const facesToNodes          = faceManager.nodeList().toViewConst();
    arrayView2d< real64 const > const faceNormal                      = faceManager.faceNormal();

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                CellElementSubRegion & elementSubRegion )
    {
      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = elementSubRegion.nodeList();
      arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();

      computeTargetNodeSet( elemsToNodes, elementSubRegion.size(), fe.getNumQuadraturePoints() );

      arrayView1d< real32 const > const density = elementSubRegion.getField< fields::ElasticDensity >();
      arrayView1d< real32 const > const velocityVp = elementSubRegion.getField< fields::ElasticVelocityVp >();
      arrayView1d< real32 const > const velocityVs = elementSubRegion.getField< fields::ElasticVelocityVs >();

      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        elasticWaveEquationSEMKernels::MassMatrixKernel< FE_TYPE > kernelM( finiteElement );
        kernelM.template launch< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                               nodeCoords,
                                                               elemsToNodes,
                                                               density,
                                                               mass );

        elasticWaveEquationSEMKernels::DampingMatrixKernel< FE_TYPE > kernelD( finiteElement );
        kernelD.template launch< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                               nodeCoords,
                                                               elemsToFaces,
                                                               facesToNodes,
                                                               facesDomainBoundaryIndicator,
                                                               freeSurfaceFaceIndicator,
                                                               faceNormal,
                                                               density,
                                                               velocityVp,
                                                               velocityVs,
                                                               dampingx,
                                                               dampingy,
                                                               dampingz );
      } );
    } );
  } );

  WaveSolverUtils::initTrace( "seismoTraceReceiver", getName(), m_outputSeismoTrace, m_receiverConstants.size( 0 ), m_receiverIsLocal );
  WaveSolverUtils::initTrace( "dasTraceReceiver", getName(), m_outputSeismoTrace, m_linearDASGeometry.size( 0 ), m_receiverIsLocal );
}


void ElasticWaveEquationSEM::applyFreeSurfaceBC( real64 const time, DomainPartition & domain )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  FaceManager & faceManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getFaceManager();
  NodeManager & nodeManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getNodeManager();

  arrayView1d< real32 > const ux_np1 = nodeManager.getField< fields::Displacementx_np1 >();
  arrayView1d< real32 > const uy_np1 = nodeManager.getField< fields::Displacementy_np1 >();
  arrayView1d< real32 > const uz_np1 = nodeManager.getField< fields::Displacementz_np1 >();
  arrayView1d< real32 > const ux_n   = nodeManager.getField< fields::Displacementx_n >();
  arrayView1d< real32 > const uy_n   = nodeManager.getField< fields::Displacementy_n >();
  arrayView1d< real32 > const uz_n   = nodeManager.getField< fields::Displacementz_n >();
  arrayView1d< real32 > const ux_nm1 = nodeManager.getField< fields::Displacementx_nm1 >();
  arrayView1d< real32 > const uy_nm1 = nodeManager.getField< fields::Displacementy_nm1 >();
  arrayView1d< real32 > const uz_nm1 = nodeManager.getField< fields::Displacementz_nm1 >();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  /// set array of indicators: 1 if a face is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceFaceIndicator = faceManager.getField< fields::ElasticFreeSurfaceFaceIndicator >();

  /// set array of indicators: 1 if a node is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceNodeIndicator = nodeManager.getField< fields::ElasticFreeSurfaceNodeIndicator >();


  fsManager.apply( time,
                   domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ),
                   WaveSolverBase::viewKeyStruct::freeSurfaceString(),
                   [&]( FieldSpecificationBase const & bc,
                        string const &,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group &,
                        string const & )
  {
    string const & functionName = bc.getFunctionName();

    if( functionName.empty() || functionManager.getGroup< FunctionBase >( functionName ).isFunctionOfTime() == 2 )
    {
      real64 const value = bc.getScale();

      for( localIndex i = 0; i < targetSet.size(); ++i )
      {
        localIndex const kf = targetSet[ i ];
        freeSurfaceFaceIndicator[kf] = 1;

        localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );
        for( localIndex a=0; a < numNodes; ++a )
        {
          localIndex const dof = faceToNodeMap( kf, a );
          freeSurfaceNodeIndicator[dof] = 1;

          ux_np1[dof] = value;
          uy_np1[dof] = value;
          uz_np1[dof] = value;
          ux_n[dof] = value;
          uy_n[dof] = value;
          uz_n[dof] = value;
          ux_nm1[dof] = value;
          uy_nm1[dof] = value;
          uz_nm1[dof] = value;
        }
      }
    }
    else
    {
      GEOS_ERROR( getDataContext() << ": This option is not supported yet" );
    }
  } );
}



real64 ElasticWaveEquationSEM::explicitStepForward( real64 const & time_n,
                                                    real64 const & dt,
                                                    integer cycleNumber,
                                                    DomainPartition & domain,
                                                    bool GEOS_UNUSED_PARAM( computeGradient ) )
{
  real64 dtOut = explicitStepInternal( time_n, dt, cycleNumber, domain );
  return dtOut;
}



real64 ElasticWaveEquationSEM::explicitStepBackward( real64 const & time_n,
                                                     real64 const & dt,
                                                     integer cycleNumber,
                                                     DomainPartition & domain,
                                                     bool GEOS_UNUSED_PARAM( computeGradient ) )
{
  GEOS_ERROR( getDataContext() << ": Backward propagation for the elastic wave propagator not yet implemented" );
  real64 dtOut = explicitStepInternal( time_n, dt, cycleNumber, domain );
  return dtOut;
}

void ElasticWaveEquationSEM::computeUnknowns( real64 const &,
                                              real64 const & dt,
                                              integer const cycleNumber,
                                              DomainPartition &,
                                              MeshLevel & mesh,
                                              arrayView1d< string const > const & regionNames )
{
  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView1d< real32 const > const mass = nodeManager.getField< fields::ElasticMassVector >();
  arrayView1d< real32 const > const dampingx = nodeManager.getField< fields::DampingVectorx >();
  arrayView1d< real32 const > const dampingy = nodeManager.getField< fields::DampingVectory >();
  arrayView1d< real32 const > const dampingz = nodeManager.getField< fields::DampingVectorz >();
  arrayView1d< real32 > const stiffnessVectorx = nodeManager.getField< fields::StiffnessVectorx >();
  arrayView1d< real32 > const stiffnessVectory = nodeManager.getField< fields::StiffnessVectory >();
  arrayView1d< real32 > const stiffnessVectorz = nodeManager.getField< fields::StiffnessVectorz >();

  arrayView1d< real32 > const ux_nm1 = nodeManager.getField< fields::Displacementx_nm1 >();
  arrayView1d< real32 > const uy_nm1 = nodeManager.getField< fields::Displacementy_nm1 >();
  arrayView1d< real32 > const uz_nm1 = nodeManager.getField< fields::Displacementz_nm1 >();
  arrayView1d< real32 > const ux_n = nodeManager.getField< fields::Displacementx_n >();
  arrayView1d< real32 > const uy_n = nodeManager.getField< fields::Displacementy_n >();
  arrayView1d< real32 > const uz_n = nodeManager.getField< fields::Displacementz_n >();
  arrayView1d< real32 > const ux_np1 = nodeManager.getField< fields::Displacementx_np1 >();
  arrayView1d< real32 > const uy_np1 = nodeManager.getField< fields::Displacementy_np1 >();
  arrayView1d< real32 > const uz_np1 = nodeManager.getField< fields::Displacementz_np1 >();

  /// get array of indicators: 1 if node on free surface; 0 otherwise
  arrayView1d< localIndex const > const freeSurfaceNodeIndicator = nodeManager.getField< fields::ElasticFreeSurfaceNodeIndicator >();

  arrayView1d< real32 > const rhsx = nodeManager.getField< fields::ForcingRHSx >();
  arrayView1d< real32 > const rhsy = nodeManager.getField< fields::ForcingRHSy >();
  arrayView1d< real32 > const rhsz = nodeManager.getField< fields::ForcingRHSz >();

  auto kernelFactory = elasticWaveEquationSEMKernels::ExplicitElasticSEMFactory( dt );

  finiteElement::
    regionBasedKernelApplication< EXEC_POLICY,
                                  constitutive::NullModel,
                                  CellElementSubRegion >( mesh,
                                                          regionNames,
                                                          getDiscretizationName(),
                                                          "",
                                                          kernelFactory );


  addSourceToRightHandSide( cycleNumber, rhsx, rhsy, rhsz );

  real64 const dt2 = pow( dt, 2 );
  SortedArrayView< localIndex const > const solverTargetNodesSet = m_solverTargetNodesSet.toViewConst();
  forAll< EXEC_POLICY >( solverTargetNodesSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const n )
  {
    localIndex const a = solverTargetNodesSet[n];
    if( freeSurfaceNodeIndicator[a] != 1 )
    {
      ux_np1[a] = ux_n[a];
      ux_np1[a] *= 2.0*mass[a];
      ux_np1[a] -= (mass[a]-0.5*dt*dampingx[a])*ux_nm1[a];
      ux_np1[a] += dt2*(rhsx[a]-stiffnessVectorx[a]);
      ux_np1[a] /= mass[a]+0.5*dt*dampingx[a];
      uy_np1[a] = uy_n[a];
      uy_np1[a] *= 2.0*mass[a];
      uy_np1[a] -= (mass[a]-0.5*dt*dampingy[a])*uy_nm1[a];
      uy_np1[a] += dt2*(rhsy[a]-stiffnessVectory[a]);
      uy_np1[a] /= mass[a]+0.5*dt*dampingy[a];
      uz_np1[a] = uz_n[a];
      uz_np1[a] *= 2.0*mass[a];
      uz_np1[a] -= (mass[a]-0.5*dt*dampingz[a])*uz_nm1[a];
      uz_np1[a] += dt2*(rhsz[a]-stiffnessVectorz[a]);
      uz_np1[a] /= mass[a]+0.5*dt*dampingz[a];
    }
  } );
}

void ElasticWaveEquationSEM::synchronizeUnknowns( real64 const & time_n,
                                                  real64 const & dt,
                                                  integer const,
                                                  DomainPartition & domain,
                                                  MeshLevel & mesh,
                                                  arrayView1d< string const > const & )
{
  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView1d< real32 > const stiffnessVectorx = nodeManager.getField< fields::StiffnessVectorx >();
  arrayView1d< real32 > const stiffnessVectory = nodeManager.getField< fields::StiffnessVectory >();
  arrayView1d< real32 > const stiffnessVectorz = nodeManager.getField< fields::StiffnessVectorz >();

  arrayView1d< real32 > const ux_nm1 = nodeManager.getField< fields::Displacementx_nm1 >();
  arrayView1d< real32 > const uy_nm1 = nodeManager.getField< fields::Displacementy_nm1 >();
  arrayView1d< real32 > const uz_nm1 = nodeManager.getField< fields::Displacementz_nm1 >();
  arrayView1d< real32 > const ux_n   = nodeManager.getField< fields::Displacementx_n >();
  arrayView1d< real32 > const uy_n   = nodeManager.getField< fields::Displacementy_n >();
  arrayView1d< real32 > const uz_n   = nodeManager.getField< fields::Displacementz_n >();
  arrayView1d< real32 > const ux_np1 = nodeManager.getField< fields::Displacementx_np1 >();
  arrayView1d< real32 > const uy_np1 = nodeManager.getField< fields::Displacementy_np1 >();
  arrayView1d< real32 > const uz_np1 = nodeManager.getField< fields::Displacementz_np1 >();

  arrayView1d< real32 > const rhsx = nodeManager.getField< fields::ForcingRHSx >();
  arrayView1d< real32 > const rhsy = nodeManager.getField< fields::ForcingRHSy >();
  arrayView1d< real32 > const rhsz = nodeManager.getField< fields::ForcingRHSz >();

  /// synchronize displacement fields
  FieldIdentifiers fieldsToBeSync;
  fieldsToBeSync.addFields( FieldLocation::Node, { fields::Displacementx_np1::key(), fields::Displacementy_np1::key(), fields::Displacementz_np1::key() } );

  CommunicationTools & syncFields = CommunicationTools::getInstance();
  syncFields.synchronizeFields( fieldsToBeSync,
                                domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ),
                                domain.getNeighbors(),
                                true );

  // compute the seismic traces since last step.
  arrayView2d< real32 > const uXReceivers = m_displacementXNp1AtReceivers.toView();
  arrayView2d< real32 > const uYReceivers = m_displacementYNp1AtReceivers.toView();
  arrayView2d< real32 > const uZReceivers = m_displacementZNp1AtReceivers.toView();

  computeAllSeismoTraces( time_n, dt, ux_np1, ux_n, uXReceivers );
  computeAllSeismoTraces( time_n, dt, uy_np1, uy_n, uYReceivers );
  computeAllSeismoTraces( time_n, dt, uz_np1, uz_n, uZReceivers );

  incrementIndexSeismoTrace( time_n );
}

void ElasticWaveEquationSEM::prepareNextTimestep( MeshLevel & mesh )
{
  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView1d< real32 > const ux_nm1 = nodeManager.getField< fields::Displacementx_nm1 >();
  arrayView1d< real32 > const uy_nm1 = nodeManager.getField< fields::Displacementy_nm1 >();
  arrayView1d< real32 > const uz_nm1 = nodeManager.getField< fields::Displacementz_nm1 >();
  arrayView1d< real32 > const ux_n   = nodeManager.getField< fields::Displacementx_n >();
  arrayView1d< real32 > const uy_n   = nodeManager.getField< fields::Displacementy_n >();
  arrayView1d< real32 > const uz_n   = nodeManager.getField< fields::Displacementz_n >();
  arrayView1d< real32 > const ux_np1 = nodeManager.getField< fields::Displacementx_np1 >();
  arrayView1d< real32 > const uy_np1 = nodeManager.getField< fields::Displacementy_np1 >();
  arrayView1d< real32 > const uz_np1 = nodeManager.getField< fields::Displacementz_np1 >();

  arrayView1d< real32 > const stiffnessVectorx = nodeManager.getField< fields::StiffnessVectorx >();
  arrayView1d< real32 > const stiffnessVectory = nodeManager.getField< fields::StiffnessVectory >();
  arrayView1d< real32 > const stiffnessVectorz = nodeManager.getField< fields::StiffnessVectorz >();

  arrayView1d< real32 > const rhsx = nodeManager.getField< fields::ForcingRHSx >();
  arrayView1d< real32 > const rhsy = nodeManager.getField< fields::ForcingRHSy >();
  arrayView1d< real32 > const rhsz = nodeManager.getField< fields::ForcingRHSz >();

  SortedArrayView< localIndex const > const solverTargetNodesSet = m_solverTargetNodesSet.toViewConst();

  forAll< EXEC_POLICY >( solverTargetNodesSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const n )
  {
    localIndex const a = solverTargetNodesSet[n];
    ux_nm1[a] = ux_n[a];
    uy_nm1[a] = uy_n[a];
    uz_nm1[a] = uz_n[a];
    ux_n[a] = ux_np1[a];
    uy_n[a] = uy_np1[a];
    uz_n[a] = uz_np1[a];

    stiffnessVectorx[a] = stiffnessVectory[a] = stiffnessVectorz[a] = 0.0;
    rhsx[a] = rhsy[a] = rhsz[a] = 0.0;
  } );

}

real64 ElasticWaveEquationSEM::explicitStepInternal( real64 const & time_n,
                                                     real64 const & dt,
                                                     integer const cycleNumber,
                                                     DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  GEOS_LOG_RANK_0_IF( dt < epsilonLoc, "Warning! Value for dt: " << dt << "s is smaller than local threshold: " << epsilonLoc );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    computeUnknowns( time_n, dt, cycleNumber, domain, mesh, regionNames );
    synchronizeUnknowns( time_n, dt, cycleNumber, domain, mesh, regionNames );
    prepareNextTimestep( mesh );
  } );

  return dt;
}

void ElasticWaveEquationSEM::cleanup( real64 const time_n,
                                      integer const cycleNumber,
                                      integer const eventCounter,
                                      real64 const eventProgress,
                                      DomainPartition & domain )
{
  // call the base class cleanup (for reporting purposes)
  SolverBase::cleanup( time_n, cycleNumber, eventCounter, eventProgress, domain );

  // compute the remaining seismic traces, if needed
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();
    arrayView1d< real32 const > const ux_n   = nodeManager.getField< fields::Displacementx_n >();
    arrayView1d< real32 const > const ux_np1 = nodeManager.getField< fields::Displacementx_np1 >();
    arrayView1d< real32 const > const uy_n   = nodeManager.getField< fields::Displacementy_n >();
    arrayView1d< real32 const > const uy_np1 = nodeManager.getField< fields::Displacementy_np1 >();
    arrayView1d< real32 const > const uz_n   = nodeManager.getField< fields::Displacementz_n >();
    arrayView1d< real32 const > const uz_np1 = nodeManager.getField< fields::Displacementz_np1 >();
    arrayView2d< real32 > const uXReceivers  = m_displacementXNp1AtReceivers.toView();
    arrayView2d< real32 > const uYReceivers  = m_displacementYNp1AtReceivers.toView();
    arrayView2d< real32 > const uZReceivers  = m_displacementZNp1AtReceivers.toView();

    computeAllSeismoTraces( time_n, 0.0, ux_np1, ux_n, uXReceivers );
    computeAllSeismoTraces( time_n, 0.0, uy_np1, uy_n, uYReceivers );
    computeAllSeismoTraces( time_n, 0.0, uz_np1, uz_n, uZReceivers );

    if( m_useDAS <= 0 )
    {
      WaveSolverUtils::writeSeismoTraceVector( "seismoTraceReceiver", getName(), m_outputSeismoTrace, m_receiverConstants.size( 0 ),
                                               m_receiverIsLocal, m_nsamplesSeismoTrace, uXReceivers, uYReceivers, uZReceivers );
    }

    /// Compute DAS data if requested
    /// Pairs or sets of receivers are assumed to be modeled ( see WaveSolverBase::initializeDAS() )
    if( m_useDAS > 0 )
    {
      computeDAS( uXReceivers, uYReceivers, uZReceivers );
      WaveSolverUtils::writeSeismoTrace( "dasTraceReceiver", getName(), m_outputSeismoTrace, m_linearDASGeometry.size( 0 ),
                                         m_receiverIsLocal, m_nsamplesSeismoTrace, uXReceivers );
    }
  } );
}

void ElasticWaveEquationSEM::initializePML()
{
  GEOS_ERROR( getDataContext() << ": PML for the elastic wave propagator not yet implemented" );
}

void ElasticWaveEquationSEM::applyPML( real64 const, DomainPartition & )
{
  GEOS_ERROR( getDataContext() << ": PML for the elastic wave propagator not yet implemented" );
}

REGISTER_CATALOG_ENTRY( SolverBase, ElasticWaveEquationSEM, string const &, dataRepository::Group * const )

} /* namespace geos */
