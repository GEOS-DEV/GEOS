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

  registerWrapper( viewKeyStruct::sourceNodeIdsString(), &m_sourceNodeIds ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Indices of the nodes (in the right order) for each source point" );

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

  registerWrapper( viewKeyStruct::sourceIsAccessibleString(), &m_sourceIsAccessible ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Flag that indicates whether the source is accessible to this MPI rank" );

  registerWrapper( viewKeyStruct::receiverNodeIdsString(), &m_receiverNodeIds ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Indices of the nodes (in the right order) for each receiver point" );

  registerWrapper( viewKeyStruct::sourceConstantsString(), &m_receiverConstants ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Constant part of the receiver for the nodes listed in m_receiverNodeIds" );

  registerWrapper( viewKeyStruct::receiverIsLocalString(), &m_receiverIsLocal ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Flag that indicates whether the receiver is local to this MPI rank" );

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

localIndex ElasticWaveEquationSEM::getNumNodesPerElem()
{
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();

  FiniteElementDiscretizationManager const &
  feDiscretizationManager = numericalMethodManager.getFiniteElementDiscretizationManager();

  FiniteElementDiscretization const * const
  feDiscretization = feDiscretizationManager.getGroupPointer< FiniteElementDiscretization >( m_discretizationName );
  GEOSX_THROW_IF_IF( feDiscretization == nullptr,
                  getName() << ": FE discretization not found: " << m_discretizationName,
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

void ElasticWaveEquationSEM::initializePreSubGroups()
{

  WaveSolverBase::initializePreSubGroups();

  localIndex const numNodesPerElem = getNumNodesPerElem();

  localIndex const numSourcesGlobal = m_sourceCoordinates.size( 0 );
  m_sourceNodeIds.resize( numSourcesGlobal, numNodesPerElem );
  m_sourceConstantsx.resize( numSourcesGlobal, numNodesPerElem );
  m_sourceConstantsy.resize( numSourcesGlobal, numNodesPerElem );
  m_sourceConstantsz.resize( numSourcesGlobal, numNodesPerElem );

  localIndex const numReceiversGlobal = m_receiverCoordinates.size( 0 );
  m_receiverNodeIds.resize( numReceiversGlobal, numNodesPerElem );
  m_receiverConstants.resize( numReceiversGlobal, numNodesPerElem );

}


void ElasticWaveEquationSEM::registerDataOnMesh( Group & meshBodies )
{

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
                               fields::MassVector,
                               fields::DampingVectorx,
                               fields::DampingVectory,
                               fields::DampingVectorz,
                               fields::StiffnessVectorx,
                               fields::StiffnessVectory,
                               fields::StiffnessVectorz,
                               fields::FreeSurfaceNodeIndicator >( this->getName() );

    FaceManager & faceManager = mesh.getFaceManager();
    faceManager.registerField< fields::FreeSurfaceFaceIndicator >( this->getName() );

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerField< fields::MediumVelocityVp >( this->getName() );
      subRegion.registerField< fields::MediumVelocityVs >( this->getName() );
      subRegion.registerField< fields::MediumDensity >( this->getName() );
    } );

  } );
}



void ElasticWaveEquationSEM::postProcessInput()
{

  WaveSolverBase::postProcessInput();

  GEOS_ERROR_IF( m_sourceCoordinates.size( 1 ) != 3,
                  "Invalid number of physical coordinates for the sources" );

  GEOS_ERROR_IF( m_receiverCoordinates.size( 1 ) != 3,
                  "Invalid number of physical coordinates for the receivers" );

  EventManager const & event = this->getGroupByPath< EventManager >( "/Problem/Events" );
  real64 const & maxTime = event.getReference< real64 >( EventManager::viewKeyStruct::maxTimeString() );
  real64 dt = 0;
  for( localIndex numSubEvent = 0; numSubEvent < event.numSubGroups(); ++numSubEvent )
  {
    EventBase const * subEvent = static_cast< EventBase const * >( event.getSubGroups()[numSubEvent] );
    if( subEvent->getEventName() == "/Solvers/" + this->getName() )
    {
      dt = subEvent->getReference< real64 >( EventBase::viewKeyStruct::forceDtString() );
    }
  }

  GEOSX_THROW_IF_IF( dt < epsilonLoc*maxTime, "Value for dt: " << dt <<" is smaller than local threshold: " << epsilonLoc, std::runtime_error );

  if( m_dtSeismoTrace > 0 )
  {
    m_nsamplesSeismoTrace = int( maxTime / m_dtSeismoTrace) + 1;
  }
  else
  {
    m_nsamplesSeismoTrace = 0;
  }
  localIndex const nsamples = int(maxTime/dt) + 1;

  localIndex const numSourcesGlobal = m_sourceCoordinates.size( 0 );
  m_sourceIsAccessible.resize( numSourcesGlobal );

  localIndex const numReceiversGlobal = m_receiverCoordinates.size( 0 );
  m_receiverIsLocal.resize( numReceiversGlobal );

  m_displacementXNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );
  m_displacementYNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );
  m_displacementZNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );
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

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X =
    nodeManager.referencePosition().toViewConst();
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

  real32 const timeSourceFrequency = this->m_timeSourceFrequency;
  localIndex const rickerOrder = this->m_rickerOrder;
  arrayView2d< real32 > const sourceValue = m_sourceValue.toView();

  real64 dt = 0;
  EventManager const & event = this->getGroupByPath< EventManager >( "/Problem/Events" );
  for( localIndex numSubEvent = 0; numSubEvent < event.numSubGroups(); ++numSubEvent )
  {
    EventBase const * subEvent = static_cast< EventBase const * >( event.getSubGroups()[numSubEvent] );
    if( subEvent->getEventName() == "/Solvers/" + this->getName() )
    {
      dt = subEvent->getReference< real64 >( EventBase::viewKeyStruct::forceDtString() );
    }
  }

  mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                        CellElementSubRegion & elementSubRegion )
  {

    GEOSX_THROW_IF_IF( elementSubRegion.getElementType() != ElementType::Hexahedron,
                    "Invalid type of element, the elastic solver is designed for hexahedral meshes only (C3D8) ",
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
        timeSourceFrequency,
        rickerOrder,
        m_sourceForce,
        m_sourceMoment );
    } );
  } );
}

void ElasticWaveEquationSEM::computeDAS ( arrayView2d< real32 > const xCompRcv,
                                          arrayView2d< real32 > const yCompRcv,
                                          arrayView2d< real32 > const zCompRcv )
{

  arrayView2d< real64 const > const linearDASGeometry = m_linearDASGeometry.toViewConst();
  arrayView1d< localIndex const > const receiverIsLocal = m_receiverIsLocal.toViewConst();

  localIndex const numReceiversGlobal = linearDASGeometry.size( 0 );
  localIndex const nsamplesSeismoTrace = m_nsamplesSeismoTrace;

  if( m_nsamplesSeismoTrace > 0 )
  {
    /// synchronize receivers across MPI ranks
    MpiWrapper::allReduce( xCompRcv.data(),
                           xCompRcv.data(),
                           2*numReceiversGlobal*m_nsamplesSeismoTrace,
                           MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
                           MPI_COMM_GEOSX );

    MpiWrapper::allReduce( yCompRcv.data(),
                           yCompRcv.data(),
                           2*numReceiversGlobal*m_nsamplesSeismoTrace,
                           MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
                           MPI_COMM_GEOSX );

    MpiWrapper::allReduce( zCompRcv.data(),
                           zCompRcv.data(),
                           2*numReceiversGlobal*m_nsamplesSeismoTrace,
                           MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
                           MPI_COMM_GEOSX );

    forAll< EXEC_POLICY >( numReceiversGlobal, [=] GEOS_HOST_DEVICE ( localIndex const ircv )
    {
      if( receiverIsLocal[ircv] == 1 )
      {
        real32 const cd = cos( linearDASGeometry[ircv][0] );
        real32 const sd = sin( linearDASGeometry[ircv][0] );
        real32 const ca = cos( linearDASGeometry[ircv][1] );
        real32 const sa = sin( linearDASGeometry[ircv][1] );

        /// convert dipole data (pairs of geophones) to average strain data and
        for( localIndex iSample = 0; iSample < nsamplesSeismoTrace; ++iSample )
        {
          // store strain data in the z-component of the receiver (copied to x after resize)
          zCompRcv[iSample][ircv] =
            cd * ca * ( xCompRcv[iSample][numReceiversGlobal+ircv] - xCompRcv[iSample][ircv] )
            + cd * sa * ( yCompRcv[iSample][numReceiversGlobal+ircv] - yCompRcv[iSample][ircv] )
            + sd * ( zCompRcv[iSample][numReceiversGlobal+ircv] - zCompRcv[iSample][ircv] );
          zCompRcv[iSample][ircv] /= linearDASGeometry[ircv][2];

        }
      }
    } );
  }

  /// temporary output to txt
  if( this->m_outputSeismoTrace == 1 )
  {
    forAll< serialPolicy >( numReceiversGlobal, [=] ( localIndex const ircv )
    {
      if( receiverIsLocal[ircv] == 1 )
      {
        std::ofstream f( GEOSX_FMT( "dasTraceReceiver{:03}.txt", ircv ), std::ios::app );
        for( localIndex iSample = 0; iSample < nsamplesSeismoTrace; ++iSample )
        {
          f<< iSample << " " << zCompRcv[iSample][ircv] << std::endl;
        }
        f.close();
      }
    } );
  }

  /// resize the receiver arrays by dropping the extra pair to avoid confusion
  /// the remaining x-component contains DAS data, the other components are set to zero
  m_displacementXNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );
  arrayView2d< real32 > const dasReceiver = m_displacementXNp1AtReceivers.toView();
  forAll< EXEC_POLICY >( numReceiversGlobal, [=] GEOS_HOST_DEVICE ( localIndex const ircv )
  {
    if( receiverIsLocal[ircv] == 1 )
    {
      /// convert dipole data (pairs of geophones) to average strain data and
      for( localIndex iSample = 0; iSample < nsamplesSeismoTrace; ++iSample )
      {
        // store strain data in the z-component of the receiver (copied to x after resize)
        dasReceiver[iSample][ircv] = zCompRcv[iSample][ircv];
      }
    }
  } );
  /// set the y and z components to zero to avoid any confusion
  m_displacementYNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );
  m_displacementYNp1AtReceivers.zero();
  m_displacementZNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );
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

  GEOSX_THROW_IF_IF( cycleNumber > sourceValue.size( 0 ), "Too many steps compared to array size", std::runtime_error );
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

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  real64 const time = 0.0;
  applyFreeSurfaceBC( time, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    precomputeSourceAndReceiverTerm( mesh, regionNames );

    NodeManager & nodeManager = mesh.getNodeManager();
    FaceManager & faceManager = mesh.getFaceManager();

    /// get the array of indicators: 1 if the face is on the boundary; 0 otherwise
    arrayView1d< integer const > const & facesDomainBoundaryIndicator = faceManager.getDomainBoundaryIndicator();
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();
    arrayView2d< real64 const > const faceNormal  = faceManager.faceNormal();

    /// get face to node map
    ArrayOfArraysView< localIndex const > const facesToNodes = faceManager.nodeList().toViewConst();

    // mass matrix to be computed in this function
    arrayView1d< real32 > const mass = nodeManager.getField< fields::MassVector >();
    mass.zero();
    /// damping matrix to be computed for each dof in the boundary of the mesh
    arrayView1d< real32 > const dampingx = nodeManager.getField< fields::DampingVectorx >();
    arrayView1d< real32 > const dampingy = nodeManager.getField< fields::DampingVectory >();
    arrayView1d< real32 > const dampingz = nodeManager.getField< fields::DampingVectorz >();
    dampingx.zero();
    dampingy.zero();
    dampingz.zero();

    /// get array of indicators: 1 if face is on the free surface; 0 otherwise
    arrayView1d< localIndex const > const freeSurfaceFaceIndicator = faceManager.getField< fields::FreeSurfaceFaceIndicator >();

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          CellElementSubRegion & elementSubRegion )
    {

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = elementSubRegion.nodeList();
      arrayView2d< localIndex const > const facesToElements = faceManager.elementList();
      arrayView1d< real32 > const density = elementSubRegion.getField< fields::MediumDensity >();
      arrayView1d< real32 > const velocityVp = elementSubRegion.getField< fields::MediumVelocityVp >();
      arrayView1d< real32 > const velocityVs = elementSubRegion.getField< fields::MediumVelocityVs >();

      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe,
                                                                               [&]
                                                                                 ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        elasticWaveEquationSEMKernels::MassMatrixKernel< FE_TYPE > kernelM( finiteElement );

        kernelM.template launch< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                               X,
                                                               elemsToNodes,
                                                               density,
                                                               mass );

        elasticWaveEquationSEMKernels::DampingMatrixKernel< FE_TYPE > kernelD( finiteElement );

        kernelD.template launch< EXEC_POLICY, ATOMIC_POLICY >( faceManager.size(),
                                                               X,
                                                               facesToElements,
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
  arrayView1d< real32 > const ux_n = nodeManager.getField< fields::Displacementx_n >();
  arrayView1d< real32 > const uy_n = nodeManager.getField< fields::Displacementy_n >();
  arrayView1d< real32 > const uz_n = nodeManager.getField< fields::Displacementz_n >();
  arrayView1d< real32 > const ux_nm1 = nodeManager.getField< fields::Displacementx_nm1 >();
  arrayView1d< real32 > const uy_nm1 = nodeManager.getField< fields::Displacementy_nm1 >();
  arrayView1d< real32 > const uz_nm1 = nodeManager.getField< fields::Displacementz_nm1 >();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  /// set array of indicators: 1 if a face is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceFaceIndicator = faceManager.getField< fields::FreeSurfaceFaceIndicator >();

  /// set array of indicators: 1 if a node is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceNodeIndicator = nodeManager.getField< fields::FreeSurfaceNodeIndicator >();


  freeSurfaceFaceIndicator.zero();
  freeSurfaceNodeIndicator.zero();

  fsManager.apply( time,
                   domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ),
                   string( "FreeSurface" ),
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
      GEOS_ERROR( "This option is not supported yet" );
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
  GEOS_ERROR( "Backward propagation for the elastic wave propagator not yet implemented" );
  real64 dtOut = explicitStepInternal( time_n, dt, cycleNumber, domain );
  return dtOut;
}

real64 ElasticWaveEquationSEM::explicitStepInternal( real64 const & time_n,
                                                     real64 const & dt,
                                                     integer const cycleNumber,
                                                     DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  GEOS_UNUSED_VAR( time_n, dt, cycleNumber );

  GEOS_LOG_RANK_0_IF( dt < epsilonLoc, "Warning! Value for dt: " << dt << "s is smaller than local threshold: " << epsilonLoc );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< real32 const > const mass = nodeManager.getField< fields::MassVector >();
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
    arrayView1d< localIndex const > const freeSurfaceNodeIndicator = nodeManager.getField< fields::FreeSurfaceNodeIndicator >();

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



    real64 const dt2 = dt*dt;
    forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
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

    /// synchronize pressure fields
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addFields( FieldLocation::Node, { fields::Displacementx_np1::key(), fields::Displacementy_np1::key(), fields::Displacementz_np1::key() } );

    CommunicationTools & syncFields = CommunicationTools::getInstance();
    syncFields.synchronizeFields( fieldsToBeSync,
                                  domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ),
                                  domain.getNeighbors(),
                                  true );

    // compute the seismic traces since last step.
    arrayView2d< real32 > const uXReceivers   = m_displacementXNp1AtReceivers.toView();
    arrayView2d< real32 > const uYReceivers   = m_displacementYNp1AtReceivers.toView();
    arrayView2d< real32 > const uZReceivers   = m_displacementZNp1AtReceivers.toView();

    computeAllSeismoTraces( time_n, dt, ux_np1, ux_n, uXReceivers );
    computeAllSeismoTraces( time_n, dt, uy_np1, uy_n, uYReceivers );
    computeAllSeismoTraces( time_n, dt, uz_np1, uz_n, uZReceivers );

    forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      ux_nm1[a] = ux_n[a];
      uy_nm1[a] = uy_n[a];
      uz_nm1[a] = uz_n[a];
      ux_n[a] = ux_np1[a];
      uy_n[a] = uy_np1[a];
      uz_n[a] = uz_np1[a];

      stiffnessVectorx[a] = 0.0;
      stiffnessVectory[a] = 0.0;
      stiffnessVectorz[a] = 0.0;
      rhsx[a] = 0.0;
      rhsy[a] = 0.0;
      rhsz[a] = 0.0;
    } );

    // increment m_indexSeismoTrace
    while( (m_dtSeismoTrace*m_indexSeismoTrace) <= (time_n + epsilonLoc) && m_indexSeismoTrace < m_nsamplesSeismoTrace )
    {
      m_indexSeismoTrace++;
    }

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
    arrayView1d< real32 const > const ux_n = nodeManager.getField< fields::Displacementx_n >();
    arrayView1d< real32 const > const ux_np1 = nodeManager.getField< fields::Displacementx_np1 >();
    arrayView1d< real32 const > const uy_n = nodeManager.getField< fields::Displacementy_n >();
    arrayView1d< real32 const > const uy_np1 = nodeManager.getField< fields::Displacementy_np1 >();
    arrayView1d< real32 const > const uz_n = nodeManager.getField< fields::Displacementz_n >();
    arrayView1d< real32 const > const uz_np1 = nodeManager.getField< fields::Displacementz_np1 >();
    arrayView2d< real32 > const uXReceivers   = m_displacementXNp1AtReceivers.toView();
    arrayView2d< real32 > const uYReceivers   = m_displacementYNp1AtReceivers.toView();
    arrayView2d< real32 > const uZReceivers   = m_displacementZNp1AtReceivers.toView();

    computeAllSeismoTraces( time_n, 0, ux_np1, ux_n, uXReceivers );
    computeAllSeismoTraces( time_n, 0, uy_np1, uy_n, uYReceivers );
    computeAllSeismoTraces( time_n, 0, uz_np1, uz_n, uZReceivers );

    /// Compute DAS data if requested
    /// Pairs of receivers are assumed to be modeled ( see WaveSolverBase::initializeDAS() )
    if( m_useDAS )
    {
      computeDAS( uXReceivers, uYReceivers, uZReceivers );
    }
  } );

  // increment m_indexSeismoTrace
  while( (m_dtSeismoTrace*m_indexSeismoTrace) <= (time_n + epsilonLoc) && m_indexSeismoTrace < m_nsamplesSeismoTrace )
  {
    m_indexSeismoTrace++;
  }
}

void ElasticWaveEquationSEM::computeAllSeismoTraces( real64 const time_n,
                                                     real64 const dt,
                                                     arrayView1d< real32 const > const var_np1,
                                                     arrayView1d< real32 const > const var_n,
                                                     arrayView2d< real32 > varAtReceivers )
{
  localIndex indexSeismoTrace = m_indexSeismoTrace;

  for( real64 timeSeismo;
       (timeSeismo = m_dtSeismoTrace*indexSeismoTrace) <= (time_n + epsilonLoc) && indexSeismoTrace < m_nsamplesSeismoTrace;
       indexSeismoTrace++ )
  {
    WaveSolverUtils::computeSeismoTrace( time_n, dt, timeSeismo, indexSeismoTrace, m_receiverNodeIds, m_receiverConstants, m_receiverIsLocal,
                                         m_nsamplesSeismoTrace, m_outputSeismoTrace, var_np1, var_n, varAtReceivers );
  }
}

void ElasticWaveEquationSEM::initializePML()
{
  GEOS_ERROR( "PML for the elastic wave propagator not yet implemented" );
}

void ElasticWaveEquationSEM::applyPML( real64 const, DomainPartition & )
{
  GEOS_ERROR( "PML for the elastic wave propagator not yet implemented" );
}

REGISTER_CATALOG_ENTRY( SolverBase, ElasticWaveEquationSEM, string const &, dataRepository::Group * const )

} /* namespace geos */
