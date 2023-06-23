/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOS Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file AcousticVTIWaveEquationSEM.cpp
 */

#include "AcousticVTIWaveEquationSEM.hpp"
#include "AcousticVTIWaveEquationSEMKernel.hpp"

#include "finiteElement/FiniteElementDiscretization.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/PerfectlyMatchedLayer.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/ElementType.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "WaveSolverUtils.hpp"

namespace geos
{

using namespace dataRepository;

AcousticVTIWaveEquationSEM::AcousticVTIWaveEquationSEM( const std::string & name,
                                                  Group * const parent ):
  WaveSolverBase( name,
                  parent )
{

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
    setDescription( "Flag that indicates whether the source is accessible to this MPI rank" );

  registerWrapper( viewKeyStruct::receiverNodeIdsString(), &m_receiverNodeIds ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Indices of the nodes (in the right order) for each receiver point" );

  registerWrapper( viewKeyStruct::receiverConstantsString(), &m_receiverConstants ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Constant part of the receiver for the nodes listed in m_receiverNodeIds" );

  registerWrapper( viewKeyStruct::receiverIsLocalString(), &m_receiverIsLocal ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Flag that indicates whether the receiver is local to this MPI rank" );

  registerWrapper( viewKeyStruct::pressureNp1AtReceiversString(), &m_pressureNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Pressure value at each receiver for each timestep" );

}

AcousticVTIWaveEquationSEM::~AcousticVTIWaveEquationSEM()
{
  // TODO Auto-generated destructor stub
}

localIndex AcousticVTIWaveEquationSEM::getNumNodesPerElem()
{
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();

  FiniteElementDiscretizationManager const &
  feDiscretizationManager = numericalMethodManager.getFiniteElementDiscretizationManager();

  FiniteElementDiscretization const * const
  feDiscretization = feDiscretizationManager.getGroupPointer< FiniteElementDiscretization >( m_discretizationName );
  GEOS_THROW_IF( feDiscretization == nullptr,
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

void AcousticVTIWaveEquationSEM::initializePreSubGroups()
{
  WaveSolverBase::initializePreSubGroups();

  localIndex numNodesPerElem = getNumNodesPerElem();

  localIndex const numSourcesGlobal = m_sourceCoordinates.size( 0 );
  m_sourceNodeIds.resize( numSourcesGlobal, numNodesPerElem );
  m_sourceConstants.resize( numSourcesGlobal, numNodesPerElem );

  localIndex const numReceiversGlobal = m_receiverCoordinates.size( 0 );
  m_receiverNodeIds.resize( numReceiversGlobal, numNodesPerElem );
  m_receiverConstants.resize( numReceiversGlobal, numNodesPerElem );

}


void AcousticVTIWaveEquationSEM::registerDataOnMesh( Group & meshBodies )
{

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    nodeManager.registerField< fields::Pressure_p_nm1,
                               fields::Pressure_p_n,
                               fields::Pressure_p_np1,
                               fields::Pressure_q_nm1,
                               fields::Pressure_q_n,
                               fields::Pressure_q_np1,
                               fields::PressureDoubleDerivative,
                               fields::ForcingRHS,
                               fields::MassVector,
                               fields::DampingVector_p,
                               fields::DampingVector_pq,
                               fields::DampingVector_q,
                               fields::DampingVector_qp,
                               fields::StiffnessVector_p,
                               fields::StiffnessVector_q,
                               fields::StiffnessVector, //Debug
                               fields::DampingVector, //Debug
                               fields::FreeSurfaceNodeIndicator,
                               fields::LateralSurfaceNodeIndicator,
                               fields::BottomSurfaceNodeIndicator >( this->getName() );


    FaceManager & faceManager = mesh.getFaceManager();
    faceManager.registerField< fields::FreeSurfaceFaceIndicator >( this->getName() );
    faceManager.registerField< fields::LateralSurfaceFaceIndicator >( this->getName() );
    faceManager.registerField< fields::BottomSurfaceFaceIndicator >( this->getName() );

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerField< fields::Delta >( this->getName() );
      subRegion.registerField< fields::Epsilon >( this->getName() );
      subRegion.registerField< fields::F >( this->getName() );
      subRegion.registerField< fields::MediumVelocity >( this->getName() );
      subRegion.registerField< fields::PartialGradient >( this->getName() );
    } );

    arrayView1d< real32 > const p_dt2 = nodeManager.getField< fields::PressureDoubleDerivative >();
    int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
    std::string lifoPrefix = GEOS_FMT( "lifo/rank_{:05}/pdt2_shot{:06}", rank, m_shotIndex );
    m_lifo = std::unique_ptr< lifoStorage< real32 > >( new lifoStorage< real32 >( lifoPrefix, p_dt2, m_lifoOnDevice, m_lifoOnHost, m_lifoSize ) );

  } );
}


void AcousticVTIWaveEquationSEM::postProcessInput()
{

  WaveSolverBase::postProcessInput();
  GEOS_THROW_IF( m_sourceCoordinates.size( 1 ) != 3,
                  "Invalid number of physical coordinates for the sources",
                  InputError );

  GEOS_THROW_IF( m_receiverCoordinates.size( 1 ) != 3,
                  "Invalid number of physical coordinates for the receivers",
                  InputError );

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

  GEOS_THROW_IF( dt < epsilonLoc*maxTime, "Value for dt: " << dt <<" is smaller than local threshold: " << epsilonLoc, std::runtime_error );

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

  m_pressureNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );
  m_sourceValue.resize( nsamples, numSourcesGlobal );

}

void AcousticVTIWaveEquationSEM::precomputeSourceAndReceiverTerm( MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
{
  NodeManager const & nodeManager = mesh.getNodeManager();
  FaceManager const & faceManager = mesh.getFaceManager();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const
  X = nodeManager.referencePosition().toViewConst();

  arrayView2d< real64 const > const faceNormal  = faceManager.faceNormal();
  arrayView2d< real64 const > const faceCenter  = faceManager.faceCenter();


  arrayView2d< real64 const > const sourceCoordinates = m_sourceCoordinates.toViewConst();
  arrayView2d< localIndex > const sourceNodeIds = m_sourceNodeIds.toView();
  arrayView2d< real64 > const sourceConstants = m_sourceConstants.toView();
  arrayView1d< localIndex > const sourceIsAccessible = m_sourceIsAccessible.toView();
  sourceNodeIds.setValues< EXEC_POLICY >( -1 );
  sourceConstants.setValues< EXEC_POLICY >( -1 );
  sourceIsAccessible.zero();

  arrayView2d< real64 const > const receiverCoordinates = m_receiverCoordinates.toViewConst();
  arrayView2d< localIndex > const receiverNodeIds = m_receiverNodeIds.toView();
  arrayView2d< real64 > const receiverConstants = m_receiverConstants.toView();
  arrayView1d< localIndex > const receiverIsLocal = m_receiverIsLocal.toView();
  receiverNodeIds.setValues< EXEC_POLICY >( -1 );
  receiverConstants.setValues< EXEC_POLICY >( -1 );
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
    GEOS_THROW_IF( elementSubRegion.getElementType() != ElementType::Hexahedron,
                    "Invalid type of element, the acoustic solver is designed for hexahedral meshes only (C3D8), using the SEM formulation",
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

      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      localIndex const numFacesPerElem = elementSubRegion.numFacesPerElement();

      acousticVTIWaveEquationSEMKernels::
        PrecomputeSourceAndReceiverKernel::
        launch< EXEC_POLICY, FE_TYPE >
        ( elementSubRegion.size(),
        numNodesPerElem,
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
        sourceConstants,
        receiverCoordinates,
        receiverIsLocal,
        receiverNodeIds,
        receiverConstants,
        sourceValue,
        dt,
        timeSourceFrequency,
        rickerOrder );
    } );
  } );
}

void AcousticVTIWaveEquationSEM::addSourceToRightHandSide( integer const & cycleNumber, arrayView1d< real32 > const rhs )
{
  arrayView2d< localIndex const > const sourceNodeIds = m_sourceNodeIds.toViewConst();
  arrayView2d< real64 const > const sourceConstants   = m_sourceConstants.toViewConst();
  arrayView1d< localIndex const > const sourceIsAccessible = m_sourceIsAccessible.toViewConst();
  arrayView2d< real32 const > const sourceValue   = m_sourceValue.toViewConst();

  GEOS_THROW_IF( cycleNumber > sourceValue.size( 0 ), "Too many steps compared to array size", std::runtime_error );
  forAll< EXEC_POLICY >( sourceConstants.size( 0 ), [=] GEOS_HOST_DEVICE ( localIndex const isrc )
  {
    if( sourceIsAccessible[isrc] == 1 )
    {
      for( localIndex inode = 0; inode < sourceConstants.size( 1 ); ++inode )
      {
        real32 const localIncrement = sourceConstants[isrc][inode] * sourceValue[cycleNumber][isrc];
        RAJA::atomicAdd< ATOMIC_POLICY >( &rhs[sourceNodeIds[isrc][inode]], localIncrement );
      }
    }
  } );
}

void AcousticVTIWaveEquationSEM::initializePostInitialConditionsPreSubGroups()
{

  WaveSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  real64 const time = 0.0;
  applyFreeSurfaceBC( time, domain );
  precomputeSurfaceFieldIndicator( time, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    precomputeSourceAndReceiverTerm( mesh, regionNames );

    NodeManager & nodeManager = mesh.getNodeManager();
    FaceManager & faceManager = mesh.getFaceManager();

    /// get the array of indicators: 1 if the face is on the boundary; 0 otherwise
    arrayView1d< integer > const & facesDomainBoundaryIndicator = faceManager.getDomainBoundaryIndicator();
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

    /// get face to node map
    ArrayOfArraysView< localIndex const > const facesToNodes = faceManager.nodeList().toViewConst();

    // mass matrix to be computed in this function
    arrayView1d< real32 > const mass = nodeManager.getField< fields::MassVector >();
    mass.zero();
    /// damping matrices to be computed for each dof in the boundary of the mesh
    arrayView1d< real32 > const damping_p  = nodeManager.getField< fields::DampingVector_p >();
    arrayView1d< real32 > const damping_pq = nodeManager.getField< fields::DampingVector_pq >();
    arrayView1d< real32 > const damping_q  = nodeManager.getField< fields::DampingVector_q >();
    arrayView1d< real32 > const damping_qp = nodeManager.getField< fields::DampingVector_qp >();
    damping_p.zero();
    damping_pq.zero();
    damping_q.zero();
    damping_qp.zero();
    //DEBUG    
    arrayView1d< real32 > const damping = nodeManager.getField< fields::DampingVector >();
    damping.zero();

    /// get array of indicators: 1 if face is on the free surface; 0 otherwise
    arrayView1d< localIndex const > const freeSurfaceFaceIndicator = faceManager.getField< fields::FreeSurfaceFaceIndicator >();
    arrayView1d< localIndex const > const lateralSurfaceFaceIndicator = faceManager.getField< fields::LateralSurfaceFaceIndicator >();
    arrayView1d< localIndex const > const bottomSurfaceFaceIndicator = faceManager.getField< fields::BottomSurfaceFaceIndicator >();

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          CellElementSubRegion & elementSubRegion )
    {

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = elementSubRegion.nodeList();
      arrayView2d< localIndex const > const facesToElements = faceManager.elementList();
      arrayView1d< real32 const > const velocity = elementSubRegion.getField< fields::MediumVelocity >();
      arrayView1d< real32 const > const epsilon  = elementSubRegion.getField< fields::Epsilon >();
      arrayView1d< real32 const > const delta    = elementSubRegion.getField< fields::Delta >();
      arrayView1d< real32 const > const vti_f    = elementSubRegion.getField< fields::F >();

      /// Partial gradient if gradient as to be computed
      arrayView1d< real32 > grad = elementSubRegion.getField< fields::PartialGradient >();
      grad.zero();

      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        acousticVTIWaveEquationSEMKernels::MassMatrixKernel< FE_TYPE > kernelM( finiteElement );

        kernelM.template launch< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                               X,
                                                               elemsToNodes,
                                                               velocity,
                                                               mass );

        acousticVTIWaveEquationSEMKernels::DampingMatrixKernel< FE_TYPE > kernelD( finiteElement );

        kernelD.template launch< EXEC_POLICY, ATOMIC_POLICY >( faceManager.size(),
                                                               X,
                                                               facesToElements,
                                                               facesToNodes,
                                                               facesDomainBoundaryIndicator,
                                                               freeSurfaceFaceIndicator,
                                                               lateralSurfaceFaceIndicator,
                                                               bottomSurfaceFaceIndicator,
                                                               velocity,
                                                               epsilon,
                                                               delta,
                                                               vti_f,
                                                               damping, //DEBUG    
                                                               damping_p,
                                                               damping_q,
                                                               damping_pq,
                                                               damping_qp );
      } );
    } );
  } );

}

void AcousticVTIWaveEquationSEM::precomputeSurfaceFieldIndicator( real64 time, DomainPartition & domain )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  FaceManager & faceManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getFaceManager();
  NodeManager & nodeManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getNodeManager();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  /// array of indicators: 1 if a face is on on lateral surface; 0 otherwise
  arrayView1d< localIndex > const lateralSurfaceFaceIndicator = faceManager.getField< fields::LateralSurfaceFaceIndicator >();
  /// array of indicators: 1 if a node is on on lateral surface; 0 otherwise
  arrayView1d< localIndex > const lateralSurfaceNodeIndicator = nodeManager.getField< fields::LateralSurfaceNodeIndicator >();

  /// array of indicators: 1 if a face is on on bottom surface; 0 otherwise
  arrayView1d< localIndex > const bottomSurfaceFaceIndicator = faceManager.getField< fields::BottomSurfaceFaceIndicator >();
  /// array of indicators: 1 if a node is on on bottom surface; 0 otherwise
  arrayView1d< localIndex > const bottomSurfaceNodeIndicator = nodeManager.getField< fields::BottomSurfaceNodeIndicator >();

  // Lateral surfaces
  fsManager.apply< FaceManager >( time,
                                  domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ),
                                  string( "LateralSurface" ),
                                  [&]( FieldSpecificationBase const & bc,
                                       string const &,
                                       SortedArrayView< localIndex const > const & targetSet,
                                       FaceManager &,
                                       string const & )
  {
    string const & functionName = bc.getFunctionName();

    if( functionName.empty() || functionManager.getGroup< FunctionBase >( functionName ).isFunctionOfTime() == 2 )
    {
//      real64 const value = bc.getScale();

      for( localIndex i = 0; i < targetSet.size(); ++i )
      {
        localIndex const kf = targetSet[ i ];
        lateralSurfaceFaceIndicator[kf] = 1;

        localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );
        for( localIndex a=0; a < numNodes; ++a )
        {
          localIndex const dof = faceToNodeMap( kf, a );
          lateralSurfaceNodeIndicator[dof] = 1;
        }
      }
    }
    else
    {
      GEOS_ERROR( "This option is not supported yet" );
    }
  } );

  // For the Bottom surface
  fsManager.apply< FaceManager >( time,
                                  domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ),
                                  string( "BottomSurface" ),
                                  [&]( FieldSpecificationBase const & bc,
                                       string const &,
                                       SortedArrayView< localIndex const > const & targetSet,
                                       FaceManager &,
                                       string const & )
  {
    string const & functionName = bc.getFunctionName();

    if( functionName.empty() || functionManager.getGroup< FunctionBase >( functionName ).isFunctionOfTime() == 2 )
    {
//      real64 const value = bc.getScale();

      for( localIndex i = 0; i < targetSet.size(); ++i )
      {
        localIndex const kf = targetSet[ i ];
        bottomSurfaceFaceIndicator[kf] = 1;

        localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );
        for( localIndex a=0; a < numNodes; ++a )
        {
          localIndex const dof = faceToNodeMap( kf, a );
          bottomSurfaceNodeIndicator[dof] = 1;
        }
      }
    }
    else
    {
      GEOS_ERROR( "This option is not supported yet" );
    }
  } );
}

void AcousticVTIWaveEquationSEM::applyFreeSurfaceBC( real64 time, DomainPartition & domain )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  FaceManager & faceManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getFaceManager();
  NodeManager & nodeManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getNodeManager();

  arrayView1d< real32 > const p_nm1 = nodeManager.getField< fields::Pressure_p_nm1 >();
  arrayView1d< real32 > const p_n = nodeManager.getField< fields::Pressure_p_n >();
  arrayView1d< real32 > const p_np1 = nodeManager.getField< fields::Pressure_p_np1 >();

  arrayView1d< real32 > const q_nm1 = nodeManager.getField< fields::Pressure_q_nm1 >();
  arrayView1d< real32 > const q_n = nodeManager.getField< fields::Pressure_q_n >();
  arrayView1d< real32 > const q_np1 = nodeManager.getField< fields::Pressure_q_np1 >();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  /// array of indicators: 1 if a face is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceFaceIndicator = faceManager.getField< fields::FreeSurfaceFaceIndicator >();

  /// array of indicators: 1 if a node is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceNodeIndicator = nodeManager.getField< fields::FreeSurfaceNodeIndicator >();

  //freeSurfaceFaceIndicator.zero();
  //freeSurfaceNodeIndicator.zero();

  fsManager.apply< FaceManager >( time,
                                  domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ),
                                  string( "FreeSurface" ),
                                  [&]( FieldSpecificationBase const & bc,
                                       string const &,
                                       SortedArrayView< localIndex const > const & targetSet,
                                       FaceManager &,
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

          p_np1[dof] = value;
          p_n[dof]   = value;
          p_nm1[dof] = value;

          q_np1[dof] = value;
          q_n[dof]   = value;
          q_nm1[dof] = value;

        }
      }
    }
    else
    {
      GEOS_ERROR( "This option is not supported yet" );
    }
  } );
}


/**
 * Checks if a directory exists.
 *
 * @param dirName Directory name to check existence of.
 * @return true is dirName exists and is a directory.
 */
bool AcousticVTIWaveEquationSEM::dirExists( const std::string & dirName )
{
  struct stat buffer;
  return stat( dirName.c_str(), &buffer ) == 0;
}

real64 AcousticVTIWaveEquationSEM::explicitStepForward( real64 const & time_n,
                                                     real64 const & dt,
                                                     integer cycleNumber,
                                                     DomainPartition & domain,
                                                     bool computeGradient )
{
  real64 dtOut = explicitStepInternal( time_n, dt, cycleNumber, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                  [&] ( string const &,
                                        MeshLevel & mesh,
                                        arrayView1d< string const > const & GEOS_UNUSED_PARAM ( regionNames ) )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< real32 > const p_nm1 = nodeManager.getField< fields::Pressure_p_nm1 >();
    arrayView1d< real32 > const p_n = nodeManager.getField< fields::Pressure_p_n >();
    arrayView1d< real32 > const p_np1 = nodeManager.getField< fields::Pressure_p_np1 >();

    arrayView1d< real32 > const q_nm1 = nodeManager.getField< fields::Pressure_q_nm1 >();
    arrayView1d< real32 > const q_n = nodeManager.getField< fields::Pressure_q_n >();
    arrayView1d< real32 > const q_np1 = nodeManager.getField< fields::Pressure_q_np1 >();

    if( computeGradient )
    {

      arrayView1d< real32 > const p_dt2 = nodeManager.getField< fields::PressureDoubleDerivative >();

      if( NULL == std::getenv( "DISABLE_LIFO" ) )
      {
        m_lifo->pushWait();
      }
      forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const nodeIdx )
      {
        p_dt2[nodeIdx] = (p_np1[nodeIdx] - 2*p_n[nodeIdx] + p_nm1[nodeIdx])/(dt*dt);
      } );

      if( NULL == std::getenv( "DISABLE_LIFO" ) )
      {
        // Need to tell LvArray data is on GPU to avoir HtoD copy
        p_dt2.move( MemorySpace::cuda, false );
        m_lifo->pushAsync( p_dt2 );
      }
      else
      {
        GEOS_MARK_SCOPE ( DirectWrite );
        p_dt2.move( MemorySpace::host, false );
        int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
        std::string fileName = GEOS_FMT( "lifo/rank_{:05}/pressuredt2_{:06}_{:08}.dat", rank, m_shotIndex, cycleNumber );
        int lastDirSeparator = fileName.find_last_of( "/\\" );
        std::string dirName = fileName.substr( 0, lastDirSeparator );
        if( string::npos != (size_t)lastDirSeparator && !dirExists( dirName ))
          makeDirsForPath( dirName );

        //std::string fileName = GEOS_FMT( "pressuredt2_{:06}_{:08}_{:04}.dat", m_shotIndex, cycleNumber, rank );
        //const int fileDesc = open( fileName.c_str(), O_CREAT | O_WRONLY | O_DIRECT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP |
        // S_IROTH | S_IWOTH );

        std::ofstream wf( fileName, std::ios::out | std::ios::binary );
        GEOS_THROW_IF( !wf,
                        "Could not open file "<< fileName << " for writting",
                        InputError );
        wf.write( (char *)&p_dt2[0], p_dt2.size()*sizeof( real32 ) );
        wf.close( );
        GEOS_THROW_IF( !wf.good(),
                        "An error occured while writting "<< fileName,
                        InputError );
      }

    }

    forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      p_nm1[a] = p_n[a];
      p_n[a]   = p_np1[a];

      q_nm1[a] = q_n[a];
      q_n[a]   = q_np1[a];

    } );
  } );

  return dtOut;
}


real64 AcousticVTIWaveEquationSEM::explicitStepBackward( real64 const & time_n,
                                                      real64 const & dt,
                                                      integer cycleNumber,
                                                      DomainPartition & domain,
                                                      bool computeGradient )
{
  real64 dtOut = explicitStepInternal( time_n, dt, cycleNumber, domain );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                  [&] ( string const &,
                                        MeshLevel & mesh,
                                        arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< real32 const > const mass = nodeManager.getField< fields::MassVector >();

    arrayView1d< real32 > const p_nm1 = nodeManager.getField< fields::Pressure_p_nm1 >();
    arrayView1d< real32 > const p_n = nodeManager.getField< fields::Pressure_p_n >();
    arrayView1d< real32 > const p_np1 = nodeManager.getField< fields::Pressure_p_np1 >();

    arrayView1d< real32 > const q_nm1 = nodeManager.getField< fields::Pressure_q_nm1 >();
    arrayView1d< real32 > const q_n = nodeManager.getField< fields::Pressure_q_n >();
    arrayView1d< real32 > const q_np1 = nodeManager.getField< fields::Pressure_q_np1 >();

    if( computeGradient )
    {
      ElementRegionManager & elemManager = mesh.getElemManager();

      arrayView1d< real32 > const p_dt2 = nodeManager.getField< fields::PressureDoubleDerivative >();

      if( NULL == std::getenv( "DISABLE_LIFO" ) )
      {
        m_lifo->pop( p_dt2 );
      }
      else
      {
        GEOS_MARK_SCOPE ( DirectRead );

        int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
        std::string fileName = GEOS_FMT( "lifo/rank_{:05}/pressuredt2_{:06}_{:08}.dat", rank, m_shotIndex, cycleNumber );
        std::ifstream wf( fileName, std::ios::in | std::ios::binary );
        GEOS_THROW_IF( !wf,
                        "Could not open file "<< fileName << " for reading",
                        InputError );
        //std::string fileName = GEOS_FMT( "pressuredt2_{:06}_{:08}_{:04}.dat", m_shotIndex, cycleNumber, rank );
        //const int fileDesc = open( fileName.c_str(), O_RDONLY | O_DIRECT );
        //GEOS_ERROR_IF( fileDesc == -1,
        //                "Could not open file "<< fileName << " for reading: " << strerror( errno ) );
        // maybe better with registerTouch()
        p_dt2.move( MemorySpace::host, true );
        wf.read( (char *)&p_dt2[0], p_dt2.size()*sizeof( real32 ) );
        wf.close( );
        remove( fileName.c_str() );
      }
      elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                  CellElementSubRegion & elementSubRegion )
      {
        arrayView1d< real32 const > const velocity = elementSubRegion.getField< fields::MediumVelocity >();
        arrayView1d< real32 > grad = elementSubRegion.getField< fields::PartialGradient >();
        arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();
        constexpr localIndex numNodesPerElem = 8;

        GEOS_MARK_SCOPE ( updatePartialGradient );
        forAll< EXEC_POLICY >( elementSubRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const eltIdx )
        {
          for( localIndex i = 0; i < numNodesPerElem; ++i )
          {
            localIndex nodeIdx = elemsToNodes[eltIdx][i];
            grad[eltIdx] += (-2/velocity[eltIdx]) * mass[nodeIdx]/8.0 * (p_dt2[nodeIdx] * p_n[nodeIdx]);
          }
        } );
      } );
    }

    forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      p_nm1[a] = p_n[a];
      p_n[a]   = p_np1[a];

      q_nm1[a] = q_n[a];
      q_n[a]   = q_np1[a];

    } );
  } );

  return dtOut;
}

real64 AcousticVTIWaveEquationSEM::explicitStepInternal( real64 const & time_n,
                                                      real64 const & dt,
                                                      integer cycleNumber,
                                                      DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  GEOS_LOG_RANK_0_IF( dt < epsilonLoc, "Warning! Value for dt: " << dt << "s is smaller than local threshold: " << epsilonLoc );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                  [&] ( string const &,
                                        MeshLevel & mesh,
                                        arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< real32 const > const mass = nodeManager.getField< fields::MassVector >();
    arrayView1d< real32 const > const damping_p = nodeManager.getField< fields::DampingVector_p >();
    arrayView1d< real32 const > const damping_q = nodeManager.getField< fields::DampingVector_q >();
    arrayView1d< real32 const > const damping_pq = nodeManager.getField< fields::DampingVector_pq >();
    arrayView1d< real32 const > const damping_qp = nodeManager.getField< fields::DampingVector_qp >();

    arrayView1d< real32 > const p_nm1 = nodeManager.getField< fields::Pressure_p_nm1 >();
    arrayView1d< real32 > const p_n = nodeManager.getField< fields::Pressure_p_n >();
    arrayView1d< real32 > const p_np1 = nodeManager.getField< fields::Pressure_p_np1 >();

    arrayView1d< real32 > const q_nm1 = nodeManager.getField< fields::Pressure_q_nm1 >();
    arrayView1d< real32 > const q_n = nodeManager.getField< fields::Pressure_q_n >();
    arrayView1d< real32 > const q_np1 = nodeManager.getField< fields::Pressure_q_np1 >();

    arrayView1d< localIndex const > const freeSurfaceNodeIndicator = nodeManager.getField< fields::FreeSurfaceNodeIndicator >();
    arrayView1d< localIndex const > const lateralSurfaceNodeIndicator = nodeManager.getField< fields::LateralSurfaceNodeIndicator >();
    arrayView1d< localIndex const > const bottomSurfaceNodeIndicator = nodeManager.getField< fields::BottomSurfaceNodeIndicator >();
    arrayView1d< real32 > const stiffnessVector_p = nodeManager.getField< fields::StiffnessVector_p >();
    arrayView1d< real32 > const stiffnessVector_q = nodeManager.getField< fields::StiffnessVector_q >();
    arrayView1d< real32 > const rhs = nodeManager.getField< fields::ForcingRHS >();

    //DEbug
    arrayView1d< real32 > const stiffnessVector = nodeManager.getField< fields::StiffnessVector >();
    arrayView1d< real32 > const damping = nodeManager.getField< fields::DampingVector >();

    auto kernelFactory = acousticVTIWaveEquationSEMKernels::ExplicitAcousticSEMFactory( dt );

    finiteElement::
      regionBasedKernelApplication< EXEC_POLICY,
                                    constitutive::NullModel,
                                    CellElementSubRegion >( mesh,
                                                            regionNames,
                                                            getDiscretizationName(),
                                                            "",
                                                            kernelFactory );

    addSourceToRightHandSide( cycleNumber, rhs );

    /// calculate your time integrators
    real64 const dt2 = dt*dt;

//    if( !usePML )
  //  {
      GEOS_MARK_SCOPE ( updateP );
      forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        //Debug
        /*
        if( freeSurfaceNodeIndicator[a] != 1 )
        {
          p_np1[a] = p_n[a];
          p_np1[a] *= 2.0*mass[a];
          p_np1[a] -= (mass[a]-0.5*dt*damping[a])*p_nm1[a];
          p_np1[a] += dt2*(rhs[a]-stiffnessVector[a]);

          q_np1[a] = 1;

          if(lateralSurfaceNodeIndicator[a] != 1 && bottomSurfaceNodeIndicator[a] != 1)
          {
            // Interior node, no boundary terms
            p_np1[a] /= mass[a];
          }
          else
          {
            p_np1[a] /= mass[a]+0.5*dt*damping[a];
          }
        }*/
        // Good ?
        if( freeSurfaceNodeIndicator[a] != 1 )
        {
          p_np1[a] = 2.0*mass[a]*p_n[a]/dt2;
          p_np1[a] -= mass[a]*p_nm1[a]/dt2;
          p_np1[a] += stiffnessVector_p[a];
          p_np1[a] += rhs[a];

          q_np1[a] = 2.0*mass[a]*q_n[a]/dt2;
          q_np1[a] -= mass[a]*q_nm1[a]/dt2;
          q_np1[a] += stiffnessVector_q[a];
          q_np1[a] += rhs[a];

          if(lateralSurfaceNodeIndicator[a] != 1 && bottomSurfaceNodeIndicator[a] != 1)
          {
            // Interior node, no boundary terms
            p_np1[a] /= mass[a]/dt2;
            q_np1[a] /= mass[a]/dt2;
          }
          else
          {
            // Boundary node
            // Hand-made Inversion of 2x2 matrix
            real32 coef_pp = mass[a]/dt2;
            coef_pp += damping_p[a]/dt/2;
            real32 coef_pq = damping_pq[a]/dt/2;

            real32 coef_qq = mass[a]/dt2;
            coef_qq += damping_q[a]/2/dt;
            real32 coef_qp = damping_qp[a]/dt/2;
            
            real32 det_pq = 1/(coef_pp * coef_qq - coef_pq*coef_qp);
            
            real32 aux_p_np1 = p_np1[a];
            p_np1[a] = det_pq*(coef_qq*p_np1[a] - coef_pq*q_np1[a]);
            q_np1[a] = det_pq*(coef_pp*q_np1[a] - coef_qp*aux_p_np1);
          }
        }
      } );
//    }

    /// synchronize pressure fields
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addFields( FieldLocation::Node, { fields::Pressure_p_np1::key() } );
    fieldsToBeSync.addFields( FieldLocation::Node, { fields::Pressure_q_np1::key() } );

    CommunicationTools & syncFields = CommunicationTools::getInstance();
    syncFields.synchronizeFields( fieldsToBeSync,
                                  mesh,
                                  domain.getNeighbors(),
                                  true );

    // compute the seismic traces since last step.
    arrayView2d< real32 > const pReceivers   = m_pressureNp1AtReceivers.toView();

    computeAllSeismoTraces( time_n, dt, p_np1, p_n, pReceivers );

    /// prepare next step
    forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      stiffnessVector[a] = 0.0; //Debug
      stiffnessVector_p[a] = 0.0;
      stiffnessVector_q[a] = 0.0;
      rhs[a] = 0.0;
    } );

  } );

  return dt;
}

void AcousticVTIWaveEquationSEM::cleanup( real64 const time_n,
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
    arrayView1d< real32 const > const p_n = nodeManager.getField< fields::Pressure_p_n >();
    arrayView1d< real32 const > const p_np1 = nodeManager.getField< fields::Pressure_p_np1 >();
    arrayView2d< real32 > const pReceivers   = m_pressureNp1AtReceivers.toView();
    computeAllSeismoTraces( time_n, 0, p_np1, p_n, pReceivers );
  } );
}

void AcousticVTIWaveEquationSEM::computeAllSeismoTraces( real64 const time_n,
                                                      real64 const dt,
                                                      arrayView1d< real32 const > const var_np1,
                                                      arrayView1d< real32 const > const var_n,
                                                      arrayView2d< real32 > varAtReceivers )
{

  /*
   * In forward case we compute seismo if time_n + dt  is the first time
   * step after the timeSeismo to write.
   *
   *  time_n        timeSeismo    time_n + dt
   *   ---|--------------|-------------|
   *
   * In backward (time_n goes decreasing) case we compute seismo if
   * time_n is the last time step before the timeSeismo to write.
   *
   *  time_n - dt    timeSeismo    time_n
   *   ---|--------------|-------------|
   */
  for( real64 timeSeismo;
       (m_forward)?((timeSeismo = m_dtSeismoTrace*m_indexSeismoTrace) <= (time_n + dt + epsilonLoc) && m_indexSeismoTrace < m_nsamplesSeismoTrace):
       ((timeSeismo = m_dtSeismoTrace*(m_nsamplesSeismoTrace-m_indexSeismoTrace-1)) >= (time_n - dt -  epsilonLoc) && m_indexSeismoTrace < m_nsamplesSeismoTrace);
       m_indexSeismoTrace++ )
  {
    WaveSolverUtils::computeSeismoTrace( time_n, (m_forward)?dt:-dt, timeSeismo, (m_forward)?m_indexSeismoTrace:(m_nsamplesSeismoTrace-m_indexSeismoTrace-1), m_receiverNodeIds, m_receiverConstants,
                                         m_receiverIsLocal,
                                         m_nsamplesSeismoTrace, m_outputSeismoTrace, var_np1, var_n, varAtReceivers );
  }
}


REGISTER_CATALOG_ENTRY( SolverBase, AcousticVTIWaveEquationSEM, string const &, dataRepository::Group * const )

} /* namespace geos */
