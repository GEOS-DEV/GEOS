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
#include "WaveSolverBaseFields.hpp"

namespace geos
{

using namespace dataRepository;

AcousticVTIWaveEquationSEM::AcousticVTIWaveEquationSEM( const std::string & name,
                                                        Group * const parent ):
  WaveSolverBase( name,
                  parent )
{

  registerWrapper( viewKeyStruct::pressureNp1AtReceiversString(), &m_pressureNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Pressure value at each receiver for each timestep" );
}

void AcousticVTIWaveEquationSEM::initializePreSubGroups()
{
  WaveSolverBase::initializePreSubGroups();

}

void AcousticVTIWaveEquationSEM::registerDataOnMesh( Group & meshBodies )
{
  WaveSolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    nodeManager.registerField< fields::wavesolverfields::Pressure_p_nm1,
                               fields::wavesolverfields::Pressure_p_n,
                               fields::wavesolverfields::Pressure_p_np1,
                               fields::wavesolverfields::Pressure_q_nm1,
                               fields::wavesolverfields::Pressure_q_n,
                               fields::wavesolverfields::Pressure_q_np1,
                               fields::wavesolverfields::ForcingRHS,
                               fields::wavesolverfields::MassVector,
                               fields::wavesolverfields::DampingVector_p,
                               fields::wavesolverfields::DampingVector_pq,
                               fields::wavesolverfields::DampingVector_q,
                               fields::wavesolverfields::DampingVector_qp,
                               fields::wavesolverfields::StiffnessVector_p,
                               fields::wavesolverfields::StiffnessVector_q,
                               fields::wavesolverfields::FreeSurfaceNodeIndicator,
                               fields::wavesolverfields::LateralSurfaceNodeIndicator,
                               fields::wavesolverfields::BottomSurfaceNodeIndicator >( getName() );


    FaceManager & faceManager = mesh.getFaceManager();
    faceManager.registerField< fields::wavesolverfields::FreeSurfaceFaceIndicator >( getName() );
    faceManager.registerField< fields::wavesolverfields::LateralSurfaceFaceIndicator >( getName() );
    faceManager.registerField< fields::wavesolverfields::BottomSurfaceFaceIndicator >( getName() );

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerField< fields::wavesolverfields::MediumDelta >( getName() );
      subRegion.registerField< fields::wavesolverfields::MediumEpsilon >( getName() );
      subRegion.registerField< fields::wavesolverfields::MediumSigma >( getName() );
      subRegion.registerField< fields::wavesolverfields::MediumVelocity >( getName() );
      subRegion.registerField< fields::wavesolverfields::MediumDensity >( getName() );
    } );
  } );
}


void AcousticVTIWaveEquationSEM::postProcessInput()
{

  WaveSolverBase::postProcessInput();

  localIndex const numReceiversGlobal = m_receiverCoordinates.size( 0 );

  m_pressureNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal + 1 );
}

void AcousticVTIWaveEquationSEM::precomputeSourceAndReceiverTerm( MeshLevel & mesh,
                                                                  arrayView1d< string const > const & regionNames )
{
  NodeManager const & nodeManager = mesh.getNodeManager();
  FaceManager const & faceManager = mesh.getFaceManager();

  arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const X32 =
    nodeManager.getField< fields::referencePosition32 >().toViewConst();

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

      localIndex const numFacesPerElem = elementSubRegion.numFacesPerElement();

      acousticVTIWaveEquationSEMKernels::
        PrecomputeSourceAndReceiverKernel::
        launch< EXEC_POLICY, FE_TYPE >
        ( elementSubRegion.size(),
        numFacesPerElem,
        X32,
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
        m_timeSourceFrequency,
        m_timeSourceDelay,
        m_rickerOrder );
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

  DomainPartition & domain = getGroupByPath< DomainPartition >( "/Problem/domain" );

  applyFreeSurfaceBC( 0.0, domain );
  precomputeSurfaceFieldIndicator( domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    precomputeSourceAndReceiverTerm( mesh, regionNames );

    NodeManager & nodeManager = mesh.getNodeManager();
    FaceManager & faceManager = mesh.getFaceManager();

    /// get the array of indicators: 1 if the face is on the boundary; 0 otherwise
    arrayView1d< integer > const & facesDomainBoundaryIndicator = faceManager.getDomainBoundaryIndicator();
    arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords =
      nodeManager.getField< fields::referencePosition32 >().toViewConst();

    /// get face to node map
    ArrayOfArraysView< localIndex const > const facesToNodes = faceManager.nodeList().toViewConst();

    // mass matrix to be computed in this function
    arrayView1d< real32 > const mass = nodeManager.getField< fields::wavesolverfields::MassVector >();
    mass.zero();
    /// damping matrices to be computed for each dof in the boundary of the mesh
    arrayView1d< real32 > const damping_p  = nodeManager.getField< fields::wavesolverfields::DampingVector_p >();
    arrayView1d< real32 > const damping_pq = nodeManager.getField< fields::wavesolverfields::DampingVector_pq >();
    arrayView1d< real32 > const damping_q  = nodeManager.getField< fields::wavesolverfields::DampingVector_q >();
    arrayView1d< real32 > const damping_qp = nodeManager.getField< fields::wavesolverfields::DampingVector_qp >();
    damping_p.zero();
    damping_pq.zero();
    damping_q.zero();
    damping_qp.zero();

    /// get array of indicators: 1 if face is on the free surface; 0 otherwise
    arrayView1d< localIndex const > const freeSurfaceFaceIndicator = faceManager.getField< fields::wavesolverfields::FreeSurfaceFaceIndicator >();
    arrayView1d< localIndex const > const lateralSurfaceFaceIndicator = faceManager.getField< fields::wavesolverfields::LateralSurfaceFaceIndicator >();
    arrayView1d< localIndex const > const bottomSurfaceFaceIndicator = faceManager.getField< fields::wavesolverfields::BottomSurfaceFaceIndicator >();

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          CellElementSubRegion & elementSubRegion )
    {

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = elementSubRegion.nodeList();
      arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();
      arrayView1d< real32 const > const velocity = elementSubRegion.getField< fields::wavesolverfields::MediumVelocity >();
      arrayView1d< real32 const > const density  = elementSubRegion.getField< fields::wavesolverfields::MediumDensity >();
      arrayView1d< real32 const > const vti_epsilon  = elementSubRegion.getField< fields::wavesolverfields::MediumEpsilon >();
      arrayView1d< real32 const > const vti_delta    = elementSubRegion.getField< fields::wavesolverfields::MediumDelta >();
      arrayView1d< real32 const > const vti_sigma    = elementSubRegion.getField< fields::wavesolverfields::MediumSigma >();

      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        acousticVTIWaveEquationSEMKernels::MassMatrixKernel< FE_TYPE > kernelM( finiteElement );

        kernelM.template launch< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                               nodeCoords,
                                                               elemsToNodes,
                                                               velocity,
                                                               density,
                                                               mass );

        acousticVTIWaveEquationSEMKernels::DampingMatrixKernel< FE_TYPE > kernelD( finiteElement );

        kernelD.template launch< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                               nodeCoords,
                                                               elemsToFaces,
                                                               facesToNodes,
                                                               facesDomainBoundaryIndicator,
                                                               freeSurfaceFaceIndicator,
                                                               lateralSurfaceFaceIndicator,
                                                               bottomSurfaceFaceIndicator,
                                                               velocity,
                                                               density,
                                                               vti_epsilon,
                                                               vti_delta,
                                                               vti_sigma,
                                                               damping_p,
                                                               damping_q,
                                                               damping_pq,
                                                               damping_qp );
      } );
    } );
  } );

  WaveSolverUtils::initTrace( "seismoTraceReceiver", getName(), m_outputSeismoTrace, m_receiverConstants.size( 0 ), m_receiverIsLocal );
}

void AcousticVTIWaveEquationSEM::precomputeSurfaceFieldIndicator( DomainPartition & domain )
{
  real64 const time = 0.0;
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  FaceManager & faceManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getFaceManager();
  NodeManager & nodeManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getNodeManager();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  /// array of indicators: 1 if a face is on on lateral surface; 0 otherwise
  arrayView1d< localIndex > const lateralSurfaceFaceIndicator = faceManager.getField< fields::wavesolverfields::LateralSurfaceFaceIndicator >();
  /// array of indicators: 1 if a node is on on lateral surface; 0 otherwise
  arrayView1d< localIndex > const lateralSurfaceNodeIndicator = nodeManager.getField< fields::wavesolverfields::LateralSurfaceNodeIndicator >();

  /// array of indicators: 1 if a face is on on bottom surface; 0 otherwise
  arrayView1d< localIndex > const bottomSurfaceFaceIndicator = faceManager.getField< fields::wavesolverfields::BottomSurfaceFaceIndicator >();
  /// array of indicators: 1 if a node is on on bottom surface; 0 otherwise
  arrayView1d< localIndex > const bottomSurfaceNodeIndicator = nodeManager.getField< fields::wavesolverfields::BottomSurfaceNodeIndicator >();

  // Lateral surfaces
  fsManager.apply< FaceManager >( time,
                                  domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ),
                                  viewKeyStruct::lateralSurfaceString(),
                                  [&]( FieldSpecificationBase const & bc,
                                       string const &,
                                       SortedArrayView< localIndex const > const & targetSet,
                                       FaceManager &,
                                       string const & )
  {
    string const & functionName = bc.getFunctionName();

    if( functionName.empty() || functionManager.getGroup< FunctionBase >( functionName ).isFunctionOfTime() == 2 )
    {
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
                                  viewKeyStruct::bottomSurfaceString(),
                                  [&]( FieldSpecificationBase const & bc,
                                       string const &,
                                       SortedArrayView< localIndex const > const & targetSet,
                                       FaceManager &,
                                       string const & )
  {
    string const & functionName = bc.getFunctionName();

    if( functionName.empty() || functionManager.getGroup< FunctionBase >( functionName ).isFunctionOfTime() == 2 )
    {
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

  arrayView1d< real32 > const p_nm1 = nodeManager.getField< fields::wavesolverfields::Pressure_p_nm1 >();
  arrayView1d< real32 > const p_n = nodeManager.getField< fields::wavesolverfields::Pressure_p_n >();
  arrayView1d< real32 > const p_np1 = nodeManager.getField< fields::wavesolverfields::Pressure_p_np1 >();

  arrayView1d< real32 > const q_nm1 = nodeManager.getField< fields::wavesolverfields::Pressure_q_nm1 >();
  arrayView1d< real32 > const q_n = nodeManager.getField< fields::wavesolverfields::Pressure_q_n >();
  arrayView1d< real32 > const q_np1 = nodeManager.getField< fields::wavesolverfields::Pressure_q_np1 >();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  /// array of indicators: 1 if a face is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceFaceIndicator = faceManager.getField< fields::wavesolverfields::FreeSurfaceFaceIndicator >();

  /// array of indicators: 1 if a node is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceNodeIndicator = nodeManager.getField< fields::wavesolverfields::FreeSurfaceNodeIndicator >();

  fsManager.apply< FaceManager >( time,
                                  domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ),
                                  WaveSolverBase::viewKeyStruct::freeSurfaceString(),
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

    arrayView1d< real32 > const p_nm1 = nodeManager.getField< fields::wavesolverfields::Pressure_p_nm1 >();
    arrayView1d< real32 > const p_n = nodeManager.getField< fields::wavesolverfields::Pressure_p_n >();
    arrayView1d< real32 > const p_np1 = nodeManager.getField< fields::wavesolverfields::Pressure_p_np1 >();

    arrayView1d< real32 > const q_nm1 = nodeManager.getField< fields::wavesolverfields::Pressure_q_nm1 >();
    arrayView1d< real32 > const q_n = nodeManager.getField< fields::wavesolverfields::Pressure_q_n >();
    arrayView1d< real32 > const q_np1 = nodeManager.getField< fields::wavesolverfields::Pressure_q_np1 >();

    if( computeGradient )
    {
      GEOS_ERROR( "This option is not supported yet" );
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

void AcousticVTIWaveEquationSEM::initializePML()
{
  GEOS_ERROR( "This option is not supported yet" );
  return;
}

void AcousticVTIWaveEquationSEM::applyPML( real64 const GEOS_UNUSED_PARAM( time ),
                                           DomainPartition & GEOS_UNUSED_PARAM( domain ))
{
  GEOS_ERROR( "This option is not supported yet" );
  return;
}

real64 AcousticVTIWaveEquationSEM::explicitStepBackward( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                                         real64 const & GEOS_UNUSED_PARAM( dt ),
                                                         integer GEOS_UNUSED_PARAM( cycleNumber ),
                                                         DomainPartition & GEOS_UNUSED_PARAM( domain ),
                                                         bool GEOS_UNUSED_PARAM( computeGradient ) )
{
  GEOS_ERROR( "This option is not supported yet" );
  return -1;
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

    arrayView1d< real32 const > const mass = nodeManager.getField< fields::wavesolverfields::MassVector >();
    arrayView1d< real32 const > const damping_p = nodeManager.getField< fields::wavesolverfields::DampingVector_p >();
    arrayView1d< real32 const > const damping_q = nodeManager.getField< fields::wavesolverfields::DampingVector_q >();
    arrayView1d< real32 const > const damping_pq = nodeManager.getField< fields::wavesolverfields::DampingVector_pq >();
    arrayView1d< real32 const > const damping_qp = nodeManager.getField< fields::wavesolverfields::DampingVector_qp >();

    arrayView1d< real32 > const p_nm1 = nodeManager.getField< fields::wavesolverfields::Pressure_p_nm1 >();
    arrayView1d< real32 > const p_n = nodeManager.getField< fields::wavesolverfields::Pressure_p_n >();
    arrayView1d< real32 > const p_np1 = nodeManager.getField< fields::wavesolverfields::Pressure_p_np1 >();

    arrayView1d< real32 > const q_nm1 = nodeManager.getField< fields::wavesolverfields::Pressure_q_nm1 >();
    arrayView1d< real32 > const q_n = nodeManager.getField< fields::wavesolverfields::Pressure_q_n >();
    arrayView1d< real32 > const q_np1 = nodeManager.getField< fields::wavesolverfields::Pressure_q_np1 >();

    arrayView1d< localIndex const > const freeSurfaceNodeIndicator = nodeManager.getField< fields::wavesolverfields::FreeSurfaceNodeIndicator >();
    arrayView1d< localIndex const > const lateralSurfaceNodeIndicator = nodeManager.getField< fields::wavesolverfields::LateralSurfaceNodeIndicator >();
    arrayView1d< localIndex const > const bottomSurfaceNodeIndicator = nodeManager.getField< fields::wavesolverfields::BottomSurfaceNodeIndicator >();
    arrayView1d< real32 > const stiffnessVector_p = nodeManager.getField< fields::wavesolverfields::StiffnessVector_p >();
    arrayView1d< real32 > const stiffnessVector_q = nodeManager.getField< fields::wavesolverfields::StiffnessVector_q >();
    arrayView1d< real32 > const rhs = nodeManager.getField< fields::wavesolverfields::ForcingRHS >();

    auto kernelFactory = acousticVTIWaveEquationSEMKernels::ExplicitAcousticVTISEMFactory( dt );

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

    GEOS_MARK_SCOPE ( updateP );
    forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
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

        if( lateralSurfaceNodeIndicator[a] != 1 && bottomSurfaceNodeIndicator[a] != 1 )
        {
          // Interior node, no boundary terms
          p_np1[a] /= mass[a]/dt2;
          q_np1[a] /= mass[a]/dt2;
        }
        else
        {
          // Boundary node
          p_np1[a] += damping_p[a]*p_nm1[a]/dt/2;
          p_np1[a] += damping_pq[a]*q_nm1[a]/dt/2;

          q_np1[a] += damping_q[a]*q_nm1[a]/dt/2;
          q_np1[a] += damping_qp[a]*p_nm1[a]/dt/2;
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

    /// synchronize pressure fields
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addFields( FieldLocation::Node, { fields::wavesolverfields::Pressure_p_np1::key() } );
    fieldsToBeSync.addFields( FieldLocation::Node, { fields::wavesolverfields::Pressure_q_np1::key() } );

    CommunicationTools & syncFields = CommunicationTools::getInstance();
    syncFields.synchronizeFields( fieldsToBeSync,
                                  mesh,
                                  domain.getNeighbors(),
                                  true );

    // compute the seismic traces since last step.
    arrayView2d< real32 > const pReceivers = m_pressureNp1AtReceivers.toView();
    computeAllSeismoTraces( time_n, dt, p_np1, p_n, pReceivers );

    incrementIndexSeismoTrace( time_n );

    /// prepare next step
    forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
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
    arrayView1d< real32 const > const p_n   = nodeManager.getField< fields::wavesolverfields::Pressure_p_n >();
    arrayView1d< real32 const > const p_np1 = nodeManager.getField< fields::wavesolverfields::Pressure_p_np1 >();
    arrayView2d< real32 > const pReceivers = m_pressureNp1AtReceivers.toView();
    computeAllSeismoTraces( time_n, 0.0, p_np1, p_n, pReceivers );

    WaveSolverUtils::writeSeismoTrace( "seismoTraceReceiver", getName(), m_outputSeismoTrace, m_receiverConstants.size( 0 ),
                                       m_receiverIsLocal, m_nsamplesSeismoTrace, pReceivers );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, AcousticVTIWaveEquationSEM, string const &, dataRepository::Group * const )

} /* namespace geos */
