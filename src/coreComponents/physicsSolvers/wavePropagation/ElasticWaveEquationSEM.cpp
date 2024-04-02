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
#include "ElasticTimeSchemeSEMKernel.hpp"
#include "ElasticMatricesSEMKernel.hpp"

namespace geos
{

using namespace dataRepository;
using namespace fields;

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

  registerWrapper( viewKeyStruct::dasSignalNp1AtReceiversString(), &m_dasSignalNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "DAS signal value at each receiver for each timestep" );

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

  localIndex const numNodesPerElem = WaveSolverBase::getNumNodesPerElem();

  localIndex const numSourcesGlobal = m_sourceCoordinates.size( 0 );
  m_sourceConstantsx.resize( numSourcesGlobal, numNodesPerElem );
  m_sourceConstantsy.resize( numSourcesGlobal, numNodesPerElem );
  m_sourceConstantsz.resize( numSourcesGlobal, numNodesPerElem );

  localIndex const numReceiversGlobal = m_receiverCoordinates.size( 0 );
  integer nsamples = m_useDAS == WaveSolverUtils::DASType::none ? 1 : m_linearDASSamples;
  m_receiverConstants.resize( numReceiversGlobal, nsamples * numNodesPerElem );
  m_receiverNodeIds.resize( numReceiversGlobal, nsamples * numNodesPerElem );
}


void ElasticWaveEquationSEM::registerDataOnMesh( Group & meshBodies )
{
  WaveSolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    nodeManager.registerField< elasticfields::Displacementx_nm1,
                               elasticfields::Displacementy_nm1,
                               elasticfields::Displacementz_nm1,
                               elasticfields::Displacementx_n,
                               elasticfields::Displacementy_n,
                               elasticfields::Displacementz_n,
                               elasticfields::Displacementx_np1,
                               elasticfields::Displacementy_np1,
                               elasticfields::Displacementz_np1,
                               elasticfields::ForcingRHSx,
                               elasticfields::ForcingRHSy,
                               elasticfields::ForcingRHSz,
                               elasticfields::ElasticMassVector,
                               elasticfields::DampingVectorx,
                               elasticfields::DampingVectory,
                               elasticfields::DampingVectorz,
                               elasticfields::StiffnessVectorx,
                               elasticfields::StiffnessVectory,
                               elasticfields::StiffnessVectorz,
                               elasticfields::ElasticFreeSurfaceNodeIndicator >( getName() );

    FaceManager & faceManager = mesh.getFaceManager();
    faceManager.registerField< elasticfields::ElasticFreeSurfaceFaceIndicator >( getName() );

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerField< elasticfields::ElasticVelocityVp >( getName() );
      subRegion.registerField< elasticfields::ElasticVelocityVs >( getName() );
      subRegion.registerField< elasticfields::ElasticDensity >( getName() );
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

  localIndex const numSourcesGlobal = m_sourceCoordinates.size( 0 );
  m_sourceIsAccessible.resize( numSourcesGlobal );
  m_sourceValue.resize( nsamples, numSourcesGlobal );

  if( m_useDAS == WaveSolverUtils::DASType::none )
  {
    localIndex const numReceiversGlobal = m_receiverCoordinates.size( 0 );
    m_displacementXNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal + 1 );
    m_displacementYNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal + 1 );
    m_displacementZNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal + 1 );
    m_displacementXNp1AtReceivers.zero();
    m_displacementYNp1AtReceivers.zero();
    m_displacementZNp1AtReceivers.zero();
  }
  else
  {
    localIndex const numReceiversGlobal = m_linearDASGeometry.size( 0 );
    m_dasSignalNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal + 1 );
    m_dasSignalNp1AtReceivers.zero();
  }
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
    arrayView1d< real32 > const mass = nodeManager.getField< elasticfields::ElasticMassVector >();
    mass.zero();
    /// damping matrix to be computed for each dof in the boundary of the mesh
    arrayView1d< real32 > const dampingx = nodeManager.getField< elasticfields::DampingVectorx >();
    arrayView1d< real32 > const dampingy = nodeManager.getField< elasticfields::DampingVectory >();
    arrayView1d< real32 > const dampingz = nodeManager.getField< elasticfields::DampingVectorz >();
    dampingx.zero();
    dampingy.zero();
    dampingz.zero();

    /// get array of indicators: 1 if face is on the free surface; 0 otherwise
    arrayView1d< localIndex const > const freeSurfaceFaceIndicator    = faceManager.getField< elasticfields::ElasticFreeSurfaceFaceIndicator >();
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

      arrayView1d< real32 const > const density = elementSubRegion.getField< elasticfields::ElasticDensity >();
      arrayView1d< real32 const > const velocityVp = elementSubRegion.getField< elasticfields::ElasticVelocityVp >();
      arrayView1d< real32 const > const velocityVs = elementSubRegion.getField< elasticfields::ElasticVelocityVs >();

      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        ElasticMatricesSEM::MassMatrix< FE_TYPE > kernelM( finiteElement );
        kernelM.template computeMassMatrix< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                                          nodeCoords,
                                                                          elemsToNodes,
                                                                          density,
                                                                          mass );

        ElasticMatricesSEM::DampingMatrix< FE_TYPE > kernelD( finiteElement );
        kernelD.template computeDampingMatrix< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
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

  arrayView1d< real32 > const ux_np1 = nodeManager.getField< elasticfields::Displacementx_np1 >();
  arrayView1d< real32 > const uy_np1 = nodeManager.getField< elasticfields::Displacementy_np1 >();
  arrayView1d< real32 > const uz_np1 = nodeManager.getField< elasticfields::Displacementz_np1 >();
  arrayView1d< real32 > const ux_n   = nodeManager.getField< elasticfields::Displacementx_n >();
  arrayView1d< real32 > const uy_n   = nodeManager.getField< elasticfields::Displacementy_n >();
  arrayView1d< real32 > const uz_n   = nodeManager.getField< elasticfields::Displacementz_n >();
  arrayView1d< real32 > const ux_nm1 = nodeManager.getField< elasticfields::Displacementx_nm1 >();
  arrayView1d< real32 > const uy_nm1 = nodeManager.getField< elasticfields::Displacementy_nm1 >();
  arrayView1d< real32 > const uz_nm1 = nodeManager.getField< elasticfields::Displacementz_nm1 >();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  /// set array of indicators: 1 if a face is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceFaceIndicator = faceManager.getField< elasticfields::ElasticFreeSurfaceFaceIndicator >();

  /// set array of indicators: 1 if a node is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceNodeIndicator = nodeManager.getField< elasticfields::ElasticFreeSurfaceNodeIndicator >();


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

  arrayView1d< real32 const > const mass = nodeManager.getField< elasticfields::ElasticMassVector >();
  arrayView1d< real32 const > const dampingx = nodeManager.getField< elasticfields::DampingVectorx >();
  arrayView1d< real32 const > const dampingy = nodeManager.getField< elasticfields::DampingVectory >();
  arrayView1d< real32 const > const dampingz = nodeManager.getField< elasticfields::DampingVectorz >();
  arrayView1d< real32 > const stiffnessVectorx = nodeManager.getField< elasticfields::StiffnessVectorx >();
  arrayView1d< real32 > const stiffnessVectory = nodeManager.getField< elasticfields::StiffnessVectory >();
  arrayView1d< real32 > const stiffnessVectorz = nodeManager.getField< elasticfields::StiffnessVectorz >();

  arrayView1d< real32 > const ux_nm1 = nodeManager.getField< elasticfields::Displacementx_nm1 >();
  arrayView1d< real32 > const uy_nm1 = nodeManager.getField< elasticfields::Displacementy_nm1 >();
  arrayView1d< real32 > const uz_nm1 = nodeManager.getField< elasticfields::Displacementz_nm1 >();
  arrayView1d< real32 > const ux_n = nodeManager.getField< elasticfields::Displacementx_n >();
  arrayView1d< real32 > const uy_n = nodeManager.getField< elasticfields::Displacementy_n >();
  arrayView1d< real32 > const uz_n = nodeManager.getField< elasticfields::Displacementz_n >();
  arrayView1d< real32 > const ux_np1 = nodeManager.getField< elasticfields::Displacementx_np1 >();
  arrayView1d< real32 > const uy_np1 = nodeManager.getField< elasticfields::Displacementy_np1 >();
  arrayView1d< real32 > const uz_np1 = nodeManager.getField< elasticfields::Displacementz_np1 >();

  /// get array of indicators: 1 if node on free surface; 0 otherwise
  arrayView1d< localIndex const > const freeSurfaceNodeIndicator = nodeManager.getField< elasticfields::ElasticFreeSurfaceNodeIndicator >();

  arrayView1d< real32 > const rhsx = nodeManager.getField< elasticfields::ForcingRHSx >();
  arrayView1d< real32 > const rhsy = nodeManager.getField< elasticfields::ForcingRHSy >();
  arrayView1d< real32 > const rhsz = nodeManager.getField< elasticfields::ForcingRHSz >();

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

  SortedArrayView< localIndex const > const solverTargetNodesSet = m_solverTargetNodesSet.toViewConst();
  ElasticTimeSchemeSEM::LeapFrog( dt, ux_np1, ux_n, ux_nm1, uy_np1, uy_n, uy_nm1, uz_np1, uz_n, uz_nm1,
                                  mass, dampingx, dampingy, dampingz, stiffnessVectorx, stiffnessVectory,
                                  stiffnessVectorz, rhsx, rhsy, rhsz, freeSurfaceNodeIndicator,
                                  solverTargetNodesSet );
}

void ElasticWaveEquationSEM::synchronizeUnknowns( real64 const & time_n,
                                                  real64 const & dt,
                                                  integer const,
                                                  DomainPartition & domain,
                                                  MeshLevel & mesh,
                                                  arrayView1d< string const > const & )
{
  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView1d< real32 > const stiffnessVectorx = nodeManager.getField< elasticfields::StiffnessVectorx >();
  arrayView1d< real32 > const stiffnessVectory = nodeManager.getField< elasticfields::StiffnessVectory >();
  arrayView1d< real32 > const stiffnessVectorz = nodeManager.getField< elasticfields::StiffnessVectorz >();

  arrayView1d< real32 > const ux_nm1 = nodeManager.getField< elasticfields::Displacementx_nm1 >();
  arrayView1d< real32 > const uy_nm1 = nodeManager.getField< elasticfields::Displacementy_nm1 >();
  arrayView1d< real32 > const uz_nm1 = nodeManager.getField< elasticfields::Displacementz_nm1 >();
  arrayView1d< real32 > const ux_n   = nodeManager.getField< elasticfields::Displacementx_n >();
  arrayView1d< real32 > const uy_n   = nodeManager.getField< elasticfields::Displacementy_n >();
  arrayView1d< real32 > const uz_n   = nodeManager.getField< elasticfields::Displacementz_n >();
  arrayView1d< real32 > const ux_np1 = nodeManager.getField< elasticfields::Displacementx_np1 >();
  arrayView1d< real32 > const uy_np1 = nodeManager.getField< elasticfields::Displacementy_np1 >();
  arrayView1d< real32 > const uz_np1 = nodeManager.getField< elasticfields::Displacementz_np1 >();

  arrayView1d< real32 > const rhsx = nodeManager.getField< elasticfields::ForcingRHSx >();
  arrayView1d< real32 > const rhsy = nodeManager.getField< elasticfields::ForcingRHSy >();
  arrayView1d< real32 > const rhsz = nodeManager.getField< elasticfields::ForcingRHSz >();

  /// synchronize displacement fields
  FieldIdentifiers fieldsToBeSync;
  fieldsToBeSync.addFields( FieldLocation::Node, { elasticfields::Displacementx_np1::key(), elasticfields::Displacementy_np1::key(), elasticfields::Displacementz_np1::key() } );

  CommunicationTools & syncFields = CommunicationTools::getInstance();
  syncFields.synchronizeFields( fieldsToBeSync,
                                domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ),
                                domain.getNeighbors(),
                                true );

  // compute the seismic traces since last step.
  if( m_useDAS == WaveSolverUtils::DASType::none )
  {
    arrayView2d< real32 > const uXReceivers = m_displacementXNp1AtReceivers.toView();
    arrayView2d< real32 > const uYReceivers = m_displacementYNp1AtReceivers.toView();
    arrayView2d< real32 > const uZReceivers = m_displacementZNp1AtReceivers.toView();
    computeAllSeismoTraces( time_n, dt, ux_np1, ux_n, uXReceivers );
    computeAllSeismoTraces( time_n, dt, uy_np1, uy_n, uYReceivers );
    computeAllSeismoTraces( time_n, dt, uz_np1, uz_n, uZReceivers );
  }
  else
  {
    arrayView2d< real32 > const dasReceivers  = m_dasSignalNp1AtReceivers.toView();
    computeAllSeismoTraces( time_n, dt, ux_np1, ux_n, dasReceivers, m_linearDASVectorX.toView(), true );
    computeAllSeismoTraces( time_n, dt, uy_np1, uy_n, dasReceivers, m_linearDASVectorY.toView(), true );
    computeAllSeismoTraces( time_n, dt, uz_np1, uz_n, dasReceivers, m_linearDASVectorZ.toView(), true );
  }

  incrementIndexSeismoTrace( time_n );
}

void ElasticWaveEquationSEM::prepareNextTimestep( MeshLevel & mesh )
{
  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView1d< real32 > const ux_nm1 = nodeManager.getField< elasticfields::Displacementx_nm1 >();
  arrayView1d< real32 > const uy_nm1 = nodeManager.getField< elasticfields::Displacementy_nm1 >();
  arrayView1d< real32 > const uz_nm1 = nodeManager.getField< elasticfields::Displacementz_nm1 >();
  arrayView1d< real32 > const ux_n   = nodeManager.getField< elasticfields::Displacementx_n >();
  arrayView1d< real32 > const uy_n   = nodeManager.getField< elasticfields::Displacementy_n >();
  arrayView1d< real32 > const uz_n   = nodeManager.getField< elasticfields::Displacementz_n >();
  arrayView1d< real32 > const ux_np1 = nodeManager.getField< elasticfields::Displacementx_np1 >();
  arrayView1d< real32 > const uy_np1 = nodeManager.getField< elasticfields::Displacementy_np1 >();
  arrayView1d< real32 > const uz_np1 = nodeManager.getField< elasticfields::Displacementz_np1 >();

  arrayView1d< real32 > const stiffnessVectorx = nodeManager.getField< elasticfields::StiffnessVectorx >();
  arrayView1d< real32 > const stiffnessVectory = nodeManager.getField< elasticfields::StiffnessVectory >();
  arrayView1d< real32 > const stiffnessVectorz = nodeManager.getField< elasticfields::StiffnessVectorz >();

  arrayView1d< real32 > const rhsx = nodeManager.getField< elasticfields::ForcingRHSx >();
  arrayView1d< real32 > const rhsy = nodeManager.getField< elasticfields::ForcingRHSy >();
  arrayView1d< real32 > const rhsz = nodeManager.getField< elasticfields::ForcingRHSz >();

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
    arrayView1d< real32 const > const ux_n   = nodeManager.getField< elasticfields::Displacementx_n >();
    arrayView1d< real32 const > const ux_np1 = nodeManager.getField< elasticfields::Displacementx_np1 >();
    arrayView1d< real32 const > const uy_n   = nodeManager.getField< elasticfields::Displacementy_n >();
    arrayView1d< real32 const > const uy_np1 = nodeManager.getField< elasticfields::Displacementy_np1 >();
    arrayView1d< real32 const > const uz_n   = nodeManager.getField< elasticfields::Displacementz_n >();
    arrayView1d< real32 const > const uz_np1 = nodeManager.getField< elasticfields::Displacementz_np1 >();

    if( m_useDAS == WaveSolverUtils::DASType::none )
    {
      arrayView2d< real32 > const uXReceivers  = m_displacementXNp1AtReceivers.toView();
      arrayView2d< real32 > const uYReceivers  = m_displacementYNp1AtReceivers.toView();
      arrayView2d< real32 > const uZReceivers  = m_displacementZNp1AtReceivers.toView();
      computeAllSeismoTraces( time_n, 0.0, ux_np1, ux_n, uXReceivers );
      computeAllSeismoTraces( time_n, 0.0, uy_np1, uy_n, uYReceivers );
      computeAllSeismoTraces( time_n, 0.0, uz_np1, uz_n, uZReceivers );
      WaveSolverUtils::writeSeismoTraceVector( "seismoTraceReceiver", getName(), m_outputSeismoTrace, m_receiverConstants.size( 0 ),
                                               m_receiverIsLocal, m_nsamplesSeismoTrace, uXReceivers, uYReceivers, uZReceivers );
    }
    else
    {
      arrayView2d< real32 > const dasReceivers  = m_dasSignalNp1AtReceivers.toView();
      computeAllSeismoTraces( time_n, 0.0, ux_np1, ux_n, dasReceivers, m_linearDASVectorX.toView(), true );
      computeAllSeismoTraces( time_n, 0.0, uy_np1, uy_n, dasReceivers, m_linearDASVectorY.toView(), true );
      computeAllSeismoTraces( time_n, 0.0, uz_np1, uz_n, dasReceivers, m_linearDASVectorZ.toView(), true );
      // sum contributions from all MPI ranks, since some receivers might be split among multiple ranks
      MpiWrapper::allReduce( dasReceivers.data(),
                             dasReceivers.data(),
                             m_linearDASGeometry.size( 0 ),
                             MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
                             MPI_COMM_GEOSX );
      WaveSolverUtils::writeSeismoTrace( "dasTraceReceiver", getName(), m_outputSeismoTrace, m_linearDASGeometry.size( 0 ),
                                         m_receiverIsLocal, m_nsamplesSeismoTrace, dasReceivers );
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
