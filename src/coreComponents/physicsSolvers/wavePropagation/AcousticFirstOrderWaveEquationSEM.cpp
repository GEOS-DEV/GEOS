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
 * @file AcousticFirstOrderWaveEquationSEM.cpp
 */

#include "AcousticFirstOrderWaveEquationSEM.hpp"
#include "AcousticFirstOrderWaveEquationSEMKernel.hpp"


#include "finiteElement/FiniteElementDiscretization.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/ElementType.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

namespace geos
{

using namespace dataRepository;
using namespace fields;
//using namespace wavesolverfields;

AcousticFirstOrderWaveEquationSEM::AcousticFirstOrderWaveEquationSEM( const std::string & name,
                                                                      Group * const parent ):
  WaveSolverBase( name,
                  parent )
{

  registerWrapper( viewKeyStruct::pressureNp1AtReceiversString(), &m_pressureNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Pressure value at each receiver for each timestep" );

  registerWrapper( viewKeyStruct::uxNp1AtReceiversString(), &m_uxNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Ux value at each receiver for each timestep" );

  registerWrapper( viewKeyStruct::uyNp1AtReceiversString(), &m_uyNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Uy value at each receiver for each timestep" );

  registerWrapper( viewKeyStruct::uzNp1AtReceiversString(), &m_uzNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Uz value at each receiver for each timestep" );

  registerWrapper( viewKeyStruct::sourceElemString(), &m_sourceElem ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Element containing the sources" );

  registerWrapper( viewKeyStruct::receiverElemString(), &m_rcvElem ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Element containing the receivers" );

  registerWrapper( viewKeyStruct::sourceRegionString(), &m_sourceRegion ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Region containing the sources" );

  registerWrapper( viewKeyStruct::receiverRegionString(), &m_receiverRegion ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Region containing the receivers" );

}

AcousticFirstOrderWaveEquationSEM::~AcousticFirstOrderWaveEquationSEM()
{
  // TODO Auto-generated destructor stub
}

void AcousticFirstOrderWaveEquationSEM::initializePreSubGroups()
{
  WaveSolverBase::initializePreSubGroups();


}


void AcousticFirstOrderWaveEquationSEM::registerDataOnMesh( Group & meshBodies )
{
  WaveSolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    nodeManager.registerField< wavesolverfields::Pressure_np1,
                               wavesolverfields::ForcingRHS,
                               wavesolverfields::MassVector,
                               wavesolverfields::DampingVector,
                               wavesolverfields::FreeSurfaceNodeIndicator >( this->getName() );

    FaceManager & faceManager = mesh.getFaceManager();
    faceManager.registerField< wavesolverfields::FreeSurfaceFaceIndicator >( this->getName() );

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerField< wavesolverfields::MediumVelocity >( this->getName() );
      subRegion.registerField< wavesolverfields::MediumDensity >( this->getName() );

      subRegion.registerField< wavesolverfields::Velocity_x >( this->getName() );
      subRegion.registerField< wavesolverfields::Velocity_y >( this->getName() );
      subRegion.registerField< wavesolverfields::Velocity_z >( this->getName() );

      finiteElement::FiniteElementBase const & fe = subRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {

        using FE_TYPE = TYPEOFREF( finiteElement );

        constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;

        subRegion.getField< wavesolverfields::Velocity_x >().resizeDimension< 1 >( numNodesPerElem );
        subRegion.getField< wavesolverfields::Velocity_y >().resizeDimension< 1 >( numNodesPerElem );
        subRegion.getField< wavesolverfields::Velocity_z >().resizeDimension< 1 >( numNodesPerElem );

      } );

    } );
  } );
}


void AcousticFirstOrderWaveEquationSEM::postProcessInput()
{
  WaveSolverBase::postProcessInput();

  localIndex const numSourcesGlobal = m_sourceCoordinates.size( 0 );
  localIndex const numReceiversGlobal = m_receiverCoordinates.size( 0 );
  m_pressureNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );
  m_uxNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );
  m_uyNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );
  m_uzNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );
  m_sourceElem.resize( numSourcesGlobal );
  m_sourceRegion.resize( numSourcesGlobal );
  m_rcvElem.resize( numReceiversGlobal );
  m_receiverRegion.resize( numReceiversGlobal );
}

void AcousticFirstOrderWaveEquationSEM::precomputeSourceAndReceiverTerm( MeshLevel & mesh, arrayView1d< string const > const & regionNames )
{
  NodeManager const & nodeManager = mesh.getNodeManager();
  FaceManager const & faceManager = mesh.getFaceManager();

  arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const
  X = nodeManager.getField< fields::referencePosition32 >().toViewConst();
  arrayView2d< real64 const > const faceNormal  = faceManager.faceNormal();
  arrayView2d< real64 const > const faceCenter  = faceManager.faceCenter();

  arrayView2d< real64 const > const sourceCoordinates = m_sourceCoordinates.toViewConst();
  arrayView2d< localIndex > const sourceNodeIds = m_sourceNodeIds.toView();
  arrayView2d< real64 > const sourceConstants = m_sourceConstants.toView();
  arrayView1d< localIndex > const sourceIsAccessible = m_sourceIsAccessible.toView();
  arrayView1d< localIndex > const sourceElem = m_sourceElem.toView();
  arrayView1d< localIndex > const sourceRegion = m_sourceRegion.toView();
  sourceNodeIds.setValues< EXEC_POLICY >( -1 );
  sourceConstants.setValues< EXEC_POLICY >( -1 );
  sourceIsAccessible.zero();

  arrayView2d< real64 const > const receiverCoordinates = m_receiverCoordinates.toViewConst();
  arrayView2d< localIndex > const receiverNodeIds = m_receiverNodeIds.toView();
  arrayView2d< real64 > const receiverConstants = m_receiverConstants.toView();
  arrayView1d< localIndex > const receiverIsLocal = m_receiverIsLocal.toView();
  arrayView1d< localIndex > const rcvElem = m_rcvElem.toView();
  arrayView1d< localIndex > const receiverRegion = m_receiverRegion.toView();
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

  mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const regionIndex,
                                                                                        CellElementSubRegion & elementSubRegion )
  {
    GEOS_THROW_IF( elementSubRegion.getElementType() != ElementType::Hexahedron,
                   "Invalid type of element, the acoustic solver is designed for hexahedral meshes only (C3D8) ",
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

      acousticFirstOrderWaveEquationSEMKernels::
        PrecomputeSourceAndReceiverKernel::
        launch< EXEC_POLICY, FE_TYPE >
        ( elementSubRegion.size(),
        regionIndex,
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
        sourceElem,
        sourceNodeIds,
        sourceConstants,
        sourceRegion,
        receiverCoordinates,
        receiverIsLocal,
        rcvElem,
        receiverNodeIds,
        receiverConstants,
        receiverRegion,
        sourceValue,
        dt,
        timeSourceFrequency,
        rickerOrder );
    } );
  } );

}


void AcousticFirstOrderWaveEquationSEM::addSourceToRightHandSide( integer const & cycleNumber, arrayView1d< real32 > const rhs )
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

void AcousticFirstOrderWaveEquationSEM::initializePostInitialConditionsPreSubGroups()
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
    arrayView1d< integer > const & facesDomainBoundaryIndicator = faceManager.getDomainBoundaryIndicator();
    arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.getField< fields::referencePosition32 >().toViewConst();

    /// get table containing face to nodes map
    ArrayOfArraysView< localIndex const > const facesToNodes = faceManager.nodeList().toViewConst();

    // mass matrix to be computed in this function
    arrayView1d< real32 > const mass = nodeManager.getField< wavesolverfields::MassVector >();

    /// damping matrix to be computed for each dof in the boundary of the mesh
    arrayView1d< real32 > const damping = nodeManager.getField< wavesolverfields::DampingVector >();
    damping.zero();
    mass.zero();

    /// get array of indicators: 1 if face is on the free surface; 0 otherwise
    arrayView1d< localIndex const > const freeSurfaceFaceIndicator = faceManager.getField< wavesolverfields::FreeSurfaceFaceIndicator >();

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          CellElementSubRegion & elementSubRegion )
    {
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = elementSubRegion.nodeList();
      arrayView2d< localIndex const > const facesToElements = faceManager.elementList();
      arrayView1d< real32 const > const velocity = elementSubRegion.getField< wavesolverfields::MediumVelocity >();
      arrayView1d< real32 const > const density = elementSubRegion.getField< wavesolverfields::MediumDensity >();

      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        acousticFirstOrderWaveEquationSEMKernels::MassMatrixKernel< FE_TYPE > kernelM( finiteElement );

        kernelM.template launch< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                               X,
                                                               elemsToNodes,
                                                               velocity,
                                                               density,
                                                               mass );

        acousticFirstOrderWaveEquationSEMKernels::DampingMatrixKernel< FE_TYPE > kernelD( finiteElement );

        kernelD.template launch< EXEC_POLICY, ATOMIC_POLICY >( faceManager.size(),
                                                               X,
                                                               facesToElements,
                                                               facesToNodes,
                                                               facesDomainBoundaryIndicator,
                                                               freeSurfaceFaceIndicator,
                                                               velocity,
                                                               damping );
      } );
    } );
  } );
}


void AcousticFirstOrderWaveEquationSEM::applyFreeSurfaceBC( real64 const time, DomainPartition & domain )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  FaceManager & faceManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getFaceManager();
  NodeManager & nodeManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getNodeManager();

  arrayView1d< real32 > const p_np1 = nodeManager.getField< wavesolverfields::Pressure_np1 >();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  /// array of indicators: 1 if a face is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceFaceIndicator = faceManager.getField< wavesolverfields::FreeSurfaceFaceIndicator >();

  /// array of indicators: 1 if a node is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceNodeIndicator = nodeManager.getField< wavesolverfields::FreeSurfaceNodeIndicator >();


  freeSurfaceFaceIndicator.zero();
  freeSurfaceNodeIndicator.zero();

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
        }
      }
    }
    else
    {
      GEOS_ERROR( "This option is not supported yet" );
    }
  } );
}

// Here for retrocompatibily
real64 AcousticFirstOrderWaveEquationSEM::explicitStepForward( real64 const & time_n,
                                                               real64 const & dt,
                                                               integer cycleNumber,
                                                               DomainPartition & domain,
                                                               bool GEOS_UNUSED_PARAM( computeGradient ) )
{
  real64 dtOut = explicitStepInternal( time_n, dt, cycleNumber, domain );
  return dtOut;
}



real64 AcousticFirstOrderWaveEquationSEM::explicitStepBackward( real64 const & time_n,
                                                                real64 const & dt,
                                                                integer cycleNumber,
                                                                DomainPartition & domain,
                                                                bool GEOS_UNUSED_PARAM( computeGradient ) )
{
  GEOS_ERROR( "Backward propagation for the first-order wave propagator not yet implemented" );
  real64 dtOut = explicitStepInternal( time_n, dt, cycleNumber, domain );
  return dtOut;
}


real64 AcousticFirstOrderWaveEquationSEM::explicitStepInternal( real64 const & time_n,
                                                                real64 const & dt,
                                                                integer const cycleNumber,
                                                                DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( time_n, dt, cycleNumber );

  arrayView2d< real64 const > const sourceConstants = m_sourceConstants.toView();
  arrayView1d< localIndex const > const sourceIsAccessible = m_sourceIsAccessible.toView();
  arrayView1d< localIndex const > const sourceElem = m_sourceElem.toView();
  arrayView1d< localIndex const > const sourceRegion = m_sourceRegion.toView();
  arrayView2d< real32 const > const sourceValue = m_sourceValue.toView();

  GEOS_LOG_RANK_0_IF( dt < epsilonLoc, "Warning! Value for dt: " << dt << "s is smaller than local threshold: " << epsilonLoc );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                  [&] ( string const &,
                                        MeshLevel & mesh,
                                        arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.getField< fields::referencePosition32 >().toViewConst();

    arrayView1d< real32 const > const mass = nodeManager.getField< wavesolverfields::MassVector >();
    arrayView1d< real32 const > const damping = nodeManager.getField< wavesolverfields::DampingVector >();

    arrayView1d< real32 > const p_np1 = nodeManager.getField< wavesolverfields::Pressure_np1 >();

    arrayView1d< real32 > const rhs = nodeManager.getField< wavesolverfields::ForcingRHS >();

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const regionIndex,
                                                                                          CellElementSubRegion & elementSubRegion )
    {
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();
      arrayView1d< real32 const > const density = elementSubRegion.getField< wavesolverfields::MediumDensity >();
      arrayView2d< real32 > const velocity_x = elementSubRegion.getField< wavesolverfields::Velocity_x >();
      arrayView2d< real32 > const velocity_y = elementSubRegion.getField< wavesolverfields::Velocity_y >();
      arrayView2d< real32 > const velocity_z = elementSubRegion.getField< wavesolverfields::Velocity_z >();
      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );



        acousticFirstOrderWaveEquationSEMKernels::
          VelocityComputation< FE_TYPE > kernel( finiteElement );
        kernel.template launch< EXEC_POLICY, ATOMIC_POLICY >
          ( elementSubRegion.size(),
          X,
          elemsToNodes,
          p_np1,
          density,
          dt,
          velocity_x,
          velocity_y,
          velocity_z );
        acousticFirstOrderWaveEquationSEMKernels::
          PressureComputation< FE_TYPE > kernel2( finiteElement );
        kernel2.template launch< EXEC_POLICY, ATOMIC_POLICY >
          ( elementSubRegion.size(),
          regionIndex,
          nodeManager.size(),
          X,
          elemsToNodes,
          velocity_x,
          velocity_y,
          velocity_z,
          mass,
          damping,
          sourceConstants,
          sourceValue,
          sourceIsAccessible,
          sourceElem,
          sourceRegion,
          dt,
          cycleNumber,
          p_np1 );
      } );
      arrayView2d< real32 > const uxReceivers   = m_uxNp1AtReceivers.toView();
      arrayView2d< real32 > const uyReceivers   = m_uyNp1AtReceivers.toView();
      arrayView2d< real32 > const uzReceivers   = m_uzNp1AtReceivers.toView();

      compute2dVariableAllSeismoTraces( regionIndex, time_n, dt, velocity_x, velocity_x, uxReceivers );
      compute2dVariableAllSeismoTraces( regionIndex, time_n, dt, velocity_y, velocity_y, uyReceivers );
      compute2dVariableAllSeismoTraces( regionIndex, time_n, dt, velocity_z, velocity_z, uzReceivers );

    } );

    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addFields( FieldLocation::Node, { wavesolverfields::Pressure_np1::key() } );
    fieldsToBeSync.addElementFields( {wavesolverfields::Velocity_x::key(), wavesolverfields::Velocity_y::key(), wavesolverfields::Velocity_z::key()}, regionNames );

    CommunicationTools & syncFields = CommunicationTools::getInstance();
    syncFields.synchronizeFields( fieldsToBeSync,
                                  mesh,
                                  domain.getNeighbors(),
                                  true );

    // compute the seismic traces since last step.
    arrayView2d< real32 > const pReceivers   = m_pressureNp1AtReceivers.toView();
    computeAllSeismoTraces( time_n, dt, p_np1, p_np1, pReceivers );

    // increment m_indexSeismoTrace
    while( (m_dtSeismoTrace*m_indexSeismoTrace) <= (time_n + epsilonLoc) && m_indexSeismoTrace < m_nsamplesSeismoTrace )
    {
      m_indexSeismoTrace++;
    }

  } );
  return dt;
}

void AcousticFirstOrderWaveEquationSEM::cleanup( real64 const time_n, integer const, integer const, real64 const, DomainPartition & domain )
{
  // compute the remaining seismic traces, if needed
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();
    arrayView1d< real32 const > const p_np1 = nodeManager.getField< wavesolverfields::Pressure_np1 >();
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const regionIndex,
                                                                                          CellElementSubRegion & elementSubRegion )
    {
      arrayView2d< real32 > const velocity_x = elementSubRegion.getField< wavesolverfields::Velocity_x >();
      arrayView2d< real32 > const velocity_y = elementSubRegion.getField< wavesolverfields::Velocity_y >();
      arrayView2d< real32 > const velocity_z = elementSubRegion.getField< wavesolverfields::Velocity_z >();

      arrayView2d< real32 > const uxReceivers   = m_uxNp1AtReceivers.toView();
      arrayView2d< real32 > const uyReceivers   = m_uyNp1AtReceivers.toView();
      arrayView2d< real32 > const uzReceivers   = m_uzNp1AtReceivers.toView();

      compute2dVariableAllSeismoTraces( regionIndex, time_n, 0, velocity_x, velocity_x, uxReceivers );
      compute2dVariableAllSeismoTraces( regionIndex, time_n, 0, velocity_y, velocity_y, uyReceivers );
      compute2dVariableAllSeismoTraces( regionIndex, time_n, 0, velocity_z, velocity_z, uzReceivers );

    } );
    arrayView2d< real32 > const pReceivers   = m_pressureNp1AtReceivers.toView();
    computeAllSeismoTraces( time_n, 0, p_np1, p_np1, pReceivers );
  } );

  // increment m_indexSeismoTrace
  while( (m_dtSeismoTrace*m_indexSeismoTrace) <= (time_n + epsilonLoc) && m_indexSeismoTrace < m_nsamplesSeismoTrace )
  {
    m_indexSeismoTrace++;
  }
}

void AcousticFirstOrderWaveEquationSEM::computeAllSeismoTraces( real64 const time_n,
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
    WaveSolverUtils::computeSeismoTrace( time_n, dt, timeSeismo, indexSeismoTrace, m_receiverNodeIds, m_receiverConstants, m_receiverIsLocal, m_nsamplesSeismoTrace,
                                         m_outputSeismoTrace,
                                         var_np1, var_n, varAtReceivers );
  }
}

void AcousticFirstOrderWaveEquationSEM::compute2dVariableAllSeismoTraces( localIndex const regionIndex,
                                                                          real64 const time_n,
                                                                          real64 const dt,
                                                                          arrayView2d< real32 const > const var_np1,
                                                                          arrayView2d< real32 const > const var_n,
                                                                          arrayView2d< real32 > varAtReceivers )
{
  localIndex indexSeismoTrace = m_indexSeismoTrace;
  for( real64 timeSeismo;
       (timeSeismo = m_dtSeismoTrace*indexSeismoTrace) <= (time_n + epsilonLoc) && indexSeismoTrace < m_nsamplesSeismoTrace;
       indexSeismoTrace++ )
  {
    WaveSolverUtils::compute2dVariableSeismoTrace( time_n, dt, regionIndex, m_receiverRegion, timeSeismo, indexSeismoTrace, m_rcvElem, m_receiverConstants, m_receiverIsLocal, m_nsamplesSeismoTrace,
                                                   m_outputSeismoTrace,
                                                   var_np1, var_n, varAtReceivers );
  }
}


void AcousticFirstOrderWaveEquationSEM::initializePML()
{
  GEOS_ERROR( "PML for the first order acoustic wave propagator not yet implemented" );
}

void AcousticFirstOrderWaveEquationSEM::applyPML( real64 const, DomainPartition & )
{
  GEOS_ERROR( "PML for the first order acoustic wave propagator not yet implemented" );
}

REGISTER_CATALOG_ENTRY( SolverBase, AcousticFirstOrderWaveEquationSEM, string const &, dataRepository::Group * const )

} /* namespace geos */
