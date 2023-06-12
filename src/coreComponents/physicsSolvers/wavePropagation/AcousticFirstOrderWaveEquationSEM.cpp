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

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    nodeManager.registerField< acousticFirstOrderSemFields::Pressure_np1,
                               matricialFields::MassVector,
                               matricialFields::DampingVector,
                               geophysicalFields::FreeSurfaceNodeIndicator >( this->getName() );

    FaceManager & faceManager = mesh.getFaceManager();
    faceManager.registerField< geophysicalFields::FreeSurfaceFaceIndicator >( this->getName() );

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerField< geophysicalFields::Pwavespeed >( this->getName() );
      subRegion.registerField< geophysicalFields::Density >( this->getName() );

      subRegion.registerField< acousticFirstOrderSemFields::Velocity_x >( this->getName() );
      subRegion.registerField< acousticFirstOrderSemFields::Velocity_y >( this->getName() );
      subRegion.registerField< acousticFirstOrderSemFields::Velocity_z >( this->getName() );

      finiteElement::FiniteElementBase const & fe = subRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {

        using FE_TYPE = TYPEOFREF( finiteElement );

        constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;

        subRegion.getField< acousticFirstOrderSemFields::Velocity_x >().resizeDimension< 1 >( numNodesPerElem );
        subRegion.getField< acousticFirstOrderSemFields::Velocity_y >().resizeDimension< 1 >( numNodesPerElem );
        subRegion.getField< acousticFirstOrderSemFields::Velocity_z >().resizeDimension< 1 >( numNodesPerElem );

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
  m_rcvElem.resize( numReceiversGlobal );


}

void AcousticFirstOrderWaveEquationSEM::precomputeSourceAndReceiverTerm( MeshLevel & mesh, arrayView1d< string const > const & regionNames )
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
  arrayView1d< localIndex > const sourceElem = m_sourceElem.toView();
  sourceNodeIds.setValues< EXEC_POLICY >( -1 );
  sourceConstants.setValues< EXEC_POLICY >( -1 );
  sourceIsAccessible.zero();

  arrayView2d< real64 const > const receiverCoordinates = m_receiverCoordinates.toViewConst();
  arrayView2d< localIndex > const receiverNodeIds = m_receiverNodeIds.toView();
  arrayView2d< real64 > const receiverConstants = m_receiverConstants.toView();
  arrayView1d< localIndex > const receiverIsLocal = m_receiverIsLocal.toView();
  arrayView1d< localIndex > const rcvElem = m_rcvElem.toView();
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
        receiverCoordinates,
        receiverIsLocal,
        rcvElem,
        receiverNodeIds,
        receiverConstants,
        sourceValue,
        dt,
        timeSourceFrequency,
        rickerOrder );
    } );
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
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

    /// get table containing face to nodes map
    ArrayOfArraysView< localIndex const > const facesToNodes = faceManager.nodeList().toViewConst();

    // mass matrix to be computed in this function
    arrayView1d< real32 > const mass = nodeManager.getField< matricialFields::MassVector >();

    /// damping matrix to be computed for each dof in the boundary of the mesh
    arrayView1d< real32 > const damping = nodeManager.getField< matricialFields::DampingVector >();
    damping.zero();
    mass.zero();

    /// get array of indicators: 1 if face is on the free surface; 0 otherwise
    arrayView1d< localIndex const > const freeSurfaceFaceIndicator = faceManager.getField< geophysicalFields::FreeSurfaceFaceIndicator >();

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          CellElementSubRegion & elementSubRegion )
    {
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = elementSubRegion.nodeList();
      arrayView2d< localIndex const > const facesToElements = faceManager.elementList();
      arrayView1d< real32 const > const Vp = elementSubRegion.getField< geophysicalFields::Pwavespeed >();
      arrayView1d< real32 const > const density = elementSubRegion.getField< geophysicalFields::Density >();

      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        acousticFirstOrderWaveEquationSEMKernels::MassMatrixKernel< FE_TYPE > kernelM( finiteElement );

        kernelM.template launch< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                               X,
                                                               elemsToNodes,
                                                               Vp,
                                                               density,
                                                               mass );

        acousticFirstOrderWaveEquationSEMKernels::DampingMatrixKernel< FE_TYPE > kernelD( finiteElement );

        kernelD.template launch< EXEC_POLICY, ATOMIC_POLICY >( faceManager.size(),
                                                               X,
                                                               facesToElements,
                                                               facesToNodes,
                                                               facesDomainBoundaryIndicator,
                                                               freeSurfaceFaceIndicator,
                                                               Vp,
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

  arrayView1d< real32 > const p_np1 = nodeManager.getField< acousticFirstOrderSemFields::Pressure_np1 >();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  /// array of indicators: 1 if a face is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceFaceIndicator = faceManager.getField< geophysicalFields::FreeSurfaceFaceIndicator >();

  /// array of indicators: 1 if a node is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceNodeIndicator = nodeManager.getField< geophysicalFields::FreeSurfaceNodeIndicator >();

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
  arrayView2d< real32 const > const sourceValue = m_sourceValue.toView();

  GEOS_LOG_RANK_0_IF( dt < epsilonLoc, "Warning! Value for dt: " << dt << "s is smaller than local threshold: " << epsilonLoc );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                  [&] ( string const &,
                                        MeshLevel & mesh,
                                        arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

    arrayView1d< real32 const > const mass = nodeManager.getField< matricialFields::MassVector >();
    arrayView1d< real32 const > const damping = nodeManager.getField< matricialFields::DampingVector >();

    arrayView1d< real32 > const p_np1 = nodeManager.getField< acousticFirstOrderSemFields::Pressure_np1 >();

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          CellElementSubRegion & elementSubRegion )
    {
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();
      arrayView1d< real32 const > const density = elementSubRegion.getField< geophysicalFields::Density >();
      arrayView2d< real32 > const velocity_x = elementSubRegion.getField< acousticFirstOrderSemFields::Velocity_x >();
      arrayView2d< real32 > const velocity_y = elementSubRegion.getField< acousticFirstOrderSemFields::Velocity_y >();
      arrayView2d< real32 > const velocity_z = elementSubRegion.getField< acousticFirstOrderSemFields::Velocity_z >();
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
          dt,
          cycleNumber,
          p_np1 );
      } );
      arrayView2d< real32 > const uxReceivers   = m_uxNp1AtReceivers.toView();
      arrayView2d< real32 > const uyReceivers   = m_uyNp1AtReceivers.toView();
      arrayView2d< real32 > const uzReceivers   = m_uzNp1AtReceivers.toView();

      compute2dVariableAllSeismoTraces( time_n, dt, velocity_x, velocity_x, uxReceivers );
      compute2dVariableAllSeismoTraces( time_n, dt, velocity_y, velocity_y, uyReceivers );
      compute2dVariableAllSeismoTraces( time_n, dt, velocity_z, velocity_z, uzReceivers );

    } );

    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addFields( FieldLocation::Node, { acousticFirstOrderSemFields::Pressure_np1::key() } );
    fieldsToBeSync.addElementFields( {acousticFirstOrderSemFields::Velocity_x::key(), acousticFirstOrderSemFields::Velocity_y::key(), acousticFirstOrderSemFields::Velocity_z::key()}, regionNames );

    CommunicationTools & syncFields = CommunicationTools::getInstance();
    syncFields.synchronizeFields( fieldsToBeSync,
                                  mesh,
                                  domain.getNeighbors(),
                                  true );

    // compute the seismic traces since last step.
    arrayView2d< real32 > const pReceivers   = m_pressureNp1AtReceivers.toView();
    computeAllSeismoTraces( time_n, dt, p_np1, p_np1, pReceivers);

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
    arrayView1d< real32 const > const p_np1 = nodeManager.getField< acousticFirstOrderSemFields::Pressure_np1 >();
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          CellElementSubRegion & elementSubRegion )
    {
      arrayView2d< real32 > const velocity_x = elementSubRegion.getField< acousticFirstOrderSemFields::Velocity_x >();
      arrayView2d< real32 > const velocity_y = elementSubRegion.getField< acousticFirstOrderSemFields::Velocity_y >();
      arrayView2d< real32 > const velocity_z = elementSubRegion.getField< acousticFirstOrderSemFields::Velocity_z >();

      arrayView2d< real32 > const uxReceivers   = m_uxNp1AtReceivers.toView();
      arrayView2d< real32 > const uyReceivers   = m_uyNp1AtReceivers.toView();
      arrayView2d< real32 > const uzReceivers   = m_uzNp1AtReceivers.toView();

      compute2dVariableAllSeismoTraces( time_n, 0, velocity_x, velocity_x, uxReceivers );
      compute2dVariableAllSeismoTraces( time_n, 0, velocity_y, velocity_y, uyReceivers );
      compute2dVariableAllSeismoTraces( time_n, 0, velocity_z, velocity_z, uzReceivers );

    } );
    arrayView2d< real32 > const pReceivers   = m_pressureNp1AtReceivers.toView();
    computeAllSeismoTraces( time_n, 0, p_np1, p_np1, pReceivers);

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

void AcousticFirstOrderWaveEquationSEM::compute2dVariableAllSeismoTraces( real64 const time_n,
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
    WaveSolverUtils::compute2dVariableSeismoTrace( time_n, dt, timeSeismo, indexSeismoTrace, m_rcvElem, m_receiverConstants, m_receiverIsLocal, m_nsamplesSeismoTrace,
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
