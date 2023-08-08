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
 * @file ElasticFirstOrderWaveEquationSEM.cpp
 */

#include "ElasticFirstOrderWaveEquationSEM.hpp"
#include "ElasticFirstOrderWaveEquationSEMKernel.hpp"

#include "finiteElement/FiniteElementDiscretization.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/ElementType.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

namespace geos
{

using namespace dataRepository;
using namespace fields;

ElasticFirstOrderWaveEquationSEM::ElasticFirstOrderWaveEquationSEM( const std::string & name,
                                                                    Group * const parent ):

  WaveSolverBase( name,
                  parent )
{

  registerWrapper( viewKeyStruct::displacementxNp1AtReceiversString(), &m_displacementxNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Displacement value at each receiver for each timestep (x-components)" );

  registerWrapper( viewKeyStruct::displacementyNp1AtReceiversString(), &m_displacementyNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Displacement value at each receiver for each timestep (y-components)" );

  registerWrapper( viewKeyStruct::displacementzNp1AtReceiversString(), &m_displacementzNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Displacement value at each receiver for each timestep (z-components)" );

  registerWrapper( viewKeyStruct::displacementzNp1AtReceiversString(), &m_sigmaxxNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Displacement value at each receiver for each timestep (z-components)" );

  registerWrapper( viewKeyStruct::displacementzNp1AtReceiversString(), &m_sigmayyNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Displacement value at each receiver for each timestep (z-components)" );

  registerWrapper( viewKeyStruct::displacementzNp1AtReceiversString(), &m_sigmazzNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Displacement value at each receiver for each timestep (z-components)" );

  registerWrapper( viewKeyStruct::displacementzNp1AtReceiversString(), &m_sigmaxyNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Displacement value at each receiver for each timestep (z-components)" );

  registerWrapper( viewKeyStruct::displacementzNp1AtReceiversString(), &m_sigmaxzNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Displacement value at each receiver for each timestep (z-components)" );

  registerWrapper( viewKeyStruct::displacementzNp1AtReceiversString(), &m_sigmayzNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Displacement value at each receiver for each timestep (z-components)" );

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

ElasticFirstOrderWaveEquationSEM::~ElasticFirstOrderWaveEquationSEM()
{
  // TODO Auto-generated destructor stub
}

void ElasticFirstOrderWaveEquationSEM::initializePreSubGroups()
{
  WaveSolverBase::initializePreSubGroups();
}


void ElasticFirstOrderWaveEquationSEM::registerDataOnMesh( Group & meshBodies )
{
  WaveSolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & )

  {
    NodeManager & nodeManager = mesh.getNodeManager();

    nodeManager.registerField< wavesolverfields::Displacementx_np1,
                               wavesolverfields::Displacementy_np1,
                               wavesolverfields::Displacementz_np1,
                               wavesolverfields::ForcingRHS,
                               wavesolverfields::MassVector,
                               wavesolverfields::DampingVectorx,
                               wavesolverfields::DampingVectory,
                               wavesolverfields::DampingVectorz,
                               wavesolverfields::FreeSurfaceNodeIndicator >( this->getName() );

    FaceManager & faceManager = mesh.getFaceManager();
    faceManager.registerField< wavesolverfields::FreeSurfaceFaceIndicator >( this->getName() );

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerField< wavesolverfields::MediumVelocityVp >( this->getName() );
      subRegion.registerField< wavesolverfields::MediumVelocityVs >( this->getName() );
      subRegion.registerField< wavesolverfields::MediumDensity >( this->getName() );
      subRegion.registerField< wavesolverfields::Lambda >( this->getName() );
      subRegion.registerField< wavesolverfields::Mu >( this->getName() );

      subRegion.registerField< wavesolverfields::Stresstensorxx >( this->getName());
      subRegion.registerField< wavesolverfields::Stresstensoryy >( this->getName());
      subRegion.registerField< wavesolverfields::Stresstensorzz >( this->getName());
      subRegion.registerField< wavesolverfields::Stresstensorxy >( this->getName());
      subRegion.registerField< wavesolverfields::Stresstensorxz >( this->getName());
      subRegion.registerField< wavesolverfields::Stresstensoryz >( this->getName());

      finiteElement::FiniteElementBase const & fe = subRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {

        using FE_TYPE = TYPEOFREF( finiteElement );

        constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;

        subRegion.getField< wavesolverfields::Stresstensorxx >().resizeDimension< 1 >( numNodesPerElem );
        subRegion.getField< wavesolverfields::Stresstensoryy >().resizeDimension< 1 >( numNodesPerElem );
        subRegion.getField< wavesolverfields::Stresstensorzz >().resizeDimension< 1 >( numNodesPerElem );
        subRegion.getField< wavesolverfields::Stresstensorxy >().resizeDimension< 1 >( numNodesPerElem );
        subRegion.getField< wavesolverfields::Stresstensorxz >().resizeDimension< 1 >( numNodesPerElem );
        subRegion.getField< wavesolverfields::Stresstensoryz >().resizeDimension< 1 >( numNodesPerElem );
      } );


    } );

  } );
}



void ElasticFirstOrderWaveEquationSEM::postProcessInput()
{

  WaveSolverBase::postProcessInput();

  localIndex const numSourcesGlobal = m_sourceCoordinates.size( 0 );
  m_sourceElem.resize( numSourcesGlobal );
  m_sourceRegion.resize( numSourcesGlobal );

  localIndex const numReceiversGlobal = m_receiverCoordinates.size( 0 );
  m_rcvElem.resize( numReceiversGlobal );
  m_receiverRegion.resize( numReceiversGlobal );

  m_displacementxNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );
  m_displacementyNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );
  m_displacementzNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );

  m_sigmaxxNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );
  m_sigmayyNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );
  m_sigmazzNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );
  m_sigmaxyNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );
  m_sigmaxzNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );
  m_sigmayzNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );


}


void ElasticFirstOrderWaveEquationSEM::precomputeSourceAndReceiverTerm( MeshLevel & mesh, arrayView1d< string const > const & regionNames )
{
  NodeManager const & nodeManager = mesh.getNodeManager();
  FaceManager const & faceManager = mesh.getFaceManager();

  arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.getField< fields::referencePosition32 >().toViewConst();
  arrayView2d< real64 const > const faceNormal  = faceManager.faceNormal();
  arrayView2d< real64 const > const faceCenter  = faceManager.faceCenter();

  arrayView2d< real64 const > const sourceCoordinates = m_sourceCoordinates.toViewConst();
  arrayView2d< localIndex > const sourceNodeIds = m_sourceNodeIds.toView();
  arrayView2d< real64 > const sourceConstants = m_sourceConstants.toView();
  arrayView1d< localIndex > const sourceIsAccessible = m_sourceIsAccessible.toView();
  arrayView1d< localIndex > const sourceElem = m_sourceElem.toView();
  arrayView1d< localIndex > const sourceRegion = m_sourceRegion.toView();
  sourceNodeIds.setValues< serialPolicy >( -1 );
  sourceConstants.setValues< serialPolicy >( -1 );
  sourceIsAccessible.zero();

  arrayView2d< real64 const > const receiverCoordinates = m_receiverCoordinates.toViewConst();
  arrayView2d< localIndex > const receiverNodeIds = m_receiverNodeIds.toView();
  arrayView2d< real64 > const receiverConstants = m_receiverConstants.toView();
  arrayView1d< localIndex > const receiverIsLocal = m_receiverIsLocal.toView();
  arrayView1d< localIndex > const rcvElem = m_rcvElem.toView();
  arrayView1d< localIndex > const receiverRegion = m_receiverRegion.toView();

  receiverNodeIds.setValues< serialPolicy >( -1 );
  receiverConstants.setValues< serialPolicy >( -1 );
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

      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      localIndex const numFacesPerElem = elementSubRegion.numFacesPerElement();

      elasticFirstOrderWaveEquationSEMKernels::
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


void ElasticFirstOrderWaveEquationSEM::addSourceToRightHandSide( integer const & cycleNumber, arrayView1d< real32 > const rhs )
{
  arrayView2d< localIndex const > const sourceNodeIds = m_sourceNodeIds.toViewConst();
  arrayView2d< real64 const > const sourceConstants   = m_sourceConstants.toViewConst();
  arrayView1d< localIndex const > const sourceIsAccessible = m_sourceIsAccessible.toViewConst();
  arrayView2d< real32 const > const sourceValue   = m_sourceValue.toViewConst();

  GEOS_THROW_IF( cycleNumber > sourceValue.size( 0 ), "Too many steps compared to array size", std::runtime_error );

  forAll< serialPolicy >( m_sourceConstants.size( 0 ), [=] ( localIndex const isrc )
  {
    if( sourceIsAccessible[isrc] == 1 )
    {
      for( localIndex inode = 0; inode < m_sourceConstants.size( 1 ); ++inode )
      {
        real32 const localIncrement = sourceConstants[isrc][inode] * sourceValue[cycleNumber][isrc];
        RAJA::atomicAdd< ATOMIC_POLICY >( &rhs[sourceNodeIds[isrc][inode]], localIncrement );
      }
    }
  } );
}

void ElasticFirstOrderWaveEquationSEM::initializePostInitialConditionsPreSubGroups()
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
    arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.getField< fields::referencePosition32 >().toViewConst();
    arrayView2d< real64 const > const faceNormal  = faceManager.faceNormal();

    /// get face to node map
    ArrayOfArraysView< localIndex const > const facesToNodes = faceManager.nodeList().toViewConst();

    // mass matrix to be computed in this function
    arrayView1d< real32 > const mass = nodeManager.getField< wavesolverfields::MassVector >();
    mass.zero();
    /// damping matrix to be computed for each dof in the boundary of the mesh
    arrayView1d< real32 > const dampingx = nodeManager.getField< wavesolverfields::DampingVectorx >();
    arrayView1d< real32 > const dampingy = nodeManager.getField< wavesolverfields::DampingVectory >();
    arrayView1d< real32 > const dampingz = nodeManager.getField< wavesolverfields::DampingVectorz >();
    dampingx.zero();
    dampingy.zero();
    dampingz.zero();

    /// get array of indicators: 1 if face is on the free surface; 0 otherwise
    arrayView1d< localIndex const > const freeSurfaceFaceIndicator = faceManager.getField< wavesolverfields::FreeSurfaceFaceIndicator >();

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          CellElementSubRegion & elementSubRegion )
    {

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = elementSubRegion.nodeList();
      arrayView2d< localIndex const > const facesToElements = faceManager.elementList();
      arrayView1d< real32 > const density = elementSubRegion.getField< wavesolverfields::MediumDensity >();
      arrayView1d< real32 > const velocityVp = elementSubRegion.getField< wavesolverfields::MediumVelocityVp >();
      arrayView1d< real32 > const velocityVs = elementSubRegion.getField< wavesolverfields::MediumVelocityVs >();

      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe,
                                                                               [&]
                                                                                 ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        elasticFirstOrderWaveEquationSEMKernels::MassMatrixKernel< FE_TYPE > kernelM( finiteElement );

        kernelM.template launch< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                               X,
                                                               elemsToNodes,
                                                               density,
                                                               mass );

        elasticFirstOrderWaveEquationSEMKernels::DampingMatrixKernel< FE_TYPE > kernelD( finiteElement );

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


void ElasticFirstOrderWaveEquationSEM::applyFreeSurfaceBC( real64 const time, DomainPartition & domain )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  FaceManager & faceManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getFaceManager();
  NodeManager & nodeManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getNodeManager();

  arrayView1d< real32 > const ux_np1 = nodeManager.getField< wavesolverfields::Displacementx_np1 >();
  arrayView1d< real32 > const uy_np1 = nodeManager.getField< wavesolverfields::Displacementy_np1 >();
  arrayView1d< real32 > const uz_np1 = nodeManager.getField< wavesolverfields::Displacementz_np1 >();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  /// set array of indicators: 1 if a face is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceFaceIndicator = faceManager.getField< wavesolverfields::FreeSurfaceFaceIndicator >();

  /// set array of indicators: 1 if a node is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceNodeIndicator = nodeManager.getField< wavesolverfields::FreeSurfaceNodeIndicator >();

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
        }
      }
    }
    else
    {
      GEOS_ERROR( "This option is not supported yet" );
    }
  } );
}

real64 ElasticFirstOrderWaveEquationSEM::explicitStepForward( real64 const & time_n,
                                                              real64 const & dt,
                                                              integer cycleNumber,
                                                              DomainPartition & domain,
                                                              bool GEOS_UNUSED_PARAM( computeGradient ) )
{
  real64 dtOut = explicitStepInternal( time_n, dt, cycleNumber, domain );
  return dtOut;
}



real64 ElasticFirstOrderWaveEquationSEM::explicitStepBackward( real64 const & time_n,
                                                               real64 const & dt,
                                                               integer cycleNumber,
                                                               DomainPartition & domain,
                                                               bool GEOS_UNUSED_PARAM( computeGradient ) )
{
  GEOS_ERROR( "Backward propagation for the first order elastic wave propagator not yet implemented" );
  real64 dtOut = explicitStepInternal( time_n, dt, cycleNumber, domain );
  return dtOut;
}

real64 ElasticFirstOrderWaveEquationSEM::explicitStepInternal( real64 const & time_n,
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


  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )

  {

    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.getField< fields::referencePosition32 >().toViewConst();

    arrayView1d< real32 const > const mass = nodeManager.getField< wavesolverfields::MassVector >();
    arrayView1d< real32 > const dampingx = nodeManager.getField< wavesolverfields::DampingVectorx >();
    arrayView1d< real32 > const dampingy = nodeManager.getField< wavesolverfields::DampingVectory >();
    arrayView1d< real32 > const dampingz = nodeManager.getField< wavesolverfields::DampingVectorz >();


    arrayView1d< real32 > const ux_np1 = nodeManager.getField< wavesolverfields::Displacementx_np1 >();
    arrayView1d< real32 > const uy_np1 = nodeManager.getField< wavesolverfields::Displacementy_np1 >();
    arrayView1d< real32 > const uz_np1 = nodeManager.getField< wavesolverfields::Displacementz_np1 >();

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const regionIndex,
                                                                                          CellElementSubRegion & elementSubRegion )
    {

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();

      arrayView1d< real32 const > const velocityVp = elementSubRegion.getField< wavesolverfields::MediumVelocityVp >();
      arrayView1d< real32 const > const velocityVs = elementSubRegion.getField< wavesolverfields::MediumVelocityVs >();
      arrayView1d< real32 const > const density = elementSubRegion.getField< wavesolverfields::MediumDensity >();

      arrayView1d< real32 > const lambda = elementSubRegion.getField< wavesolverfields::Lambda >();
      arrayView1d< real32 > const mu = elementSubRegion.getField< wavesolverfields::Mu >();

      arrayView2d< real32 > const stressxx = elementSubRegion.getField< wavesolverfields::Stresstensorxx >();
      arrayView2d< real32 > const stressyy = elementSubRegion.getField< wavesolverfields::Stresstensoryy >();
      arrayView2d< real32 > const stresszz = elementSubRegion.getField< wavesolverfields::Stresstensorzz >();
      arrayView2d< real32 > const stressxy = elementSubRegion.getField< wavesolverfields::Stresstensorxy >();
      arrayView2d< real32 > const stressxz = elementSubRegion.getField< wavesolverfields::Stresstensorxz >();
      arrayView2d< real32 > const stressyz = elementSubRegion.getField< wavesolverfields::Stresstensoryz >();


      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        elasticFirstOrderWaveEquationSEMKernels::
          StressComputation< FE_TYPE > kernel( finiteElement );
        kernel.template launch< EXEC_POLICY, ATOMIC_POLICY >
          ( elementSubRegion.size(),
          regionIndex,
          X,
          elemsToNodes,
          ux_np1,
          uy_np1,
          uz_np1,
          density,
          velocityVp,
          velocityVs,
          lambda,
          mu,
          sourceConstants,
          sourceIsAccessible,
          sourceElem,
          sourceRegion,
          sourceValue,
          dt,
          cycleNumber,
          stressxx,
          stressyy,
          stresszz,
          stressxy,
          stressxz,
          stressyz );

        elasticFirstOrderWaveEquationSEMKernels::
          VelocityComputation< FE_TYPE > kernel2( finiteElement );
        kernel2.template launch< EXEC_POLICY, ATOMIC_POLICY >
          ( elementSubRegion.size(),
          nodeManager.size(),
          X,
          elemsToNodes,
          stressxx,
          stressyy,
          stresszz,
          stressxy,
          stressxz,
          stressyz,
          mass,
          dampingx,
          dampingy,
          dampingz,
          dt,
          ux_np1,
          uy_np1,
          uz_np1 );


      } );

      arrayView2d< real32 > const sigmaxxReceivers   = m_sigmaxxNp1AtReceivers.toView();
      arrayView2d< real32 > const sigmayyReceivers   = m_sigmayyNp1AtReceivers.toView();
      arrayView2d< real32 > const sigmazzReceivers   = m_sigmazzNp1AtReceivers.toView();
      arrayView2d< real32 > const sigmaxyReceivers   = m_sigmaxyNp1AtReceivers.toView();
      arrayView2d< real32 > const sigmaxzReceivers   = m_sigmaxzNp1AtReceivers.toView();
      arrayView2d< real32 > const sigmayzReceivers   = m_sigmayzNp1AtReceivers.toView();

      compute2dVariableAllSeismoTraces( regionIndex, time_n, dt, stressxx, stressxx, sigmaxxReceivers );
      compute2dVariableAllSeismoTraces( regionIndex, time_n, dt, stressyy, stressyy, sigmayyReceivers );
      compute2dVariableAllSeismoTraces( regionIndex, time_n, dt, stresszz, stresszz, sigmazzReceivers );
      compute2dVariableAllSeismoTraces( regionIndex, time_n, dt, stressxy, stressxy, sigmaxyReceivers );
      compute2dVariableAllSeismoTraces( regionIndex, time_n, dt, stressxz, stressxz, sigmaxzReceivers );
      compute2dVariableAllSeismoTraces( regionIndex, time_n, dt, stressyz, stressyz, sigmayzReceivers );


    } );

    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addFields( FieldLocation::Node, { wavesolverfields::Displacementx_np1::key(), wavesolverfields::Displacementy_np1::key(), wavesolverfields::Displacementz_np1::key()} );
    fieldsToBeSync.addElementFields( {wavesolverfields::Stresstensorxx::key(), wavesolverfields::Stresstensoryy::key(), wavesolverfields::Stresstensorzz::key(),
                                      wavesolverfields::Stresstensorxy::key(),
                                      wavesolverfields::Stresstensorxz::key(), wavesolverfields::Stresstensoryz::key()}, regionNames );


    CommunicationTools & syncFields = CommunicationTools::getInstance();
    syncFields.synchronizeFields( fieldsToBeSync,
                                  domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ),
                                  domain.getNeighbors(),
                                  true );

    // compute the seismic traces since last step.
    arrayView2d< real32 > const uxReceivers   = m_displacementxNp1AtReceivers.toView();
    arrayView2d< real32 > const uyReceivers   = m_displacementyNp1AtReceivers.toView();
    arrayView2d< real32 > const uzReceivers   = m_displacementzNp1AtReceivers.toView();

    computeAllSeismoTraces( time_n, dt, ux_np1, ux_np1, uxReceivers );
    computeAllSeismoTraces( time_n, dt, uy_np1, uy_np1, uyReceivers );
    computeAllSeismoTraces( time_n, dt, uz_np1, uz_np1, uzReceivers );

    // increment m_indexSeismoTrace
    while( (m_dtSeismoTrace*m_indexSeismoTrace) <= (time_n + epsilonLoc) && m_indexSeismoTrace < m_nsamplesSeismoTrace )
    {
      m_indexSeismoTrace++;
    }

  } );

  return dt;

}

void ElasticFirstOrderWaveEquationSEM::cleanup( real64 const time_n,
                                                integer const cycleNumber,
                                                integer const eventCounter,
                                                real64 const eventProgress,
                                                DomainPartition & domain )
{
  SolverBase::cleanup( time_n, cycleNumber, eventCounter, eventProgress, domain );

  // compute the remaining seismic traces, if needed
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const regionIndex,
                                                                                          CellElementSubRegion & elementSubRegion )
    {
      arrayView2d< real32 const > const stressxx = elementSubRegion.getField< wavesolverfields::Stresstensorxx >();
      arrayView2d< real32 const > const stressyy = elementSubRegion.getField< wavesolverfields::Stresstensoryy >();
      arrayView2d< real32 const > const stresszz = elementSubRegion.getField< wavesolverfields::Stresstensorzz >();
      arrayView2d< real32 const > const stressxy = elementSubRegion.getField< wavesolverfields::Stresstensorxy >();
      arrayView2d< real32 const > const stressxz = elementSubRegion.getField< wavesolverfields::Stresstensorxz >();
      arrayView2d< real32 const > const stressyz = elementSubRegion.getField< wavesolverfields::Stresstensoryz >();

      arrayView2d< real32 > const sigmaxxReceivers   = m_sigmaxxNp1AtReceivers.toView();
      arrayView2d< real32 > const sigmayyReceivers   = m_sigmayyNp1AtReceivers.toView();
      arrayView2d< real32 > const sigmazzReceivers   = m_sigmazzNp1AtReceivers.toView();
      arrayView2d< real32 > const sigmaxyReceivers   = m_sigmaxyNp1AtReceivers.toView();
      arrayView2d< real32 > const sigmaxzReceivers   = m_sigmaxzNp1AtReceivers.toView();
      arrayView2d< real32 > const sigmayzReceivers   = m_sigmayzNp1AtReceivers.toView();

      compute2dVariableAllSeismoTraces( regionIndex, time_n, 0, stressxx, stressxx, sigmaxxReceivers );
      compute2dVariableAllSeismoTraces( regionIndex, time_n, 0, stressyy, stressyy, sigmayyReceivers );
      compute2dVariableAllSeismoTraces( regionIndex, time_n, 0, stresszz, stresszz, sigmazzReceivers );
      compute2dVariableAllSeismoTraces( regionIndex, time_n, 0, stressxy, stressxy, sigmaxyReceivers );
      compute2dVariableAllSeismoTraces( regionIndex, time_n, 0, stressxz, stressxz, sigmaxzReceivers );
      compute2dVariableAllSeismoTraces( regionIndex, time_n, 0, stressyz, stressyz, sigmayzReceivers );
    } );
    arrayView1d< real32 > const ux_np1 = nodeManager.getField< wavesolverfields::Displacementx_np1 >();
    arrayView1d< real32 > const uy_np1 = nodeManager.getField< wavesolverfields::Displacementy_np1 >();
    arrayView1d< real32 > const uz_np1 = nodeManager.getField< wavesolverfields::Displacementz_np1 >();

    // compute the seismic traces since last step.
    arrayView2d< real32 > const uxReceivers   = m_displacementxNp1AtReceivers.toView();
    arrayView2d< real32 > const uyReceivers   = m_displacementyNp1AtReceivers.toView();
    arrayView2d< real32 > const uzReceivers   = m_displacementzNp1AtReceivers.toView();

    computeAllSeismoTraces( time_n, 0, ux_np1, ux_np1, uxReceivers );
    computeAllSeismoTraces( time_n, 0, uy_np1, uy_np1, uyReceivers );
    computeAllSeismoTraces( time_n, 0, uz_np1, uz_np1, uzReceivers );

  } );

// increment m_indexSeismoTrace
  while( (m_dtSeismoTrace*m_indexSeismoTrace) <= (time_n + epsilonLoc) && m_indexSeismoTrace < m_nsamplesSeismoTrace )
  {
    m_indexSeismoTrace++;
  }

}

void ElasticFirstOrderWaveEquationSEM::computeAllSeismoTraces( real64 const time_n,
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

void ElasticFirstOrderWaveEquationSEM::compute2dVariableAllSeismoTraces( localIndex const regionIndex,
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
    WaveSolverUtils::compute2dVariableSeismoTrace( time_n, dt, regionIndex, m_receiverRegion, timeSeismo, indexSeismoTrace, m_rcvElem, m_receiverConstants, m_receiverIsLocal,
                                                   m_nsamplesSeismoTrace, m_outputSeismoTrace, var_np1, var_n, varAtReceivers );
  }
}

void ElasticFirstOrderWaveEquationSEM::initializePML()
{
  GEOS_ERROR( "PML for the first order elastic wave propagator not yet implemented" );
}

void ElasticFirstOrderWaveEquationSEM::applyPML( real64 const, DomainPartition & )
{
  GEOS_ERROR( "PML for the first order elastic wave propagator not yet implemented" );
}

REGISTER_CATALOG_ENTRY( SolverBase, ElasticFirstOrderWaveEquationSEM, string const &, dataRepository::Group * const )

} /* namespace geos */
