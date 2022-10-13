/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
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

#include "dataRepository/KeyNames.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/ElementType.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

namespace geosx
{

using namespace dataRepository;

ElasticFirstOrderWaveEquationSEM::ElasticFirstOrderWaveEquationSEM( const std::string & name,
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

}

ElasticFirstOrderWaveEquationSEM::~ElasticFirstOrderWaveEquationSEM()
{
  // TODO Auto-generated destructor stub
}

localIndex ElasticFirstOrderWaveEquationSEM::getNumNodesPerElem()
{
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();

  FiniteElementDiscretizationManager const &
  feDiscretizationManager = numericalMethodManager.getFiniteElementDiscretizationManager();

  FiniteElementDiscretization const * const
  feDiscretization = feDiscretizationManager.getGroupPointer< FiniteElementDiscretization >( m_discretizationName );
  GEOSX_THROW_IF( feDiscretization == nullptr,
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

void ElasticFirstOrderWaveEquationSEM::initializePreSubGroups()
{

  WaveSolverBase::initializePreSubGroups();

  localIndex numNodesPerElem = getNumNodesPerElem();

  localIndex const numSourcesGlobal = m_sourceCoordinates.size( 0 );
  m_sourceNodeIds.resize( numSourcesGlobal, numNodesPerElem );
  m_sourceConstants.resize( numSourcesGlobal, numNodesPerElem );
  //m_sourceIsAccessible.resize( numSourcesGlobal );

  localIndex const numReceiversGlobal = m_receiverCoordinates.size( 0 );
  m_receiverNodeIds.resize( numReceiversGlobal, numNodesPerElem );
  m_receiverConstants.resize( numReceiversGlobal, numNodesPerElem );
  m_receiverIsLocal.resize( numReceiversGlobal );

  m_displacementxNp1AtReceivers.resizeDimension< 1 >( numReceiversGlobal );
  m_displacementyNp1AtReceivers.resizeDimension< 1 >( numReceiversGlobal );
  m_displacementzNp1AtReceivers.resizeDimension< 1 >( numReceiversGlobal );


}


void ElasticFirstOrderWaveEquationSEM::registerDataOnMesh( Group & meshBodies )
{

 forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                   MeshLevel & mesh,
                                                   arrayView1d< string const > const & )

  {
    NodeManager & nodeManager = mesh.getNodeManager();

    nodeManager.registerExtrinsicData< extrinsicMeshData::Displacementx_np1,
                                       extrinsicMeshData::Displacementy_np1,
                                       extrinsicMeshData::Displacementz_np1,
                                       extrinsicMeshData::ForcingRHS,
                                       extrinsicMeshData::MassVector,
                                       extrinsicMeshData::DampingVectorx,
                                       extrinsicMeshData::DampingVectory,
                                       extrinsicMeshData::DampingVectorz,
                                       extrinsicMeshData::FreeSurfaceNodeIndicator >( this->getName() );

    FaceManager & faceManager = mesh.getFaceManager();
    faceManager.registerExtrinsicData< extrinsicMeshData::FreeSurfaceFaceIndicator >( this->getName() );

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerExtrinsicData< extrinsicMeshData::MediumVelocityVp >( this->getName() );
      subRegion.registerExtrinsicData< extrinsicMeshData::MediumVelocityVs >( this->getName() );
      subRegion.registerExtrinsicData< extrinsicMeshData::MediumDensity >( this->getName() );
      subRegion.registerExtrinsicData< extrinsicMeshData::Lambda >( this->getName() );
      subRegion.registerExtrinsicData< extrinsicMeshData::Mu >( this->getName() );

      subRegion.registerExtrinsicData< extrinsicMeshData::Stresstensorxx >(this->getName());
      subRegion.registerExtrinsicData< extrinsicMeshData::Stresstensoryy >(this->getName());
      subRegion.registerExtrinsicData< extrinsicMeshData::Stresstensorzz >(this->getName());
      subRegion.registerExtrinsicData< extrinsicMeshData::Stresstensorxy >(this->getName());
      subRegion.registerExtrinsicData< extrinsicMeshData::Stresstensorxz >(this->getName());
      subRegion.registerExtrinsicData< extrinsicMeshData::Stresstensoryz >(this->getName());

      finiteElement::FiniteElementBase const & fe = subRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      finiteElement::dispatch3D( fe,
                                 [&]
                                   ( auto const finiteElement )
      {

       using FE_TYPE = TYPEOFREF( finiteElement );

       constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;

      subRegion.getExtrinsicData< extrinsicMeshData::Stresstensorxx >().resizeDimension < 1 > ( numNodesPerElem ) ;
      subRegion.getExtrinsicData< extrinsicMeshData::Stresstensoryy >().resizeDimension < 1 > ( numNodesPerElem ) ;
      subRegion.getExtrinsicData< extrinsicMeshData::Stresstensorzz >().resizeDimension < 1 > ( numNodesPerElem ) ;
      subRegion.getExtrinsicData< extrinsicMeshData::Stresstensorxy >().resizeDimension < 1 > ( numNodesPerElem ) ;
      subRegion.getExtrinsicData< extrinsicMeshData::Stresstensorxz >().resizeDimension < 1 > ( numNodesPerElem ) ;
      subRegion.getExtrinsicData< extrinsicMeshData::Stresstensoryz >().resizeDimension < 1 > ( numNodesPerElem ) ;
      } );


    } );

  } );
}



void ElasticFirstOrderWaveEquationSEM::postProcessInput()
{

  WaveSolverBase::postProcessInput();
  GEOSX_THROW_IF( m_sourceCoordinates.size( 1 ) != 3,
                  "Invalid number of physical coordinates for the sources",
                  InputError );

  GEOSX_THROW_IF( m_receiverCoordinates.size( 1 ) != 3,
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

  GEOSX_THROW_IF( dt < epsilonLoc*maxTime, "Value for dt: " << dt <<" is smaller than local threshold: " << epsilonLoc, std::runtime_error );

  if( m_dtSeismoTrace > 0 )
  {
    m_nsamplesSeismoTrace = int( maxTime / m_dtSeismoTrace) + 1;
  }
  else
  {
    m_nsamplesSeismoTrace = 0;
  }
  localIndex const nsamples = int(maxTime/dt) + 1;

  localIndex const numNodesPerElem = 8;

  localIndex const numSourcesGlobal = m_sourceCoordinates.size( 0 );
  m_sourceNodeIds.resize( numSourcesGlobal, numNodesPerElem );
  m_sourceConstants.resize( numSourcesGlobal, numNodesPerElem );
  m_sourceIsAccessible.resize( numSourcesGlobal );
  m_sourceElem.resize(numSourcesGlobal);

  localIndex const numReceiversGlobal = m_receiverCoordinates.size( 0 );
  m_receiverNodeIds.resize( numReceiversGlobal, numNodesPerElem );
  m_receiverConstants.resize( numReceiversGlobal, numNodesPerElem );
  m_receiverIsLocal.resize( numReceiversGlobal );

  m_displacementxNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );
  m_displacementyNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );
  m_displacementzNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );
  m_sourceValue.resize( nsamples, numSourcesGlobal );

}


void ElasticFirstOrderWaveEquationSEM::precomputeSourceAndReceiverTerm( MeshLevel & mesh, arrayView1d< string const > const & regionNames )
{
  NodeManager const & nodeManager = mesh.getNodeManager();
  FaceManager const & faceManager = mesh.getFaceManager();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();
  arrayView2d< real64 const > const faceNormal  = faceManager.faceNormal();
  arrayView2d< real64 const > const faceCenter  = faceManager.faceCenter();

  arrayView2d< real64 const > const sourceCoordinates = m_sourceCoordinates.toViewConst();
  arrayView2d< localIndex > const sourceNodeIds = m_sourceNodeIds.toView();
  arrayView2d< real64 > const sourceConstants = m_sourceConstants.toView();
  arrayView1d< localIndex > const sourceIsAccessible = m_sourceIsAccessible.toView();
  arrayView1d< localIndex > const sourceElem = m_sourceElem.toView();
  sourceNodeIds.setValues< serialPolicy >( -1 );
  sourceConstants.setValues< serialPolicy >( -1 );
  sourceIsAccessible.zero();

  arrayView2d< real64 const > const receiverCoordinates = m_receiverCoordinates.toViewConst();
  arrayView2d< localIndex > const receiverNodeIds = m_receiverNodeIds.toView();
  arrayView2d< real64 > const receiverConstants = m_receiverConstants.toView();
  arrayView1d< localIndex > const receiverIsLocal = m_receiverIsLocal.toView();
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

  mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                     CellElementSubRegion & elementSubRegion )
  {

    GEOSX_THROW_IF( elementSubRegion.getElementType() != ElementType::Hexahedron,
                    "Invalid type of element, the elastic solver is designed for hexahedral meshes only (C3D8) ",
                    InputError );

    arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();
    arrayView2d< real64 const > const elemCenter = elementSubRegion.getElementCenter();
    arrayView1d< integer const > const elemGhostRank = elementSubRegion.ghostRank();

    finiteElement::FiniteElementBase const &
    fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
    finiteElement::dispatch3D( fe,
                               [&]
                                 ( auto const finiteElement )
    {
      using FE_TYPE = TYPEOFREF( finiteElement );

      localIndex const numFacesPerElem = elementSubRegion.numFacesPerElement();

      elasticFirstOrderWaveEquationSEMKernels::
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
          sourceElem,
          sourceNodeIds,
          sourceConstants,
          receiverCoordinates,
          receiverIsLocal,
          receiverNodeIds,
          receiverConstants,
          sourceValue,
          dt,
          timeSourceFrequency,
          rickerOrder);
    } );
  } );
}


void ElasticFirstOrderWaveEquationSEM::addSourceToRightHandSide( integer const & cycleNumber, arrayView1d< real32 > const rhs )
{
  arrayView2d< localIndex const > const sourceNodeIds = m_sourceNodeIds.toViewConst();
  arrayView2d< real64 const > const sourceConstants   = m_sourceConstants.toViewConst();
  arrayView1d< localIndex const > const sourceIsAccessible = m_sourceIsAccessible.toViewConst();
  arrayView2d< real32 const > const sourceValue   = m_sourceValue.toViewConst();
 
  GEOSX_THROW_IF( cycleNumber > sourceValue.size( 0 ), "Too many steps compared to array size", std::runtime_error );
  
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

void ElasticFirstOrderWaveEquationSEM::computeSeismoTrace( real64 const time_n,
                                                  real64 const dt,
                                                  real64 const timeSeismo,
                                                  localIndex iSeismo,
                                                  arrayView1d< real32 const > const var_np1,
                                                  arrayView1d< real32 const > const var_n,
                                                  arrayView2d< real32 > varAtReceivers )
{
  real64 const time_np1 = time_n+dt;
  arrayView2d< localIndex const > const receiverNodeIds = m_receiverNodeIds.toViewConst();
  arrayView2d< real64 const > const receiverConstants   = m_receiverConstants.toViewConst();
  arrayView1d< localIndex const > const receiverIsLocal = m_receiverIsLocal.toViewConst();

  real64 const a1 = (dt < epsilonLoc) ? 1.0 : (time_np1 - timeSeismo)/dt;
  real64 const a2 = 1.0 - a1;

  if( m_nsamplesSeismoTrace > 0 )
  {
    forAll< EXEC_POLICY >( receiverConstants.size( 0 ), [=] GEOSX_HOST_DEVICE ( localIndex const ircv )
    {
      if( receiverIsLocal[ircv] == 1 )
      {
        varAtReceivers[iSeismo][ircv] = 0.0;
        real32 vtmp_np1 = 0.0;
        real32 vtmp_n = 0.0;
        for( localIndex inode = 0; inode < receiverConstants.size( 1 ); ++inode )
        {
          vtmp_np1 += var_np1[receiverNodeIds[ircv][inode]] * receiverConstants[ircv][inode];
          vtmp_n += var_n[receiverNodeIds[ircv][inode]] * receiverConstants[ircv][inode];
        }
        // linear interpolation between the pressure value at time_n and time_(n+1)
        varAtReceivers[iSeismo][ircv] = a1*vtmp_n + a2*vtmp_np1;
      }
    } );
  }

  // TODO DEBUG: the following output is only temporary until our wave propagation kernels are finalized.
  // Output will then only be done via the previous code.
  if( iSeismo == m_nsamplesSeismoTrace - 1 )
  {
    forAll< serialPolicy >( receiverConstants.size( 0 ), [=] ( localIndex const ircv )
    {
      if( this->m_outputSeismoTrace == 1 )
      {
        if( receiverIsLocal[ircv] == 1 )
        {
          // Note: this "manual" output to file is temporary
          //       It should be removed as soon as we can use TimeHistory to output data not registered on the mesh
          // TODO: remove saveSeismo and replace with TimeHistory
          for( localIndex iSample = 0; iSample < m_nsamplesSeismoTrace; ++iSample )
          {
            this->saveSeismo( iSample, varAtReceivers[iSample][ircv], GEOSX_FMT( "seismoTraceReceiver{:03}.txt", ircv ) );
          }
        }
      }
    } );
  }

}

void ElasticFirstOrderWaveEquationSEM::saveSeismo( localIndex iseismo, real32 val, string const & filename)
{
  std::ofstream f( filename, std::ios::app );
  f<< iseismo << " " << val << std::endl;
  f.close();
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

    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

    /// Get table containing all the face normals
    arrayView2d< real64 const > const faceNormal  = faceManager.faceNormal();

    arrayView1d< integer > const & facesDomainBoundaryIndicator = faceManager.getDomainBoundaryIndicator();
    arrayView1d< localIndex const > const freeSurfaceFaceIndicator = faceManager.getExtrinsicData< extrinsicMeshData::FreeSurfaceFaceIndicator >();

    ArrayOfArraysView< localIndex const > const facesToNodes = faceManager.nodeList().toViewConst();

    arrayView1d< real32 > const mass = nodeManager.getExtrinsicData< extrinsicMeshData::MassVector >();

    arrayView1d< real32 > const dampingx = nodeManager.getExtrinsicData< extrinsicMeshData::DampingVectorx >();
    arrayView1d< real32 > const dampingy = nodeManager.getExtrinsicData< extrinsicMeshData::DampingVectory >();
    arrayView1d< real32 > const dampingz = nodeManager.getExtrinsicData< extrinsicMeshData::DampingVectorz >();

    dampingx.zero();
    dampingy.zero();
    dampingz.zero();

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          CellElementSubRegion & elementSubRegion )
    {

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();
      arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();

      arrayView1d< real32 > const density = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumDensity >();
      arrayView1d< real32 > const velocityVp = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumVelocityVp >();
      arrayView1d< real32 > const velocityVs = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumVelocityVs >();

      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::dispatch3D( fe,
                                 [&]
                                   ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        localIndex const numFacesPerElem = elementSubRegion.numFacesPerElement();
        localIndex const numNodesPerFace = facesToNodes.sizeOfArray( 0 );

        elasticFirstOrderWaveEquationSEMKernels::MassAndDampingMatrixKernel< FE_TYPE > kernel( finiteElement );
        kernel.template launch< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                              numFacesPerElem,
                                                              numNodesPerFace,
                                                              X,
                                                              elemsToNodes,
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
                                                              dampingz,
                                                              mass );
      } );
    } );
  } );
}


void ElasticFirstOrderWaveEquationSEM::applyFreeSurfaceBC( real64 const time, DomainPartition & domain )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  FaceManager & faceManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName).getFaceManager();
  NodeManager & nodeManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getNodeManager();

  arrayView1d< real32 > const ux_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementx_np1 >();
  arrayView1d< real32 > const uy_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementy_np1 >();
  arrayView1d< real32 > const uz_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementz_np1 >();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  /// set array of indicators: 1 if a face is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceFaceIndicator = faceManager.getExtrinsicData< extrinsicMeshData::FreeSurfaceFaceIndicator >();

  /// set array of indicators: 1 if a node is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceNodeIndicator = nodeManager.getExtrinsicData< extrinsicMeshData::FreeSurfaceNodeIndicator >();

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
      GEOSX_ERROR( "This option is not supported yet" );
    }
  } );
}

real64 ElasticFirstOrderWaveEquationSEM::solverStep( real64 const & time_n,
                                            real64 const & dt,
                                            integer const cycleNumber,
                                            DomainPartition & domain )
{
  return explicitStep( time_n, dt, cycleNumber, domain );
}

real64 ElasticFirstOrderWaveEquationSEM::explicitStep( real64 const & time_n,
                                              real64 const & dt,
                                              integer const cycleNumber,
                                              DomainPartition & domain )
{

  GEOSX_MARK_FUNCTION;

  GEOSX_UNUSED_VAR( time_n, dt, cycleNumber );

   
  arrayView2d< real64 const > const sourceConstants = m_sourceConstants.toView();
  arrayView1d< localIndex const > const sourceIsAccessible = m_sourceIsAccessible.toView();
  arrayView1d< localIndex const > const sourceElem = m_sourceElem.toView();
  arrayView2d< real32 const > const sourceValue = m_sourceValue.toView();


 forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )

  {

    NodeManager & nodeManager = mesh.getNodeManager();
  
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

    arrayView1d< real32 const > const mass = nodeManager.getExtrinsicData< extrinsicMeshData::MassVector >();
    // arrayView1d< real64 const > const dampingx = nodeManager.getExtrinsicData< extrinsicMeshData::DampingVectorx >();
    // arrayView1d< real64 const > const dampingy = nodeManager.getExtrinsicData< extrinsicMeshData::DampingVectory >();
    // arrayView1d< real64 const > const dampingz = nodeManager.getExtrinsicData< extrinsicMeshData::DampingVectorz >();

    arrayView1d< real32 > const ux_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementx_np1 >();
    arrayView1d< real32 > const uy_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementy_np1 >();
    arrayView1d< real32 > const uz_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementz_np1 >();

    /// get array of indicators: 1 if node on free surface; 0 otherwise
    //arrayView1d< localIndex const > const freeSurfaceNodeIndicator = nodeManager.getExtrinsicData< extrinsicMeshData::FreeSurfaceNodeIndicator >();

    //arrayView1d< real64 > const rhs = nodeManager.getExtrinsicData< extrinsicMeshData::ForcingRHS >();

  
   mesh.getElemManager().forElementSubRegions < CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          CellElementSubRegion & elementSubRegion )
    {

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();

      arrayView1d< real32 const > const velocityVp = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumVelocityVp >();
      arrayView1d< real32 const > const velocityVs = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumVelocityVs >();
      arrayView1d< real32 const > const density = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumDensity >();

      arrayView1d< real32 > const lambda = elementSubRegion.getExtrinsicData< extrinsicMeshData::Lambda >();
      arrayView1d< real32 > const mu = elementSubRegion.getExtrinsicData< extrinsicMeshData::Mu >();

      arrayView2d< real32 > const stressxx = elementSubRegion.getExtrinsicData< extrinsicMeshData::Stresstensorxx >();
      arrayView2d< real32 > const stressyy = elementSubRegion.getExtrinsicData< extrinsicMeshData::Stresstensoryy >();
      arrayView2d< real32 > const stresszz = elementSubRegion.getExtrinsicData< extrinsicMeshData::Stresstensorzz >();
      arrayView2d< real32 > const stressxy = elementSubRegion.getExtrinsicData< extrinsicMeshData::Stresstensorxy >();
      arrayView2d< real32 > const stressxz = elementSubRegion.getExtrinsicData< extrinsicMeshData::Stresstensorxz >();
      arrayView2d< real32 > const stressyz = elementSubRegion.getExtrinsicData< extrinsicMeshData::Stresstensoryz >();

      //addSourceToRightHandSide( cycleNumber, rhs );

      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::dispatch3D( fe,
                                 [&]
                                      ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        elasticFirstOrderWaveEquationSEMKernels::
          StressComputation< FE_TYPE > kernel( finiteElement );
        kernel.template launch< EXEC_POLICY, ATOMIC_POLICY >
        (elementSubRegion.size(),
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
         sourceValue,
         dt,
         cycleNumber,
         stressxx,
         stressyy,
         stresszz,
         stressxy,
         stressxz,
         stressyz);

         elasticFirstOrderWaveEquationSEMKernels::
           VelocityComputation< FE_TYPE > kernel2( finiteElement );
         kernel2.template launch< EXEC_POLICY, ATOMIC_POLICY >
         (elementSubRegion.size(),
          X,
          elemsToNodes,
          stressxx,
          stressyy,
          stresszz,
          stressxy,
          stressxz,
          stressyz,
          mass,
          dt,
          ux_np1,
          uy_np1,
          uz_np1);
        

      } );

      // auto kernelFactory = elasticFirstOrderWaveEquationSEMKernels::ExplicitElasticDisplacementSEMFactory( dt );
  
      // finiteElement::
      // regionBasedKernelApplication< EXEC_POLICY,
      //                               constitutive::NullModel,
      //                               CellElementSubRegion >( mesh,
      //                                                       regionNames,
      //                                                       getDiscretizationName(),
      //                                                       "",
      //                                                       kernelFactory );

      // auto kernelFactory2 = elasticFirstOrderWaveEquationSEMKernels::ExplicitElasticStressSEMFactory(stressxx,
      //                                                                                      stressyy,
      //                                                                                      stresszz,
      //                                                                                      stressxy,
      //                                                                                      stressxz,
      //                                                                                      stressyz,
      //                                                                                      sourceElem,
      //                                                                                      sourceIsAccessible,
      //                                                                                      sourceConstants,
      //                                                                                      sourceValue,
      //                                                                                      dt,
      //                                                                                      cycleNumber );

      // finiteElement::
      // regionBasedKernelApplication< EXEC_POLICY,
      //                               constitutive::NullModel,
      //                               CellElementSubRegion >( mesh,
      //                                                       regionNames,
      //                                                       getDiscretizationName(),
      //                                                       "",
      //                                                       kernelFactory2 );                                         
                      
      } );

      FieldIdentifiers fieldsToBeSync;
      fieldsToBeSync.addFields( FieldLocation::Node, { extrinsicMeshData::Displacementx_np1::key(), extrinsicMeshData::Displacementy_np1::key(),  extrinsicMeshData::Displacementz_np1::key()} );
      fieldsToBeSync.addElementFields( {extrinsicMeshData::Stresstensorxx::key(), extrinsicMeshData::Stresstensoryy::key(), extrinsicMeshData::Stresstensorzz::key(), extrinsicMeshData::Stresstensorxy::key(), 
                                        extrinsicMeshData::Stresstensorxz::key(), extrinsicMeshData::Stresstensoryz::key()}, regionNames );


      CommunicationTools & syncFields = CommunicationTools::getInstance();
      syncFields.synchronizeFields( fieldsToBeSync,
                                domain.getMeshBody( 0 ).getMeshLevel(  m_discretizationName),
                                domain.getNeighbors(),
                                true );

      // compute the seismic traces since last step.
      arrayView2d< real32 > const uxReceivers   = m_displacementxNp1AtReceivers.toView();
      arrayView2d< real32 > const uyReceivers   = m_displacementyNp1AtReceivers.toView();
      arrayView2d< real32 > const uzReceivers   = m_displacementzNp1AtReceivers.toView();

      computeAllSeismoTraces( time_n, dt, ux_np1, ux_np1, uxReceivers );
      computeAllSeismoTraces( time_n, dt, uy_np1, uy_np1, uxReceivers );
      computeAllSeismoTraces( time_n, dt, uz_np1, uz_np1, uxReceivers );

    } );
  
  return dt;

}

void ElasticFirstOrderWaveEquationSEM::cleanup( real64 const time_n,
                                                integer const cycleNumber,
                                                integer const eventCounter,
                                                real64 const eventProgress,
                                                DomainPartition & domain )
{
  // compute the remaining seismic traces, if needed
 forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< real32 > const ux_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementx_np1 >();
    arrayView1d< real32 > const uy_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementy_np1 >();
    arrayView1d< real32 > const uz_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementz_np1 >();
    // compute the seismic traces since last step.
    arrayView2d< real32 > const uxReceivers   = m_displacementxNp1AtReceivers.toView();
    arrayView2d< real32 > const uyReceivers   = m_displacementyNp1AtReceivers.toView();
    arrayView2d< real32 > const uzReceivers   = m_displacementzNp1AtReceivers.toView();

    computeAllSeismoTraces( time_n, 0, ux_np1, ux_np1, uxReceivers );
    computeAllSeismoTraces( time_n, 0, uy_np1, uy_np1, uxReceivers );
    computeAllSeismoTraces( time_n, 0, uz_np1, uz_np1, uxReceivers );

  } );
}

void ElasticFirstOrderWaveEquationSEM::computeAllSeismoTraces( real64 const time_n,
                                                      real64 const dt,
                                                      arrayView1d< real32 const > const var_np1,
                                                      arrayView1d< real32 const > const var_n,
                                                      arrayView2d< real32 > varAtReceivers )
{
  for( real64 timeSeismo;
       (timeSeismo = m_dtSeismoTrace*m_indexSeismoTrace) <= (time_n + epsilonLoc) && m_indexSeismoTrace < m_nsamplesSeismoTrace;
       m_indexSeismoTrace++ )
  {
    computeSeismoTrace( time_n, dt, timeSeismo, m_indexSeismoTrace, var_np1, var_n, varAtReceivers );
  }
}

void ElasticFirstOrderWaveEquationSEM::initializePML()
{
  GEOSX_ERROR( "PML for the first order elastic wave propagator not yet implemented" );
}

void ElasticFirstOrderWaveEquationSEM::applyPML( real64 const, DomainPartition & )
{
  GEOSX_ERROR( "PML for the first order elastic wave propagator not yet implemented" );
}

REGISTER_CATALOG_ENTRY( SolverBase, ElasticFirstOrderWaveEquationSEM, string const &, dataRepository::Group * const )

} /* namespace geosx */
