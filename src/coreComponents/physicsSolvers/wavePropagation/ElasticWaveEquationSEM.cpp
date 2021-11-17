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
 * @file ElasticWaveEquationSEM.cpp
 */

#include "ElasticWaveEquationSEM.hpp"
#include "ElasticWaveEquationSEMKernel.hpp"

#include "dataRepository/KeyNames.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/CellBlock.hpp"
#include "mesh/ElementType.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

namespace geosx
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

  registerWrapper( viewKeyStruct::sourceConstantsString(), &m_sourceConstants ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Constant part of the source for the nodes listed in m_sourceNodeIds" );
  
  registerWrapper( viewKeyStruct::sourceIsLocalString(), &m_sourceIsLocal ).
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

  registerWrapper( viewKeyStruct::displacementNp1AtReceiversString(), &m_displacementNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Displacement value at each receiver for each timestep" );

}

ElasticWaveEquationSEM::~ElasticWaveEquationSEM()
{
  // TODO Auto-generated destructor stub
}

void ElasticWaveEquationSEM::initializePreSubGroups()
{

  WaveSolverBase::initializePreSubGroups();
  
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();

  FiniteElementDiscretizationManager const &
  feDiscretizationManager = numericalMethodManager.getFiniteElementDiscretizationManager();

  FiniteElementDiscretization const * const
  feDiscretization = feDiscretizationManager.getGroupPointer< FiniteElementDiscretization >( m_discretizationName );
  GEOSX_THROW_IF( feDiscretization == nullptr,
                  getName() << ": FE discretization not found: " << m_discretizationName,
                  InputError );
}


void ElasticWaveEquationSEM::registerDataOnMesh( Group & meshBodies )
{

  meshBodies.forSubGroups< MeshBody >( [&]( MeshBody & meshBody )
  {

    MeshLevel & meshLevel =  meshBody.getMeshLevel( 0 );

    NodeManager & nodeManager = meshLevel.getNodeManager();

    nodeManager.registerExtrinsicData< extrinsicMeshData::Displacementx_np1,
                                       extrinsicMeshData::Displacementy_np1,
                                       extrinsicMeshData::Displacementz_np1,
                                       extrinsicMeshData::ForcingRHS,
                                       extrinsicMeshData::MassVector,
                                       extrinsicMeshData::DampingVector_x,
                                       extrinsicMeshData::DampingVector_y,
                                       extrinsicMeshData::DampingVector_z,
                                       extrinsicMeshData::FreeSurfaceNodeIndicator >( this->getName() );

    FaceManager & faceManager = meshLevel.getFaceManager();
    faceManager.registerExtrinsicData< extrinsicMeshData::FreeSurfaceFaceIndicator >( this->getName() );

    ElementRegionManager & elemManager = meshLevel.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerExtrinsicData< extrinsicMeshData::MediumVelocityVp >( this->getName() );
      subRegion.registerExtrinsicData< extrinsicMeshData::MediumVelocityVs >( this->getName() );
      subRegion.registerExtrinsicData< extrinsicMeshData::MediumDensity >( this->getName() );

      subRegion.registerExtrinsicData< extrinsicMeshData::LameCoefficientLambda >(this->getName());
      subRegion.registerExtrinsicData< extrinsicMeshData::LameCoefficientMu >(this->getName());

      subRegion.registerExtrinsicData< extrinsicMeshData::Stresstensor_xx >(this->getName());
      subRegion.registerExtrinsicData< extrinsicMeshData::Stresstensor_yy >(this->getName());
      subRegion.registerExtrinsicData< extrinsicMeshData::Stresstensor_zz >(this->getName());
      subRegion.registerExtrinsicData< extrinsicMeshData::Stresstensor_xy >(this->getName());
      subRegion.registerExtrinsicData< extrinsicMeshData::Stresstensor_xz >(this->getName());
      subRegion.registerExtrinsicData< extrinsicMeshData::Stresstensor_yz >(this->getName());

      finiteElement::FiniteElementBase const & fe = subRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      finiteElement::dispatch3D( fe,
                                 [&]
                                   ( auto const finiteElement )
      {

       using FE_TYPE = TYPEOFREF( finiteElement );

       constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;

      subRegion.getExtrinsicData< extrinsicMeshData::Stresstensor_xx >().resizeDimension < 1 > ( numNodesPerElem ) ;
      subRegion.getExtrinsicData< extrinsicMeshData::Stresstensor_yy >().resizeDimension < 1 > ( numNodesPerElem ) ;
      subRegion.getExtrinsicData< extrinsicMeshData::Stresstensor_zz >().resizeDimension < 1 > ( numNodesPerElem ) ;
      subRegion.getExtrinsicData< extrinsicMeshData::Stresstensor_xy >().resizeDimension < 1 > ( numNodesPerElem ) ;
      subRegion.getExtrinsicData< extrinsicMeshData::Stresstensor_xz >().resizeDimension < 1 > ( numNodesPerElem ) ;
      subRegion.getExtrinsicData< extrinsicMeshData::Stresstensor_yz >().resizeDimension < 1 > ( numNodesPerElem ) ;
      } );


    } );

  } );
}



void ElasticWaveEquationSEM::postProcessInput()
{

  GEOSX_ERROR_IF( m_sourceCoordinates.size( 1 ) != 3,
                  "Invalid number of physical coordinates for the sources" );

  GEOSX_ERROR_IF( m_receiverCoordinates.size( 1 ) != 3,
                  "Invalid number of physical coordinates for the receivers" );

  localIndex const numNodesPerElem = 8;

  localIndex const numSourcesGlobal = m_sourceCoordinates.size( 0 );

  m_sourceNodeIds.resizeDimension< 0, 1 >( numSourcesGlobal, numNodesPerElem );
  m_sourceConstants.resizeDimension< 0, 1 >( numSourcesGlobal, numNodesPerElem );
  m_sourceIsLocal.resizeDimension< 0 >( numSourcesGlobal );

  localIndex const numReceiversGlobal = m_receiverCoordinates.size( 0 );
  m_receiverNodeIds.resizeDimension< 0, 1 >( numReceiversGlobal, numNodesPerElem );
  m_receiverConstants.resizeDimension< 0, 1 >( numReceiversGlobal, numNodesPerElem );
  m_receiverIsLocal.resizeDimension< 0 >( numReceiversGlobal );

  m_displacementNp1AtReceivers.resize( numReceiversGlobal );

}


void ElasticWaveEquationSEM::precomputeSourceAndReceiverTerm( MeshLevel & mesh )
{
  NodeManager & nodeManager = mesh.getNodeManager();
  FaceManager const & faceManager = mesh.getFaceManager();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();
  ArrayOfArraysView< localIndex const > const & facesToNodes =
    faceManager.nodeList().toViewConst();

  arrayView2d< real64 const > const sourceCoordinates = m_sourceCoordinates.toViewConst();
  arrayView2d< localIndex > const sourceNodeIds = m_sourceNodeIds.toView();
  arrayView2d< real64 > const sourceConstants = m_sourceConstants.toView();
  arrayView1d< localIndex > const sourceIsLocal = m_sourceIsLocal.toView();
  sourceNodeIds.setValues< serialPolicy >( -1 );
  sourceConstants.setValues< serialPolicy >( -1 );
  sourceIsLocal.setValues< serialPolicy >( 0 );

  arrayView2d< real64 const > const receiverCoordinates = m_receiverCoordinates.toViewConst();
  arrayView2d< localIndex > const receiverNodeIds = m_receiverNodeIds.toView();
  arrayView2d< real64 > const receiverConstants = m_receiverConstants.toView();
  arrayView1d< localIndex > const receiverIsLocal = m_receiverIsLocal.toView();
  receiverNodeIds.setValues< serialPolicy >( -1 );
  receiverConstants.setValues< serialPolicy >( -1 );
  receiverIsLocal.setValues< serialPolicy >( 0 );

  forTargetRegionsComplete( mesh, [&]( localIndex const,
                                       localIndex const,
                                       ElementRegionBase & elemRegion )
  {
    elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const,
                                                                       CellElementSubRegion & elementSubRegion )
    {

      GEOSX_THROW_IF( elementSubRegion.getElementType() != ElementType::Hexahedron,
                      "Invalid type of element, the elastic solver is designed for hexahedral meshes only (C3D8) ",
                      InputError );

      arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();
      arrayView2d< real64 const > const elemCenter = elementSubRegion.getElementCenter();

      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::dispatch3D( fe,
                                 [&]
                                   ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
        ElasticWaveEquationSEMKernels::
          PrecomputeSourceAndReceiverKernel::
          launch< EXEC_POLICY, FE_TYPE >
          ( elementSubRegion.size(),
          numNodesPerElem,
          X,
          elemsToNodes,
          elemsToFaces,
          facesToNodes,
          elemCenter,
          sourceCoordinates,
          sourceIsLocal,
          sourceNodeIds,
          sourceConstants,
          receiverCoordinates,
          receiverIsLocal,
          receiverNodeIds,
          receiverConstants );
      } );
    } );
  } );
}


void ElasticWaveEquationSEM::addSourceToRightHandSide( real64 const & time_n, arrayView1d< real64 > const rhs )
{
  arrayView2d< localIndex const > const sourceNodeIds = m_sourceNodeIds.toViewConst();
  arrayView2d< real64 const > const sourceConstants   = m_sourceConstants.toViewConst();
  arrayView1d< localIndex const > const sourceIsLocal = m_sourceIsLocal.toViewConst();

  real64 const fi = evaluateRicker( time_n, this->m_timeSourceFrequency, this->m_rickerOrder );

  
  forAll< serialPolicy >( m_sourceConstants.size( 0 ), [=] ( localIndex const isrc )
  {
    if( sourceIsLocal[isrc] == 1 )
    {
      for( localIndex inode = 0; inode < m_sourceConstants.size( 1 ); ++inode )
      {
        rhs[sourceNodeIds[isrc][inode]] = sourceConstants[isrc][inode] * fi;
      }
    }
  } );
}

void ElasticWaveEquationSEM::computeSeismoTrace( localIndex const iseismo, arrayView1d< real64 > const displacement_np1 )
{
  arrayView2d< localIndex const > const receiverNodeIds = m_receiverNodeIds.toViewConst();
  arrayView2d< real64 const > const receiverConstants   = m_receiverConstants.toViewConst();
  arrayView1d< localIndex const > const receiverIsLocal = m_receiverIsLocal.toViewConst();

  arrayView1d< real64 > const u_rcvs   = m_displacementNp1AtReceivers.toView();

  forAll< EXEC_POLICY >( receiverConstants.size( 0 ), [=] GEOSX_HOST_DEVICE ( localIndex const ircv )
  {
    if( receiverIsLocal[ircv] == 1 )
    {
      u_rcvs[ircv] = 0.0;
      for( localIndex inode = 0; inode < receiverConstants.size( 1 ); ++inode )
      {
        real64 const localIncrement = displacement_np1[receiverNodeIds[ircv][inode]] * receiverConstants[ircv][inode];
        RAJA::atomicAdd< ATOMIC_POLICY >( &u_rcvs[ircv], localIncrement );
      }
    }
  } );

  forAll< serialPolicy >( receiverConstants.size( 0 ), [=] ( localIndex const ircv )
  {
    if( this->m_outputSeismoTrace == 1 )
    {
      if( receiverIsLocal[ircv] == 1 )
      {
        // Note: this "manual" output to file is temporary
        //       It should be removed as soon as we can use TimeHistory to output data not registered on the mesh
        // TODO: remove saveSeismo and replace with TimeHistory
        this->saveSeismo( iseismo, u_rcvs[ircv], GEOSX_FMT( "seismoTraceReceiver{:03}.txt", ircv ) );
      }
    }
  } );
}

void ElasticWaveEquationSEM::saveSeismo( localIndex iseismo, real64 valDisplacement, string const & filename)
{
  std::ofstream f( filename, std::ios::app );
  f<< iseismo << " " << valDisplacement << std::endl;
  f.close();
}




void ElasticWaveEquationSEM::initializePostInitialConditionsPreSubGroups()
{

  WaveSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  //real64 const time = 0.0;
  // applyFreeSurfaceBC( time, domain );
  // applyABC( time, domain );
  precomputeSourceAndReceiverTerm( mesh );

  NodeManager & nodeManager = mesh.getNodeManager();
  FaceManager & faceManager = mesh.getFaceManager();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

  /// Get table containing all the face normals

  arrayView1d< real64 > const mass = nodeManager.getExtrinsicData< extrinsicMeshData::MassVector >();

  /// get array of indicators: 1 if face is on the free surface; 0 otherwise
  arrayView1d< localIndex const > const freeSurfaceFaceIndicator = faceManager.getExtrinsicData< extrinsicMeshData::FreeSurfaceFaceIndicator >();


  forTargetRegionsComplete( mesh, [&]( localIndex const,
                                       localIndex const,
                                       ElementRegionBase & elemRegion )
  {
    elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const,
                                                                       CellElementSubRegion & elementSubRegion )
    {

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();

      arrayView1d< real64 > const density = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumDensity>();

      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::dispatch3D( fe,
                                 [&]
                                   ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );
         
        ElasticWaveEquationSEMKernels::MassMatrixKernel< FE_TYPE > kernel( finiteElement );
        kernel.template launch< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),                                                            
                                                              X,
                                                              elemsToNodes,
                                                              density,
                                                              mass );
      } );
    } );
  } );
}


void ElasticWaveEquationSEM::applyFreeSurfaceBC( real64 const time, DomainPartition & domain )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  FaceManager & faceManager = domain.getMeshBody( 0 ).getMeshLevel( 0 ).getFaceManager();
  NodeManager & nodeManager = domain.getMeshBody( 0 ).getMeshLevel( 0 ).getNodeManager();

  arrayView1d< real64 > const ux_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementx_np1 >();
  arrayView1d< real64 > const uy_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementy_np1 >();
  arrayView1d< real64 > const uz_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementz_np1 >();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  /// set array of indicators: 1 if a face is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceFaceIndicator = faceManager.getExtrinsicData< extrinsicMeshData::FreeSurfaceFaceIndicator >();

  /// set array of indicators: 1 if a node is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceNodeIndicator = nodeManager.getExtrinsicData< extrinsicMeshData::FreeSurfaceNodeIndicator >();

  fsManager.apply( time,
                   domain,
                   "faceManager",
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

      freeSurfaceFaceIndicator.zero();
      freeSurfaceNodeIndicator.zero();

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

void ElasticWaveEquationSEM::applyABC( real64 const time, DomainPartition & domain )
{

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  NodeManager & nodeManager = mesh.getNodeManager();
  FaceManager & faceManager = mesh.getFaceManager();

  /// get the array of indicators: 1 if the face is on the boundary; 0 otherwise
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();
  
  /// Get table containing all the face normals
  arrayView2d< real64 const > const faceNormal  = faceManager.faceNormal();

  FaceManager::ElemMapType const & faceToElem = faceManager.toElementRelation();
  arrayView2d< localIndex const > const & faceToElemIndex = faceToElem.m_toElementIndex;


  ArrayOfArraysView< localIndex const > const facesToNodes = faceManager.nodeList().toViewConst();

  /// damping matrix to be computed for each dof in the boundary of the mesh
  arrayView1d< real64 > const damping_x = nodeManager.getExtrinsicData< extrinsicMeshData::DampingVector_x >();
  arrayView1d< real64 > const damping_y = nodeManager.getExtrinsicData< extrinsicMeshData::DampingVector_y >();
  arrayView1d< real64 > const damping_z = nodeManager.getExtrinsicData< extrinsicMeshData::DampingVector_z >();

  damping_x.zero();
  damping_y.zero();
  damping_z.zero();

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  fsManager.apply( time,
                   domain,
                   "faceManager",
                   string( "ABC" ),
                   [&]( FieldSpecificationBase const & bc,
                        string const &,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group &,
                        string const & )
  {
    string const & functionName = bc.getFunctionName();

    if( functionName.empty() || functionManager.getGroup< FunctionBase >( functionName ).isFunctionOfTime() == 2 )
    {

      forTargetRegionsComplete( mesh, [&]( localIndex const,
                                       localIndex const,
                                       ElementRegionBase & elemRegion )
      {

        elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const,
                                                                       CellElementSubRegion & elementSubRegion )
        {

          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();

          arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();

          arrayView1d< real64 > const density = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumDensity>();
          arrayView1d< real64 > const velocityVp = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumVelocityVp >();
          arrayView1d< real64 > const velocityVs = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumVelocityVs >();
         

          finiteElement::FiniteElementBase const &
          fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

          finiteElement::dispatch3D( fe,
                                     [&]
                                     ( auto const finiteElement )
          {
            using FE_TYPE = TYPEOFREF( finiteElement );

           
             ElasticWaveEquationSEMKernels::DampingMatrixKernel< FE_TYPE > kernel( finiteElement );
              kernel.template launch< EXEC_POLICY, ATOMIC_POLICY >(X,
                                                                   elemsToNodes,
                                                                   targetSet,
                                                                   faceToElemIndex,
                                                                   facesToNodes,
                                                                   faceNormal,
                                                                   density,
                                                                   velocityVp,
                                                                   velocityVs,
                                                                   damping_x,
                                                                   damping_y,
                                                                   damping_z );
                                                            

          } );
        } );
      } );  
    }
  } );
}


real64 ElasticWaveEquationSEM::solverStep( real64 const & time_n,
                                            real64 const & dt,
                                            integer const cycleNumber,
                                            DomainPartition & domain )
{
  return explicitStep( time_n, dt, cycleNumber, domain );
}

real64 ElasticWaveEquationSEM::explicitStep( real64 const & time_n,
                                              real64 const & dt,
                                              integer const cycleNumber,
                                              DomainPartition & domain )
{

  GEOSX_MARK_FUNCTION;

  GEOSX_UNUSED_VAR( time_n, dt, cycleNumber );

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView1d< real64 const > const mass = nodeManager.getExtrinsicData< extrinsicMeshData::MassVector >();
  arrayView1d< real64 const > const damping_x = nodeManager.getExtrinsicData< extrinsicMeshData::DampingVector_x >();
  arrayView1d< real64 const > const damping_y = nodeManager.getExtrinsicData< extrinsicMeshData::DampingVector_y >();
  arrayView1d< real64 const > const damping_z = nodeManager.getExtrinsicData< extrinsicMeshData::DampingVector_z >();

  arrayView1d< real64 > const ux_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementx_np1 >();
  arrayView1d< real64 > const uy_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementy_np1 >();
  arrayView1d< real64 > const uz_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacementz_np1 >();

  /// get array of indicators: 1 if node on free surface; 0 otherwise
  arrayView1d< localIndex const > const freeSurfaceNodeIndicator = nodeManager.getExtrinsicData< extrinsicMeshData::FreeSurfaceNodeIndicator >();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

  arrayView1d< real64 > const rhs = nodeManager.getExtrinsicData< extrinsicMeshData::ForcingRHS >();

  addSourceToRightHandSide( time_n, rhs );

  forTargetRegionsComplete( mesh, [&]( localIndex const,
                                       localIndex const,
                                       ElementRegionBase & elemRegion )

  {
    elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const,
                                                                        CellElementSubRegion & elementSubRegion )
    {

      arrayView1d< real64 const > const velocityVp = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumVelocityVp >();
      arrayView1d< real64 const > const velocityVs = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumVelocityVs >();
      arrayView1d< real64 const > const density = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumDensity >();

      arrayView1d< real64 > const lambda = elementSubRegion.getExtrinsicData< extrinsicMeshData::LameCoefficientLambda >();
      arrayView1d< real64 > const mu = elementSubRegion.getExtrinsicData< extrinsicMeshData::LameCoefficientMu >();

      arrayView2d< real64 > const stressxx = elementSubRegion.getExtrinsicData< extrinsicMeshData::Stresstensor_xx >();
      arrayView2d< real64 > const stressyy = elementSubRegion.getExtrinsicData< extrinsicMeshData::Stresstensor_yy >();
      arrayView2d< real64 > const stresszz = elementSubRegion.getExtrinsicData< extrinsicMeshData::Stresstensor_zz >();
      arrayView2d< real64 > const stressxy = elementSubRegion.getExtrinsicData< extrinsicMeshData::Stresstensor_xy >();
      arrayView2d< real64 > const stressxz = elementSubRegion.getExtrinsicData< extrinsicMeshData::Stresstensor_xz >();
      arrayView2d< real64 > const stressyz = elementSubRegion.getExtrinsicData< extrinsicMeshData::Stresstensor_yz >();

      auto kernelFactory = ElasticWaveEquationSEMKernels::ExplicitElasticDisplacementSEMFactory( dt );
  
      finiteElement::
      regionBasedKernelApplication< EXEC_POLICY,
                                    constitutive::NullModel,
                                    CellElementSubRegion >( mesh,
                                                            targetRegionNames(),
                                                            getDiscretizationName(),
                                                            arrayView1d< string const >(),
                                                            kernelFactory );

      auto kernelFactory2 = ElasticWaveEquationSEMKernels::ExplicitElasticStressSEMFactory( mu,
                                                                                           lambda,
                                                                                           stressxx,
                                                                                           stressyy,
                                                                                           stresszz,
                                                                                           stressxy,
                                                                                           stressxz,
                                                                                           stressyz,
                                                                                           dt );

      finiteElement::
      regionBasedKernelApplication< EXEC_POLICY,
                                    constitutive::NullModel,
                                    CellElementSubRegion >( mesh,
                                                            targetRegionNames(),
                                                            getDiscretizationName(),
                                                            arrayView1d< string const >(),
                                                            kernelFactory2 );                                         
                      
      } );

    } );



  /// Synchronize pressure fields
  std::map< string, string_array > fieldNames;
  fieldNames["node"].emplace_back( "displacementx_np1" );
  fieldNames["node"].emplace_back( "displacementy_np1" );
  fieldNames["node"].emplace_back( "displacementz_np1" );
  fieldNames["elems"].emplace_back( "stresstensor_xx" );
  fieldNames["elems"].emplace_back( "stresstensor_yy" );
  fieldNames["elems"].emplace_back( "stresstensor_zz" );
  fieldNames["elems"].emplace_back( "stresstensor_xy" );
  fieldNames["elems"].emplace_back( "stresstensor_xz" );
  fieldNames["elems"].emplace_back( "stresstensor_yz" );

  CommunicationTools & syncFields = CommunicationTools::getInstance();
  syncFields.synchronizeFields( fieldNames,
                                domain.getMeshBody( 0 ).getMeshLevel( 0 ),
                                domain.getNeighbors(),
                                true );

  if( this->m_outputSeismoTrace == 1 )
  {
    computeSeismoTrace( cycleNumber, ux_np1 );
    //computeSeismoTrace( cycleNumber, uy_np1 );
    //computeSeismoTrace( cycleNumber, uz_np1 );
  }

  return dt;

}

REGISTER_CATALOG_ENTRY( SolverBase, ElasticWaveEquationSEM, string const &, dataRepository::Group * const )

} /* namespace geosx */
