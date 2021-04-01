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

#include "dataRepository/KeyNames.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "managers/ProblemManager.hpp"
#include "mpiCommunications/CommunicationTools.hpp"

namespace geosx
{

using namespace dataRepository;

ElasticWaveEquationSEM::ElasticWaveEquationSEM( const std::string & name,
                                                  Group * const parent ):
  SolverBase( name,
              parent )
{

  registerWrapper( viewKeyStruct::sourceCoordinatesString(), &m_sourceCoordinates ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Coordinates (x,y,z) of the sources" );

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


  registerWrapper( viewKeyStruct::timeSourceFrequencyString(), &m_timeSourceFrequency ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Central frequency for the time source" );


  registerWrapper( viewKeyStruct::receiverCoordinatesString(), &m_receiverCoordinates ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Coordinates (x,y,z) of the receivers" );

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

    registerWrapper( viewKeyStruct::rickerOrderString(), &m_rickerOrder ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 2 ).
    setDescription( "Flag that indicates the order of the Ricker to be used o, 1 or 2. Order 2 by default" );

     registerWrapper( viewKeyStruct::outputSismoTraceString(), &m_outputSismoTrace ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag that indicates if we write the sismo trace in a file .txt, 0 no output, 1 otherwise" );

  registerWrapper( viewKeyStruct::displacementNp1AtReceiversString(), &m_displacementNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Displacement value at each receiver for each timestep" );


}

ElasticWaveEquationSEM::~ElasticWaveEquationSEM()
{
  // TODO Auto-generated destructor stub
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

  m_displacementNp1AtReceivers.resizeDimension< 0 >( numReceiversGlobal );

}

void ElasticWaveEquationSEM::registerDataOnMesh( Group & meshBodies )
{

  meshBodies.forSubGroups< MeshBody >( [&]( MeshBody & meshBody )
  {

    MeshLevel & meshLevel =  meshBody.getMeshLevel( 0 );

    NodeManager & nodeManager = meshLevel.getNodeManager();

    nodeManager.registerExtrinsicData< extrinsicMeshData::Displacement_nm1,
                                       extrinsicMeshData::Displacement_n,
                                       extrinsicMeshData::Displacement_np1,
                                       extrinsicMeshData::ForcingRHS,
                                       extrinsicMeshData::MassVector,
                                       extrinsicMeshData::DampingVector,
                                       extrinsicMeshData::StiffnessVector,
                                       extrinsicMeshData::FreeSurfaceNodeIndicator >( this->getName() );

    FaceManager & faceManager = meshLevel.getFaceManager();
    faceManager.registerExtrinsicData< extrinsicMeshData::FreeSurfaceFaceIndicator >( this->getName() );

    ElementRegionManager & elemManager = meshLevel.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerExtrinsicData< extrinsicMeshData::MediumVelocityVp >( this->getName() );
      subRegion.registerExtrinsicData< extrinsicMeshData::MediumVelocityVs >( this->getName() );
      subRegion.registerExtrinsicData< extrinsicMeshData::MediumDensity >( this->getName() );
    } );

  } );
}


void ElasticWaveEquationSEM::precomputeSourceAndReceiverTerm( MeshLevel & mesh )
{
  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

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

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();

      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::dispatch3D( fe,
                                 [&]
                                   ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
        localIndex const numFacesPerElem = elementSubRegion.numFacesPerElement();
        array1d< array1d< localIndex > > faceNodes( numFacesPerElem );

        for( localIndex k = 0; k < elementSubRegion.size(); ++k )
        {

          for( localIndex kf = 0; kf < numFacesPerElem; ++kf )
          {
            elementSubRegion.getFaceNodes( k, kf, faceNodes[kf] );
          }

          /// loop over all the source that haven't been found yet
          for( localIndex isrc = 0; isrc < sourceCoordinates.size( 0 ); ++isrc )
          {
            if( sourceIsLocal[isrc] == 0 )
            {
              real64 const coords[3] = { sourceCoordinates[0][isrc],
                                         sourceCoordinates[1][isrc],
                                         sourceCoordinates[2][isrc] };

              if( computationalGeometry::IsPointInsidePolyhedron( X, faceNodes, coords ) )
              {
                sourceIsLocal[isrc] = 1;

                real64 xLocal[numNodesPerElem][3];
                for( localIndex a=0; a< numNodesPerElem; ++a )
                {
                  for( localIndex i=0; i<3; ++i )
                  {
                    xLocal[a][i] = X( elemsToNodes( k, a ), i );
                  }
                }

                /// coordsOnRefElem = invJ*(coords-coordsNode_0)
                real64 coordsOnRefElem[3];
                localIndex q=0;

                real64 invJ[3][3]={{0}};
                FE_TYPE::invJacobianTransformation( q, xLocal, invJ );

                real64 coordsRef[3]={0};
                for( localIndex i=0; i<3; ++i )
                {
                  coordsRef[i] = coords[i] - xLocal[q][i];
                }

                for( localIndex i=0; i<3; ++i )
                {
                  // Init at (-1,-1,-1) as the origin of the referential elem
                  coordsOnRefElem[i] =-1.0;
                  for( localIndex j=0; j<3; ++j )
                  {
                    coordsOnRefElem[i] += invJ[i][j]*coordsRef[j];
                  }
                }

                real64 Ntest[8];
                finiteElement::LagrangeBasis1::TensorProduct3D::value( coordsOnRefElem, Ntest );

                for( localIndex a=0; a< numNodesPerElem; ++a )
                {
                  sourceNodeIds[isrc][a] = elemsToNodes[k][a];
                  sourceConstants[isrc][a] = Ntest[a];
                }
              }
            }
          } // End loop over all source


          /// loop over all the receiver that haven't been found yet
          for( localIndex ircv = 0; ircv < receiverCoordinates.size( 0 ); ++ircv )
          {
            if( receiverIsLocal[ircv] == 0 )
            {
              real64 const coords[3] = { receiverCoordinates[0][ircv],
                                         receiverCoordinates[1][ircv],
                                         receiverCoordinates[2][ircv] };

              if( computationalGeometry::IsPointInsidePolyhedron( X, faceNodes, coords ) )
              {
                receiverIsLocal[ircv] = 1;

                real64 xLocal[numNodesPerElem][3];
                for( localIndex a=0; a< numNodesPerElem; ++a )
                {
                  for( localIndex i=0; i<3; ++i )
                  {
                    xLocal[a][i] = X( elemsToNodes( k, a ), i );
                  }
                }

                real64 coordsOnRefElem[3];
                localIndex q=0;

                real64 invJ[3][3]={{0}};
                FE_TYPE::invJacobianTransformation( q, xLocal, invJ );

                real64 coordsRef[3]={0};
                for( localIndex i=0; i<3; ++i )
                {
                  coordsRef[i] = coords[i] - xLocal[q][i];
                }

                for( localIndex i=0; i<3; ++i )
                {
                  /// Init at (-1,-1,-1) as the origin of the referential elem
                  coordsOnRefElem[i] =-1.0;
                  for( localIndex j=0; j<3; ++j )
                  {
                    coordsOnRefElem[i] += invJ[i][j]*coordsRef[j];
                  }
                }

                real64 Ntest[8];
                finiteElement::LagrangeBasis1::TensorProduct3D::value( coordsOnRefElem, Ntest );

                for( localIndex a=0; a< numNodesPerElem; ++a )
                {
                  receiverNodeIds[ircv][a] = elemsToNodes[k][a];
                  receiverConstants[ircv][a] = Ntest[a];
                }
              }
            }
          } // End loop over reciever

        } // End loop over elements
      } );
    } );
  } );
}


void ElasticWaveEquationSEM::addSourceToRightHandSide( real64 const & time, arrayView2d< real64 > const rhs )
{
  arrayView2d< localIndex const > const sourceNodeIds = m_sourceNodeIds.toViewConst();
  arrayView2d< real64 const > const sourceConstants   = m_sourceConstants.toViewConst();
  arrayView1d< localIndex const > const sourceIsLocal = m_sourceIsLocal.toViewConst();

  real64 const fi = evaluateRickerOrder2( time, this->m_timeSourceFrequency );

  for( localIndex isrc = 0; isrc < sourceConstants.size( 0 ); ++isrc )
  {
    if( sourceIsLocal[isrc] == 1 )
    {
      for( localIndex inode = 0; inode < sourceConstants.size( 1 ); ++inode )
      {
        rhs[0][sourceNodeIds[isrc][inode]] = sourceConstants[isrc][inode] * fi;
        rhs[1][sourceNodeIds[isrc][inode]] = sourceConstants[isrc][inode] * fi;
        rhs[2][sourceNodeIds[isrc][inode]] = sourceConstants[isrc][inode] * fi;
      }
    }
  }
}


void ElasticWaveEquationSEM::computeSismoTrace( localIndex const isismo, arrayView2d< real64 > const u_np1 )
{
  arrayView2d< localIndex const > const receiverNodeIds = m_receiverNodeIds.toViewConst();
  arrayView2d< real64 const > const receiverConstants   = m_receiverConstants.toViewConst();
  arrayView1d< localIndex const > const receiverIsLocal = m_receiverIsLocal.toViewConst();

  arrayView2d< real64 > const u_rcvs   = m_displacementNp1AtReceivers.toView();


  char filename[50];

  for( localIndex ircv = 0; ircv < receiverConstants.size( 0 ); ++ircv )
  {
    if( receiverIsLocal[ircv] == 1 )
    {
      u_rcvs[0][ircv] = 0.0;
      u_rcvs[1][ircv] = 0.0;
      u_rcvs[2][ircv] = 0.0;

      for( localIndex inode = 0; inode < receiverConstants.size( 1 ); ++inode )
      {
        u_rcvs[0][ircv] += u_np1[0][receiverNodeIds[ircv][inode]]*receiverConstants[ircv][inode];
        u_rcvs[1][ircv] += u_np1[1][receiverNodeIds[ircv][inode]]*receiverConstants[ircv][inode];
        u_rcvs[2][ircv] += u_np1[2][receiverNodeIds[ircv][inode]]*receiverConstants[ircv][inode];
      }

      sprintf( filename, "sismoTraceUxReceiver%0ld.txt", ircv );
      this->saveSismo( isismo, u_rcvs[0][ircv], filename );
      sprintf( filename, "sismoTraceUyReceiver%0ld.txt", ircv );
      this->saveSismo( isismo, u_rcvs[1][ircv], filename );
      sprintf( filename, "sismoTraceUzReceiver%0ld.txt", ircv );
      this->saveSismo( isismo, u_rcvs[2][ircv], filename );
    }
  }
}


void ElasticWaveEquationSEM::saveSismo( localIndex isismo, real64 val_displacement, char *filename )
{
  std::ofstream f( filename, std::ios::app );
  f<< isismo << " " << val_displacement << std::endl;
  f.close();
}


void ElasticWaveEquationSEM::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  DomainPartition & domain = getGlobalState().getProblemManager().getDomainPartition();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();

  FiniteElementDiscretizationManager const &
  feDiscretizationManager = numericalMethodManager.getFiniteElementDiscretizationManager();

  FiniteElementDiscretization const * const
  feDiscretization = feDiscretizationManager.getGroupPointer< FiniteElementDiscretization >( m_discretizationName );
  GEOSX_ERROR_IF( feDiscretization == nullptr, getName() << ": FE discretization not found: " << m_discretizationName );
}


void ElasticWaveEquationSEM::initializePostInitialConditionsPreSubGroups()
{
  DomainPartition & domain = getGlobalState().getProblemManager().getDomainPartition();
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  real64 const time = 0.0;
  applyFreeSurfaceBC( time, domain );
  precomputeSourceAndReceiverTerm( mesh );

  NodeManager & nodeManager = mesh.getNodeManager();
  FaceManager & faceManager = mesh.getFaceManager();

  /// get the array of indicators: 1 if the face is on the boundary; 0 otherwise
  arrayView1d< integer > const & facesDomainBoundaryIndicator = faceManager.getDomainBoundaryIndicator();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

  /// Get table containing all the face normals
  arrayView2d< real64 const > const faceNormal  = faceManager.faceNormal();
  ArrayOfArraysView< localIndex const > const facesToNodes = faceManager.nodeList().toViewConst();

  arrayView1d< real64 > const mass = nodeManager.getExtrinsicData< extrinsicMeshData::MassVector >();

  /// damping matrix to be computed for each dof in the boundary of the mesh
  arrayView1d< real64 > const damping = nodeManager.getExtrinsicData< extrinsicMeshData::DampingVector >();

  damping.setValues< serialPolicy >( 0.0 );

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

      arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();

      arrayView1d< real64 > const rho = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumDensity>();
      arrayView1d< real64 > const vp = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumVelocityVp >();

      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::dispatch3D( fe,
                                 [&]
                                   ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
        constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
        localIndex const numFacesPerElem = elementSubRegion.numFacesPerElement();
        localIndex const numNodesPerFace = 4;

        real64 N[numNodesPerElem];
        real64 gradN[ numNodesPerElem ][ 3 ];

        /// Loop over elements
        for( localIndex k=0; k < elemsToNodes.size( 0 ); ++k )
        {
          real64 xLocal[numNodesPerElem][3];
          for( localIndex a=0; a< numNodesPerElem; ++a )
          {
            for( localIndex i=0; i<3; ++i )
            {
              xLocal[a][i] = X( elemsToNodes( k, a ), i );
            }
          }

          for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
          {
            FE_TYPE::calcN( q, N );
            real64 const detJ = finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, gradN );

            for( localIndex a=0; a< numNodesPerElem; ++a )
            {
              mass[elemsToNodes[k][a]] += rho[k] * detJ * N[a];
            }
          }

          
          real64 const alpha = 1.0/ vp[k];

          for( localIndex kfe=0; kfe< numFacesPerElem; ++kfe )
          {
            localIndex const numFaceGl = elemsToFaces[k][kfe];

            /// Face on the domain boundary and not on free surface
            if( facesDomainBoundaryIndicator[numFaceGl]==1 && freeSurfaceFaceIndicator[numFaceGl]!=1 )
            {
              for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
              {
                FE_TYPE::calcN( q, N );
                real64 const detJ = finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, gradN );

                real64 invJ[3][3]={{0}};
                FE_TYPE::invJacobianTransformation( q, xLocal, invJ );

                for( localIndex a=0; a < numNodesPerFace; ++a )
                {
                  /// compute ds=||detJ*invJ*normalFace_{kfe}||
                  real64 tmp[3]={0};
                  real64 ds = 0.0;
                  for( localIndex i=0; i<3; ++i )
                  {
                    for( localIndex j = 0; j < 3; ++j )
                    {
                      tmp[i] += invJ[j][i]*faceNormal[numFaceGl][j];
                    }
                    ds +=tmp[i]*tmp[i];
                  }
                  ds = std::sqrt( ds );

                  localIndex numNodeGl = facesToNodes[numFaceGl][a];
                  damping[numNodeGl] += alpha*detJ*ds*N[a];
                }
              }
            }
          } // end loop over element
        }
      } );
    } );
  } );

}


void ElasticWaveEquationSEM::applyFreeSurfaceBC( real64 const time, DomainPartition & domain )
{
  FieldSpecificationManager & fsManager = getGlobalState().getFieldSpecificationManager();
  FunctionManager const & functionManager = getGlobalState().getFunctionManager();

  FaceManager & faceManager = domain.getMeshBody( 0 ).getMeshLevel( 0 ).getFaceManager();
  NodeManager & nodeManager = domain.getMeshBody( 0 ).getMeshLevel( 0 ).getNodeManager();

  arrayView2d< real64 > const u_nm1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacement_nm1 >();
  arrayView2d< real64 > const u_n = nodeManager.getExtrinsicData< extrinsicMeshData::Displacement_n >();
  arrayView2d< real64 > const u_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacement_np1 >();

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

      freeSurfaceFaceIndicator.setValues< serialPolicy >( 0 );
      freeSurfaceNodeIndicator.setValues< serialPolicy >( 0 );

      for( localIndex i = 0; i < targetSet.size(); ++i )
      {
        localIndex const kf = targetSet[ i ];
        freeSurfaceFaceIndicator[kf] = 1;

        localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );
        for( localIndex a=0; a < numNodes; ++a )
        {
          localIndex const dof = faceToNodeMap( kf, a );
          freeSurfaceNodeIndicator[dof] = 1;

          u_np1[0][dof] = value;
          u_np1[1][dof] = value;
          u_np1[2][dof] = value;
          u_n[0][dof]   = value;
          u_n[1][dof]   = value;
          u_n[2][dof]   = value;
          u_nm1[0][dof] = value;
          u_nm1[1][dof] = value;
          u_nm1[2][dof] = value;
        }
      }
    }
    else
    {
      GEOSX_ERROR( "This option is not supported yet" );
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


real64 ElasticWaveEquationSEM::evaluateRicker( real64 const & t0, real64 const & f0 )
{
  real64 o_tpeak = 1.0/f0;
  real64 pulse = 0.0;
  if((t0 <= -0.9*o_tpeak) || (t0 >= 2.9*o_tpeak))
    return pulse;

  real64 pi = 2.0*acos( 0 );
  real64 tmp = f0*t0-1.0;
  real64 f0tm1_2 = 2*(tmp*pi)*(tmp*pi);
  real64 gaussian_term = exp( -f0tm1_2 );
  pulse = -(t0-1)*gaussian_term;

  return pulse;
}


real64 ElasticWaveEquationSEM::evaluateRickerOrder1( real64 const & t, real64 const & f0 )
{
  real64 o_tpeak = 1.0/f0;
  real64 pulse = 0.0;

  if((t <= -0.9*o_tpeak) || (t >= 2.9*o_tpeak))
    return pulse;

  real64 pi = 2.0*acos( 0 );
  real64 lam = (f0*pi)*(f0*pi);
  pulse = -2.0*lam*(t-o_tpeak)*exp( -lam*(t-o_tpeak)*(t-o_tpeak));

  return pulse;
}


real64 ElasticWaveEquationSEM::evaluateRickerOrder2( real64 const & t, real64 const & f0 )
{
  real64 o_tpeak = 1.0/f0;
  real64 pulse = 0.0;
  if((t <= -0.9*o_tpeak) || (t >= 2.9*o_tpeak))
    return pulse;

  real64 pi = 2.0*acos( 0 );
  real64 lam = (f0*pi)*(f0*pi);
  pulse = 2.0*lam*(2.0*lam*(t-o_tpeak)*(t-o_tpeak)-1.0)*exp( -lam*(t-o_tpeak)*(t-o_tpeak));

  return pulse;
}



real64 ElasticWaveEquationSEM::explicitStep( real64 const & time_n,
                                              real64 const & dt,
                                              integer const cycleNumber,
                                              DomainPartition & domain )
{

  GEOSX_UNUSED_VAR( time_n, dt, cycleNumber );

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView1d< real64 const > const mass = nodeManager.getExtrinsicData< extrinsicMeshData::MassVector >();
  arrayView1d< real64 const > const damping = nodeManager.getExtrinsicData< extrinsicMeshData::DampingVector >();

  arrayView2d< real64 > const u_nm1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacement_nm1 >();
  arrayView2d< real64 > const u_n = nodeManager.getExtrinsicData< extrinsicMeshData::Displacement_n >();
  arrayView2d< real64 > const u_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Displacement_np1 >();

  /// get array of indicators: 1 if node on free surface; 0 otherwise
  arrayView1d< localIndex const > const freeSurfaceNodeIndicator = nodeManager.getExtrinsicData< extrinsicMeshData::FreeSurfaceNodeIndicator >();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

  arrayView2d< real64 > const stiffnessVector = nodeManager.getExtrinsicData< extrinsicMeshData::StiffnessVector >();

  arrayView2d< real64 > const rhs = nodeManager.getExtrinsicData< extrinsicMeshData::ForcingRHS >();

  forTargetRegionsComplete( mesh, [&]( localIndex const,
                                       localIndex const,
                                       ElementRegionBase & elemRegion )
  {
    elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const,
                                                                       CellElementSubRegion & elementSubRegion )
    {
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();

      arrayView1d< real64 > const vp = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumVelocityVp >();
      arrayView1d< real64 > const vs = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumVelocityVs >();
      arrayView1d< real64 > const rho = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumDensity >();


      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::dispatch3D( fe,
                                 [&]
                                   ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
        constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

        real64 N[numNodesPerElem];
        real64 gradN[ numNodesPerElem ][ 3 ];

        //Declaration of the first derivatives in space of the displacement V2: no arrayView (acollade)
        real64 dux_dx[numNodesPerElem] = {{0.0}};
        real64 duy_dx[numNodesPerElem] = {{0.0}};
        real64 duz_dx[numNodesPerElem] = {{0.0}};
        real64 dux_dy[numNodesPerElem] = {{0.0}};
        real64 duy_dy[numNodesPerElem] = {{0.0}};
        real64 duz_dy[numNodesPerElem] = {{0.0}};
        real64 dux_dz[numNodesPerElem] = {{0.0}};
        real64 duy_dz[numNodesPerElem] = {{0.0}};
        real64 duz_dz[numNodesPerElem] = {{0.0}};

        //Declaration of the nine components of matrix sigma V2: no arrayview
        real64 sigmaxx[numNodesPerElem] = {{0.0}};
        real64 sigmaxy[numNodesPerElem] = {{0.0}};
        real64 sigmaxz[numNodesPerElem] = {{0.0}};
        real64 sigmayx[numNodesPerElem] = {{0.0}};
        real64 sigmayy[numNodesPerElem] = {{0.0}};
        real64 sigmayz[numNodesPerElem] = {{0.0}};
        real64 sigmazx[numNodesPerElem] = {{0.0}};
        real64 sigmazy[numNodesPerElem] = {{0.0}};
        real64 sigmazz[numNodesPerElem] = {{0.0}};

        //Declaration of the Lame coefficients
        real64 lambda[numNodesPerElem] = {{0.0}};
        real64 mu[numNodesPerElem] = {{0.0}};

        // Declaration of the stiffness matrix 'line'
        real64 Rh_ij[3] = {{0.0}};

        for( localIndex k=0; k<elemsToNodes.size( 0 ); ++k )
        {
          real64 xLocal[numNodesPerElem][3] = {{0.0}};

          for( localIndex a=0; a<numNodesPerElem; ++a )
          {
            for( localIndex i=0; i<3; ++i )
            {
              xLocal[a][i] = X( elemsToNodes( k, a ), i );
            }
          }

          for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
          {
            FE_TYPE::calcN( q, N );

            real64 const detJ = finiteElement.template getGradN< FE_TYPE >( k, q, xLocal, gradN );
    

            for( localIndex i=0; i<numNodesPerElem; ++i )
            {
              for(localIndex j=0; j<numNodesPerElem; ++j)
              {

                // Computation of all derivatives of u on the reference element
                dux_dx[i] += u_n[0][j] * gradN[j][0];
                duy_dx[i] += u_n[1][j] * gradN[j][0];
                duz_dx[i] += u_n[2][j] * gradN[j][0];
                dux_dy[i] += u_n[0][j] * gradN[j][1];
                duy_dy[i] += u_n[1][j] * gradN[j][1];
                duz_dy[i] += u_n[2][j] * gradN[j][1];
                dux_dz[i] += u_n[0][j] * gradN[j][2];
                duy_dz[i] += u_n[1][j] * gradN[j][2];
                duz_dz[i] += u_n[2][j] * gradN[j][2];
               
              }

              // Computation of the Lamé coefficients lambda and mu
              mu[i] = rho[i] * vs[i] * vs[i];
              lambda[i] = rho[i] * vp[i] * vp[i] - 2*mu[i];

              // Computation of the stress tensor sigma
              sigmaxx[i] = ((lambda[i] + 2*mu[i]) * dux_dx[i] + lambda[i] * (duy_dy[i] + duz_dz[i]));
              sigmayy[i] = ((lambda[i] + 2*mu[i]) * duy_dy[i] + lambda[i] * (dux_dx[i] + duz_dz[i]));
              sigmazz[i] = ((lambda[i] + 2*mu[i]) * duz_dz[i] + lambda[i] * (duy_dy[i] + dux_dx[i]));
              sigmaxy[i] = (mu[i] * (dux_dy[i] + duy_dx[i]));
              sigmaxz[i] = (mu[i] * (dux_dz[i] + duz_dx[i]));
              sigmayz[i] = (mu[i] * (duz_dy[i] + duy_dz[i]));
              sigmayx[i] = sigmaxy[i];
              sigmazy[i] = sigmayz[i];
              sigmazx[i] = sigmaxz[i];

           
                 
            }
           
            // Computation of the stiffness matrix coefficients
            for ( localIndex i=0; i<numNodesPerElem; i++)
            {
              for ( localIndex j=0; j<numNodesPerElem; j++)
              {
                for ( localIndex a=0; a<3; a++)
                {
                  Rh_ij[a] = detJ * gradN[j][a] * N[i];
                }
               
                // Computation of the stiffness vector coefficients
                stiffnessVector[0][elemsToNodes[k][i]] += Rh_ij[0] * sigmaxx[j] + Rh_ij[1] * sigmaxy[j] + Rh_ij[2] * sigmaxz[i]; 
                stiffnessVector[1][elemsToNodes[k][i]] += Rh_ij[0] * sigmayx[j] + Rh_ij[1] * sigmayy[j] + Rh_ij[2] * sigmayz[i]; 
                stiffnessVector[2][elemsToNodes[k][i]] += Rh_ij[0] * sigmazx[j] + Rh_ij[1] * sigmazy[j] + Rh_ij[2] * sigmazz[i]; 

             }
             
            }
            

          }

        }
      
      } );
    } );
  } );

  addSourceToRightHandSide( time_n, rhs );

  /// Calculate your time integrators
  real64 dt2 = dt*dt;
  for( localIndex a=0; a<nodeManager.size(); ++a )
  {
    if( freeSurfaceNodeIndicator[a]!=1 )
    {
      for ( localIndex k=0; k<3; k++)
      {
        u_np1[k][a] = (1.0/(mass[a] + 0.5*dt*damping[a]))*(2*mass[a] * u_n[k][a] - dt2*stiffnessVector[k][a] - (mass[a] - 0.5*dt*damping[a]) * u_nm1[k][a] + dt2*rhs[k][a] );
      }
      
    }
  }

  /// Synchronize pressure fields
  std::map< string, string_array > fieldNames;
  fieldNames["node"].emplace_back( "displacement_np1" );

  CommunicationTools syncFields;
  syncFields.synchronizeFields( fieldNames,
                                domain.getMeshBody( 0 ).getMeshLevel( 0 ),
                                domain.getNeighbors(),
                                true );

  for( localIndex a=0; a<nodeManager.size(); ++a )
  {
    for ( localIndex k=0; k<3; k++)
    {
      u_nm1[k][a] = u_n[k][a];
      u_n[k][a] = u_np1[k][a];

      stiffnessVector[k][a] = 0.0;
      rhs[k][a] = 0.0;
    }
    
  }

  computeSismoTrace( cycleNumber, u_np1 );


  return dt;
  
}

REGISTER_CATALOG_ENTRY( SolverBase, ElasticWaveEquationSEM, string const &, dataRepository::Group * const )

} /* namespace geosx */
