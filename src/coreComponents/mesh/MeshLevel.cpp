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
 * @file MeshLevel.cpp
 */

#include "MeshLevel.hpp"

#include "EdgeManager.hpp"
#include "ElementRegionManager.hpp"
#include "NodeManager.hpp"
#include "FaceManager.hpp"

namespace geosx
{
using namespace dataRepository;

MeshLevel::MeshLevel( string const & name,
                      Group * const parent ):
  Group( name, parent ),
  m_nodeManager( groupStructKeys::nodeManagerString, this ),
  m_edgeManager( groupStructKeys::edgeManagerString, this ),
  m_faceManager( groupStructKeys::faceManagerString, this ),
  m_elementManager( groupStructKeys::elemManagerString, this ),
  m_embSurfNodeManager( groupStructKeys::embSurfNodeManagerString, this ),
  m_embSurfEdgeManager( groupStructKeys::embSurfEdgeManagerString, this )

{

  registerGroup( groupStructKeys::nodeManagerString, &m_nodeManager );

  registerGroup( groupStructKeys::edgeManagerString, &m_edgeManager );


  registerGroup< FaceManager >( groupStructKeys::faceManagerString, &m_faceManager );
  m_faceManager.nodeList().setRelatedObject( m_nodeManager );


  registerGroup< ElementRegionManager >( groupStructKeys::elemManagerString, &m_elementManager );

  registerGroup< EdgeManager >( groupStructKeys::embSurfEdgeManagerString, &m_embSurfEdgeManager );

  registerGroup< EmbeddedSurfaceNodeManager >( groupStructKeys::embSurfNodeManagerString, &m_embSurfNodeManager );

  registerWrapper< integer >( viewKeys.meshLevel );
}

MeshLevel::MeshLevel( string const & name,
                      Group * const parent,
                      MeshLevel const & source,
                      int const order ):
  MeshLevel( name, parent )
{

  localIndex const numBasisSupportPoints = order+1;


  // Edges
  localIndex const numNodesPerEdge = numBasisSupportPoints;
  localIndex const numNonVertexNodesPerEdge = numNodesPerEdge - 2;

  m_edgeManager.resize( source.m_edgeManager.size() );
  EdgeManager::NodeMapType const & edgesToNodesSource = source.m_edgeManager.nodeList();
  EdgeManager::NodeMapType const & edgesToNodes = m_edgeManager.nodeList();

  localIndex const numInternalEdgeNodes = m_edgeManager.size() * numNonVertexNodesPerEdge;

  m_faceManager.resize( source.m_faceManager.size() );
  m_faceManager.edgeList() = source.m_faceManager.edgeList();

  // Faces
  ArrayOfArraysView< localIndex const > const facesToNodesMapSource = m_faceManager.nodeList().toViewConst();
  ArrayOfArrays< localIndex > & faceToNodeMapNew = m_faceManager.nodeList();
  ArrayOfArraysView< localIndex const > const & facesToEdges = m_faceManager.edgeList().toViewConst();
  localIndex const estimatedNumNodesPerFace = pow( order+1, 2 );
  faceToNodeMapNew.resize( faceToNodeMapNew.size(), estimatedNumNodesPerFace );

  // add the number of non-edge face nodes
  localIndex numInternalFaceNodes = 0;
  for( localIndex kf=0; kf<m_faceManager.size(); ++kf )
  {
    localIndex const numEdgesPerFace = facesToEdges.sizeOfArray( kf );
    localIndex const numVertexNodesPerFace = facesToNodesMapSource.sizeOfArray( kf );
    localIndex const numEdgeNodesPerFace = numVertexNodesPerFace + numEdgesPerFace * numNonVertexNodesPerEdge;

    if( numEdgesPerFace==4 )
    {
      localIndex const numNonEdgeNodesPerFace = pow( order-1, 2 );
      numInternalFaceNodes += numNonEdgeNodesPerFace;

      localIndex const numNodesPerFace = pow( order+1, 2 );
      faceToNodeMapNew.resizeArray( kf, numNodesPerFace );
    }
    else
    {
      GEOSX_ERROR( "need more support for face geometry" );
    }
  }

  // add the number of non-face element nodes
  localIndex numInternalElementNodes = 0;
  source.m_elementManager.forElementRegions< CellElementRegion >( [&]( CellElementRegion const & sourceRegion )
  {
    sourceRegion.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & sourceSubRegion )
    {
      if( sourceSubRegion.getElementType() == ElementType::Hexahedron )
      {
        numInternalElementNodes += sourceSubRegion.size() * pow( order-1, 3 );
      }
    } );
  } );



  localIndex const numNodes = source.m_nodeManager.size()
                              + numInternalEdgeNodes
                              + numInternalFaceNodes
                              + numInternalElementNodes;

  m_nodeManager.resize( numNodes );

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const refPosSource = source.m_nodeManager.referencePosition();
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const refPosNew = m_nodeManager.referencePosition().toView();

  {
    Group & nodeSets = m_nodeManager.sets();
    SortedArray< localIndex > & allNodes  = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "all" ) ).reference();
    allNodes.reserve( m_nodeManager.size() );

    for( localIndex a=0; a<m_nodeManager.size(); ++a )
    {
      allNodes.insert( a );
    }

  }



  ArrayOfArraysView< localIndex const > const & faceToNodeMapSource = source.m_faceManager.nodeList().toViewConst();

  FaceManager::ElemMapType const & faceToElem = source.m_faceManager.toElementRelation();

  arrayView2d< localIndex const > const & faceToElemIndex = faceToElem.m_toElementIndex.toViewConst();


  source.m_elementManager.forElementRegions< CellElementRegion >( [&]( CellElementRegion const & sourceRegion )
  {
    CellElementRegion & region = *(dynamic_cast< CellElementRegion * >( m_elementManager.createChild( sourceRegion.getCatalogName(),
                                                                                                      sourceRegion.getName() ) ) );

    region.addCellBlockNames( sourceRegion.getCellBlockNames() );

    sourceRegion.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & sourceSubRegion )
    {
      localIndex const numNodesPerElem = pow( order+1, 3 );

      CellElementSubRegion & newSubRegion = region.getSubRegions().registerGroup< CellElementSubRegion >( sourceSubRegion.getName() );
      newSubRegion.setElementType( sourceSubRegion.getElementType() );

      newSubRegion.resizePerElementValues( numNodesPerElem,
                                           sourceSubRegion.numEdgesPerElement(),
                                           sourceSubRegion.numFacesPerElement() );

      newSubRegion.resize( sourceSubRegion.size() );

      arrayView2d< localIndex const > const & elemToFaces = sourceSubRegion.faceList();

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodesSource = sourceSubRegion.nodeList().toViewConst();
      array2d< localIndex, cells::NODE_MAP_PERMUTATION > & elemsToNodesNew = newSubRegion.nodeList();

      array2d< localIndex > & elemToFacesNew = newSubRegion.faceList();



      //Copy a new elemToFaces map from the old one
      for( localIndex elem = 0; elem < elemsToNodesNew.size( 0 ); ++elem )
      {
        for( localIndex face = 0; face < 6; ++face )
        {
          elemToFacesNew[elem][face] = elemToFaces[elem][face];
        }
      }



//      elemsToNodesNew.resize( elemsToNodesSource.size(0), numNodesPerElem );

      // Fill a temporary table which knowing the global number of a degree of freedom and a face, gives you the local number of this degree
      // of freedom on the considering face
      array2d< localIndex > localElemToLocalFace( 6, numNodesPerElem );

      //Init arrays
      for( localIndex i = 0; i < 6; ++i )
      {
        for( localIndex j = 0; j < pow( order+1, 3 ); ++j )
        {
          localElemToLocalFace[i][j]=-1;
        }

      }



      for( localIndex i = 0; i < faceToNodeMapSource.size(); ++i )
      {
        for( localIndex j = 0; j < pow( order+1, 2 ); ++j )
        {
          faceToNodeMapNew[i][j] = -1;
        }

      }

      //Face 0
      for( localIndex i = 0; i < order+1; i++ )
      {
        for( localIndex j = 0; j <order+1; j++ )
        {
          localElemToLocalFace[0][i + pow( order+1, 2 )*j] = i + (order+1)*j;
        }

      }

      //Face 1
      for( localIndex i = 0; i < order+1; i++ )
      {
        for( localIndex k = 0; k < order+1; k++ )
        {
          localElemToLocalFace[1][k+(order+1)*i] = k + (order+1)*i;
        }

      }

      //Face 2
      for( localIndex k = 0; k < order+1; k++ )
      {
        for( localIndex j = 0; j < order+1; j++ )
        {
          localElemToLocalFace[2][k*pow( order+1, 2 )+j*(order+1)] = j + (order+1)*k;
        }

      }

      //Face 3
      for( localIndex j = 0; j < order+1; j++ )
      {
        for( localIndex k = 0; k < order+1; k++ )
        {
          localElemToLocalFace[3][order +k*(order+1)+j*pow( order+1, 2 )] = k + (order+1)*j;
        }

      }

      //Face 4
      for( localIndex j = 0; j < order+1; j++ )
      {
        for( localIndex i = 0; i < order+1; i++ )
        {
          localElemToLocalFace[4][order*(order+1)+i+j*pow( order+1, 2 )] = i + (order+1)*j;
        }

      }

      //Face 5
      for( localIndex k = 0; k < order+1; k++ )
      {
        for( localIndex i = 0; i < order+1; i++ )
        {
          localElemToLocalFace[5][order*pow( order+1, 2 )+i+k*(order+1)] = i + (order+1)*k;
        }

      }

      //Initialisation of elemToNodes
      for( localIndex e = 0; e < elemsToNodesNew.size( 0 ); ++e )
      {
        for( localIndex i = 0; i < numNodesPerElem; i++ )
        {
          elemsToNodesNew[e][i] =-1;
        }

      }

      localIndex count=0;

      for( localIndex elem = 0; elem < elemsToNodesNew.size( 0 ); ++elem )
      {
        for( localIndex k = 0; k < order+1; ++k )
        {
          for( localIndex j = 0; j< order+1; ++j )
          {
            for( localIndex i = 0; i <order+1; ++i )
            {

              localIndex face = 0;
              localIndex foundFace = 0;
              while( face<6 && foundFace<1 )
              {
                localIndex m = localElemToLocalFace[face][i +(order+1)*j + pow( order+1, 2 )*k];

                if( m != -1 )
                {
                  if( faceToNodeMapNew[elemToFaces[elem][face]][m] != -1 )
                  {
                    foundFace = 1;
                    for( localIndex l = 0; l < 2; ++l )
                    {
                      localIndex elemNeighbour = faceToElemIndex[elemToFaces[elem][face]][l];
                      if( elemNeighbour != elem && elemNeighbour != -1 )
                      {
                        for( localIndex node = 0; node < pow( order+1, 3 ); ++node )
                        {
                          if( elemsToNodesNew[elemNeighbour][node] == faceToNodeMapNew[elemToFaces[elem][face]][m] )
                          {
                            elemsToNodesNew[elem][i + (order+1)*j + pow( order+1, 2 )*k] = elemsToNodesNew[elemNeighbour][node];
                            break;
                          }
                        }
                        break;
                      }
                    }
                    for( localIndex face2 = 0; face2 < 6; ++face2 )
                    {
                      localIndex m2 = localElemToLocalFace[face2][i +(order+1)*j + pow( order+1, 2 )*k];

                      if( m2 != -1 && face2!=face )
                      {
                        faceToNodeMapNew[elemToFaces[elem][face2]][m2] = faceToNodeMapNew[elemToFaces[elem][face]][m];
                      }
                    }
                  }
                  else
                  {
                    faceToNodeMapNew[elemToFaces[elem][face]][m] = count;
                  }
                }
                face++;
              }
              if( face > 5 && foundFace < 1 )
              {
                elemsToNodesNew[elem][i + (order+1)*j + pow( order+1, 2 )*k] = count;
                count++;
              }
            }
          }
        }
      }

      //Fill a temporary array which contains the Gauss-Lobatto points depending on the order
      array1d< real64 > GaussLobattoPts( 4 );

      if( order==1 )
      {
        GaussLobattoPts[0] = -1.0;
        GaussLobattoPts[1] = 1.0;
      }

      if( order==3 )
      {
        GaussLobattoPts[0] = -1.0;
        GaussLobattoPts[1] = -1./sqrt( 5 );
        GaussLobattoPts[2] = 1./sqrt( 5 );
        GaussLobattoPts[3] = 1.;
      }

      if( order==5 )
      {
        static constexpr real64 sqrt__7_plus_2sqrt7__ = 3.50592393273573196;
        static constexpr real64 sqrt__7_mins_2sqrt7__ = 1.30709501485960033;

        static constexpr real64 sqrt_inv21 = 0.218217890235992381;

        GaussLobattoPts[0] = -1.0;
        GaussLobattoPts[1] = -sqrt_inv21*sqrt__7_plus_2sqrt7__;
        GaussLobattoPts[2] = -sqrt_inv21*sqrt__7_mins_2sqrt7__;
        GaussLobattoPts[3] = sqrt_inv21*sqrt__7_mins_2sqrt7__;
        GaussLobattoPts[4] = sqrt_inv21*sqrt__7_plus_2sqrt7__;
        GaussLobattoPts[5] = 1.0;
      }

      //Three 1D arrays to contains the GL points in the new coordinates knowing the mesh nodes
      array1d< real64 > x( order+1 );
      array1d< real64 > y( order+1 );
      array1d< real64 > z( order+1 );

      for( localIndex e = 0; e < elemsToNodesNew.size( 0 ); e++ )
      {
        //Fill the three 1D array
        for( localIndex i = 0; i < order+1; i++ )
        {
          x[i] = refPosSource[elemsToNodesSource[e][0]][0] + ((refPosSource[elemsToNodesSource[e][1]][0]-refPosSource[elemsToNodesSource[e][0]][0])/2.)*GaussLobattoPts[i]
                 + (refPosSource[elemsToNodesSource[e][1]][0]-refPosSource[elemsToNodesSource[e][0]][0])/2;
          y[i] = refPosSource[elemsToNodesSource[e][0]][1] + ((refPosSource[elemsToNodesSource[e][2]][1]-refPosSource[elemsToNodesSource[e][0]][1])/2.)*GaussLobattoPts[i]
                 + (refPosSource[elemsToNodesSource[e][2]][1]-refPosSource[elemsToNodesSource[e][0]][1])/2;
          z[i] = refPosSource[elemsToNodesSource[e][0]][2] + ((refPosSource[elemsToNodesSource[e][4]][2]-refPosSource[elemsToNodesSource[e][0]][2])/2.)*GaussLobattoPts[i]
                 + (refPosSource[elemsToNodesSource[e][4]][2]-refPosSource[elemsToNodesSource[e][0]][2])/2;
        }


        //Fill new refPos array
        for( localIndex k = 0; k< order+1; k++ )
        {
          for( localIndex j = 0; j < order+1; j++ )
          {
            for( localIndex i = 0; i < order+1; i++ )
            {
              localIndex const nodeIndex = elemsToNodesNew( e, i+j*(order+1)+k*pow( order+1, 2 ) );

              refPosNew( nodeIndex, 0 ) = x[i];
              refPosNew( nodeIndex, 1 ) = y[j];
              refPosNew( nodeIndex, 2 ) = z[k];
            }

          }

        }

      }
    } );
  } );
}



MeshLevel::~MeshLevel()
{}

void MeshLevel::initializePostInitialConditionsPostSubGroups()
{
  m_elementManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {
    subRegion.calculateElementGeometricQuantities( m_nodeManager, m_faceManager );
  } );
}

void MeshLevel::generateAdjacencyLists( arrayView1d< localIndex const > const & seedNodeList,
                                        localIndex_array & nodeAdjacencyList,
                                        localIndex_array & edgeAdjacencyList,
                                        localIndex_array & faceAdjacencyList,
                                        ElementRegionManager::ElementViewAccessor< ReferenceWrapper< localIndex_array > > & elementAdjacencyList,
                                        integer const depth )
{
  NodeManager const & nodeManager = getNodeManager();

  ArrayOfArraysView< localIndex const > const & nodeToElementRegionList = nodeManager.elementRegionList().toViewConst();
  ArrayOfArraysView< localIndex const > const & nodeToElementSubRegionList = nodeManager.elementSubRegionList().toViewConst();
  ArrayOfArraysView< localIndex const > const & nodeToElementList = nodeManager.elementList().toViewConst();

  FaceManager const & faceManager = this->getFaceManager();
  ArrayOfArraysView< localIndex const > const & faceToEdges = faceManager.edgeList().toViewConst();

  ElementRegionManager const & elemManager = this->getElemManager();

  std::set< localIndex > nodeAdjacencySet;
  std::set< localIndex > edgeAdjacencySet;
  std::set< localIndex > faceAdjacencySet;
  std::vector< std::vector< std::set< localIndex > > > elementAdjacencySet( elemManager.numRegions() );

  for( localIndex a=0; a<elemManager.numRegions(); ++a )
  {
    elementAdjacencySet[a].resize( elemManager.getRegion( a ).numSubRegions() );
  }

  nodeAdjacencySet.insert( seedNodeList.begin(), seedNodeList.end() );

  for( integer d=0; d<depth; ++d )
  {
    for( localIndex const nodeIndex : nodeAdjacencySet )
    {
      for( localIndex b=0; b<nodeToElementRegionList.sizeOfArray( nodeIndex ); ++b )
      {
        localIndex const regionIndex = nodeToElementRegionList[nodeIndex][b];
        localIndex const subRegionIndex = nodeToElementSubRegionList[nodeIndex][b];
        localIndex const elementIndex = nodeToElementList[nodeIndex][b];
        elementAdjacencySet[regionIndex][subRegionIndex].insert( elementIndex );
      }
    }

    for( typename dataRepository::indexType kReg=0; kReg<elemManager.numRegions(); ++kReg )
    {
      ElementRegionBase const & elemRegion = elemManager.getRegion( kReg );

      elemRegion.forElementSubRegionsIndex< CellElementSubRegion,
                                            WellElementSubRegion >( [&]( localIndex const kSubReg,
                                                                         auto const & subRegion )
      {
        arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = subRegion.nodeList();
        arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList();
        for( auto const elementIndex : elementAdjacencySet[kReg][kSubReg] )
        {
          for( localIndex a=0; a<elemsToNodes.size( 1 ); ++a )
          {
            nodeAdjacencySet.insert( elemsToNodes[elementIndex][a] );
          }

          for( localIndex a=0; a<elemsToFaces.size( 1 ); ++a )
          {
            faceAdjacencySet.insert( elemsToFaces[elementIndex][a] );

            localIndex const faceIndex = elemsToFaces[elementIndex][a];
            localIndex const numEdges = faceToEdges.sizeOfArray( faceIndex );
            for( localIndex b=0; b<numEdges; ++b )
            {
              edgeAdjacencySet.insert( faceToEdges( faceIndex, b ));
            }

          }

        }
      } );
    }
  }

  nodeAdjacencyList.resize( LvArray::integerConversion< localIndex >( nodeAdjacencySet.size()));
  std::copy( nodeAdjacencySet.begin(), nodeAdjacencySet.end(), nodeAdjacencyList.begin() );

  edgeAdjacencyList.resize( LvArray::integerConversion< localIndex >( edgeAdjacencySet.size()));
  std::copy( edgeAdjacencySet.begin(), edgeAdjacencySet.end(), edgeAdjacencyList.begin() );

  faceAdjacencyList.resize( LvArray::integerConversion< localIndex >( faceAdjacencySet.size()));
  std::copy( faceAdjacencySet.begin(), faceAdjacencySet.end(), faceAdjacencyList.begin() );

  for( localIndex kReg=0; kReg<elemManager.numRegions(); ++kReg )
  {
    ElementRegionBase const & elemRegion = elemManager.getRegion( kReg );

    for( localIndex kSubReg = 0; kSubReg < elemRegion.numSubRegions(); ++kSubReg )
    {
      elementAdjacencyList[kReg][kSubReg].get().clear();
      elementAdjacencyList[kReg][kSubReg].get().resize( LvArray::integerConversion< localIndex >( elementAdjacencySet[kReg][kSubReg].size()) );
      std::copy( elementAdjacencySet[kReg][kSubReg].begin(),
                 elementAdjacencySet[kReg][kSubReg].end(),
                 elementAdjacencyList[kReg][kSubReg].get().begin() );

    }
  }

}

void MeshLevel::generateSets()
{
  GEOSX_MARK_FUNCTION;

  NodeManager const & nodeManager = m_nodeManager;

  dataRepository::Group const & nodeSets = nodeManager.sets();

  map< string, array1d< bool > > nodeInSet; // map to contain indicator of whether a node is in a set.
  string_array setNames; // just a holder for the names of the sets

  // loop over all wrappers and fill the nodeIndSet arrays for each set
  for( auto & wrapper : nodeSets.wrappers() )
  {
    string const & name = wrapper.second->getName();
    nodeInSet[name].resize( nodeManager.size() );
    nodeInSet[name].setValues< serialPolicy >( false );

    if( nodeSets.hasWrapper( name ) )
    {
      setNames.emplace_back( name );
      SortedArrayView< localIndex const > const & set = nodeSets.getReference< SortedArray< localIndex > >( name );
      for( localIndex const a : set )
      {
        nodeInSet[name][a] = true;
      }
    }
  }


  ElementRegionManager & elementRegionManager = m_elementManager;
  elementRegionManager.forElementSubRegions( [&]( auto & subRegion )
  {
    dataRepository::Group & elementSets = subRegion.sets();

    auto const & elemToNodeMap = subRegion.nodeList();

    for( string const & setName : setNames )
    {
      arrayView1d< bool const > const nodeInCurSet = nodeInSet[setName];

      SortedArray< localIndex > & targetSet = elementSets.registerWrapper< SortedArray< localIndex > >( setName ).reference();
      for( localIndex k = 0; k < subRegion.size(); ++k )
      {
        localIndex const numNodes = subRegion.numNodesPerElement( k );

        localIndex elementInSet = true;
        for( localIndex i = 0; i < numNodes; ++i )
        {
          if( !nodeInCurSet( elemToNodeMap[ k ][ i ] ) )
          {
            elementInSet = false;
            break;
          }
        }

        if( elementInSet )
        {
          targetSet.insert( k );
        }
      }
    }
  } );
}


} /* namespace geosx */
