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


  localIndex numNodes = source.m_nodeManager.size()+source.m_edgeManager.size()*(order-1)+pow(order,2)*source.m_faceManager.size()+pow(order-1,3)*source.m_elementManager.size();
  // find out how many node there must be on this rank
  m_nodeManager.resize(numNodes);

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const refPosSource = source.m_nodeManager.referencePosition();
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const refPosNew = m_nodeManager.referencePosition().toView();

  // {
  //   Group & nodeSets = m_nodeManager.sets();
  //   SortedArray< localIndex > & allNodes  = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "all" ) ).reference();
  //   allNodes.reserve( m_nodeManager.size() );

  //   for( localIndex a=0; a<m_nodeManager.size(); ++a )
  //   {
  //     allNodes.insert( a );
  //   }

  // }

  // forAll<parallelDevicePolicy<>>( m_nodeManager.size(),
  //                                   [=]( localIndex const a )
  // {
  //   for( localIndex i=0; i<3; ++i )
  //   {
  //     refPosNew(a,i) = refPosSource(a,i); // this needs to be another loop with a linear combination of values.
  //   }
  // });
  
  //ArrayOfArraysView< localIndex const > const faceToNodeMapSource = m_faceManager.nodeList().toViewConst();
  ArrayOfArrays< localIndex > faceToNodeMapNew = m_faceManager.nodeList();
  
  localIndex const numNodesPerFace = pow(order+1,2);
  faceToNodeMapNew.resize(faceToNodeMapNew.size(),numNodesPerFace);
  

  FaceManager::ElemMapType const & faceToElem = m_faceManager.toElementRelation();
  arrayView2d< localIndex const > const & faceToElemIndex = faceToElem.m_toElementIndex.toViewConst();

  
  source.m_elementManager.forElementRegions<CellElementRegion>([&]( CellElementRegion const & sourceRegion )
  {
    CellElementRegion & region = *(dynamic_cast<CellElementRegion *>( m_elementManager.createChild( sourceRegion.getCatalogName(),
                                                                                                        sourceRegion.getName() ) ) );

    region.addCellBlockNames( sourceRegion.getCellBlockNames() );
    // std::exit(2);
    sourceRegion.forElementSubRegions<CellElementSubRegion>( [&]( CellElementSubRegion const & sourceSubRegion )
    {
//      std::exit(2);
      std::cout << "hello" << std::endl;
      CellElementSubRegion & newSubRegion = region.getSubRegions().registerGroup< CellElementSubRegion >( sourceSubRegion.getName() );
      newSubRegion.setElementType( sourceSubRegion.getElementType() );

      newSubRegion.resize( sourceSubRegion.size() );
      
      
      // Group & sets = newSubRegion.sets();
      // SortedArray< localIndex > & allElems  = sets.registerWrapper< SortedArray< localIndex > >( string( "all" ) ).reference();


      arrayView2d< localIndex const > const & elemToFaces = sourceSubRegion.faceList();

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodesSource = sourceSubRegion.nodeList().toViewConst();
      array2d< localIndex, cells::NODE_MAP_PERMUTATION > elemsToNodesNew = newSubRegion.nodeList();


      localIndex const numNodesPerElem = pow(order+1,3);// change to pow( 2+( order>1 ? order-1 : 0 ), 3 );
      elemsToNodesNew.resize( elemsToNodesSource.size(0), numNodesPerElem );
    
      // Fill a temporary table which knowing the global number of a degree of freedom and a face, gives you the local number of this degree
      // of freedom on the considering face
      array2d < localIndex > localElemToLocalFace(6,numNodesPerElem);


      for (localIndex i = 0; i < order+1; ++i)//x
      {
        for (localIndex j = 0; j < order+1; ++j)//y
        {
          for (localIndex k = 0; k < order+1; ++k)//z
          {
            localElemToLocalFace[0][(order+1)*order - (order+1)*j + pow(order+1,2)*k] = j + (order+1)*k;

            localElemToLocalFace[1][(order+1)*order + (order+1)*j + pow(order+1,2)*k] = j + (order+1)*k;

            localElemToLocalFace[2][i + pow(order+1,2)*j] = i + (order+1)*j;

            localElemToLocalFace[3][(pow(order+1,2)-1) - i + (order+1)*j] = i + (order+1)*j;

            localElemToLocalFace[4][i + (order+1)*(order-k)] = i + (order+1)*k;

            localElemToLocalFace[5][pow(order+1,2)*order + i + (order+1)*k] = i + (order+1)*k;
          }
        }
      }

      localIndex count=0;

      for (localIndex e = 0; e < elemsToNodesNew.size(0); e++)
      {
        for (localIndex k = 0; k < order+1; ++k)
        {
          for (localIndex j = 0; j< order+1; ++j)
          {
            for (localIndex i = 0; i < order+1; ++i)
            {
              for (localIndex face = 0; face < 6; ++face)
              {
                localIndex m = localElemToLocalFace[face][i +(order+1)*j + pow(order+1,2)*k];
                if (m != -1)
                {
                  if ( faceToNodeMapNew[elemToFaces[e][face]][m] != -1)
                  {
                    for (localIndex l = 0; i < 2; ++l)
                    {
                      localIndex neighE = faceToElemIndex[elemToFaces[e][face]][l];
                      if (neighE != e && neighE != -1)
                      {
                        for (localIndex n = 0; n < pow(order+1,3); ++n)
                        {
                          if ( elemsToNodesNew[neighE][n] == faceToNodeMapNew[elemToFaces[e][face]][i + (order+1)*j + pow(order+1,2)*k])
                          {
                            elemsToNodesNew[e][i + (order+1)*j + pow(order+1,2)*k] = elemsToNodesNew[neighE][n];
                          }
                          break;
                        }
                        break;
                      }
                    }
                  }
                  else
                  {
                    elemsToNodesNew[e][i + order*j + pow(order,2)*k] = count;

                    faceToNodeMapNew[elemToFaces[e][face]][m] = count;

                    count++;
                  }
                }
                else
                {
                  elemsToNodesNew[e][i + (order+1)*j + pow(order+1,2)*k]=count;

                  count++;
                }
              }
            }
          }
        }
      }

  
    // Group & nodeSets = m_nodeManager.sets();
    // SortedArray< localIndex > & allNodes  = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "all" ) ).reference();
    // allNodes.reserve( m_nodeManager.size() );
    // std::exit(2);
    // for( localIndex a=0; a<m_nodeManager.size(); ++a )
    // {
    //   allNodes.insert( a );
    // }
    

      //FIll a temporary array which contains the Gauss-Lobatto points depending on the order
      array1d< real64 > GaussLobattoPts(order+1);

      if(order==1)
      {
        GaussLobattoPts[0] = -1.0;
        GaussLobattoPts[1] = 1.0;
      }
  
      if(order==3)
      {
        GaussLobattoPts[0] = -1.0;
        GaussLobattoPts[1] = -1./sqrt(5);
        GaussLobattoPts[2] = 1./sqrt(5);
        GaussLobattoPts[3] = 1.;
      }

      if (order==5)
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

      std::cout << elemsToNodesNew.size(0) << std::endl;
      std::exit(2);
      //Fill Position of the nodes: interpolation of the mesh nodes
      for (localIndex e = 0; e < elemsToNodesNew.size(0); e++)
      {
        //Get the space step of the mesh

        std::cout << refPosSource[elemsToNodesNew[e][order+1]][0] << std::endl;
        std::cout << refPosSource[elemsToNodesNew[e][0]][0] << std::endl;
        std::exit(2);
        real64 hx = refPosSource[elemsToNodesNew[e][order+1]][0] - refPosSource[elemsToNodesNew[e][0]][0];
        real64 hy = refPosSource[elemsToNodesNew[e][order*(order+1)]][1] - refPosSource[elemsToNodesNew[e][0]][1];
        real64 hz = refPosSource[elemsToNodesNew[e][pow(order+1,2)]][2] - refPosSource[elemsToNodesNew[e][0]][2];
        //Fill the array
        for (localIndex k = 0; k < order+1; k++)
        {
          for (localIndex j = 0; j < order+1; j++)
          {
            for (localIndex i = 0; i < order+1; i++)
            {
              refPosNew[elemsToNodesNew[e][i+j*(order+1)+k*pow(order+1,2)]][0] = refPosSource[elemsToNodesNew[e][0]][0] + i*hx*GaussLobattoPts[i];
              refPosNew[elemsToNodesNew[e][i+j*(order+1)+k*pow(order+1,2)]][1] = refPosSource[elemsToNodesNew[e][0]][1] + i*hy*GaussLobattoPts[j];
              refPosNew[elemsToNodesNew[e][i+j*(order+1)+k*pow(order+1,2)]][0] = refPosSource[elemsToNodesNew[e][0]][2] + i*hz*GaussLobattoPts[k];
            }
          }
        }
      }

    });
  });
  // localIndex numNodes = source.m_nodeManager.size()+source.m_edgeManager.size()*(order-1)+pow(order,2)*source.m_faceManager.size()+pow(order-1,3)*source.m_elementManager.size();
  // // find out how many node there must be on this rank
  // m_nodeManager.resize(numNodes);
  // array2d< localIndex, cells::NODE_MAP_PERMUTATION > elemsToNodesNew = newSubRegion.nodeList();
  // arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const refPosSource = source.m_nodeManager.referencePosition();
  // arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const refPosNew = m_nodeManager.referencePosition().toView();

  // {
  //   Group & nodeSets = m_nodeManager.sets();
  //   SortedArray< localIndex > & allNodes  = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "all" ) ).reference();
  //   allNodes.reserve( m_nodeManager.size() );

  //   for( localIndex a=0; a<m_nodeManager.size(); ++a )
  //   {
  //     allNodes.insert( a );
  //   }

  // }
  // //FIll a temporary array which contains the Gauss-Lobatto points depending on the order
  //     array1d< real64 > GaussLobattoPts(order+1);

  //     if(order==3)
  //     {
  //       GaussLobattoPts[0] = -1.0;
  //       GaussLobattoPts[1] = -1./sqrt(5);
  //       GaussLobattoPts[2] = 1./sqrt(5);
  //       GaussLobattoPts[3] = 1.;
  //     }

  //     if (order==5)
  //     {
  //       static constexpr real64 sqrt__7_plus_2sqrt7__ = 3.50592393273573196;
  //       static constexpr real64 sqrt__7_mins_2sqrt7__ = 1.30709501485960033;

  //       static constexpr real64 sqrt_inv21 = 0.218217890235992381;

  //       GaussLobattoPts[0] = -1.0;
  //       GaussLobattoPts[1] = -sqrt_inv21*sqrt__7_plus_2sqrt7__;
  //       GaussLobattoPts[2] = -sqrt_inv21*sqrt__7_mins_2sqrt7__;
  //       GaussLobattoPts[3] = sqrt_inv21*sqrt__7_mins_2sqrt7__;
  //       GaussLobattoPts[4] = sqrt_inv21*sqrt__7_plus_2sqrt7__;
  //       GaussLobattoPts[5] = 1.0;
  //     }


  //     //Fill Position of the nodes: interpolation of the mesh nodes
  //     for (localIndex e = 0; e < elemsToNodesNew.size(0); e++)
  //     {
  //       //Get the space step of the mesh
  //       real64 hx = refPosSource[elemsToNodesNew[e][order+1]][0] - refPosSource[elemsToNodesNew[e][0]][0];
  //       real64 hy = refPosSource[elemsToNodesNew[e][order*(order+1)]][1] - refPosSource[elemsToNodesNew[e][0]][1];
  //       real64 hz = refPosSource[elemsToNodesNew[e][pow(order+1,2)]][2] - refPosSource[elemsToNodesNew[e][0]][2];
  //       //Fill the array
  //       for (localIndex k = 0; k < order+1; k++)
  //       {
  //         for (localIndex j = 0; j < order+1; j++)
  //         {
  //           for (localIndex i = 0; i < order+1; i++)
  //           {
  //             refPosNew[elemsToNodesNew[e][i+j*(order+1)+k*pow(order+1,2)]][0] = refPosSource[elemsToNodesNew[e][0]][0] + i*hx*GaussLobattoPts[i];
  //             refPosNew[elemsToNodesNew[e][i+j*(order+1)+k*pow(order+1,2)]][1] = refPosSource[elemsToNodesNew[e][0]][1] + i*hy*GaussLobattoPts[j];
  //             refPosNew[elemsToNodesNew[e][i+j*(order+1)+k*pow(order+1,2)]][0] = refPosSource[elemsToNodesNew[e][0]][2] + i*hz*GaussLobattoPts[k];
  //           }
  //         }
  //       }
  //     }
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
  NodeManager & nodeManager = getNodeManager();

  ArrayOfArraysView< localIndex const > const & nodeToElementRegionList = nodeManager.elementRegionList().toViewConst();

  ArrayOfArraysView< localIndex const > const & nodeToElementSubRegionList = nodeManager.elementSubRegionList().toViewConst();

  ArrayOfArraysView< localIndex const > const & nodeToElementList = nodeManager.elementList().toViewConst();


  FaceManager & faceManager = this->getFaceManager();
  ArrayOfArraysView< localIndex const > const & faceToEdges = faceManager.edgeList().toViewConst();

  ElementRegionManager & elemManager = this->getElemManager();

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

            localIndex const faceID = elemsToFaces[elementIndex][a];
            localIndex const numEdges = faceToEdges.sizeOfArray( faceID );
            for( localIndex b=0; b<numEdges; ++b )
            {
              edgeAdjacencySet.insert( faceToEdges( faceID, b ));
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
