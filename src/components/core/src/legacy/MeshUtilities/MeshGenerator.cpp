// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
/*
 * MeshGenerator.cpp
 *
 *  Created on: Nov 19, 2012
 *      Author: settgast1
 */

#include "MeshGenerator.h"
#include "SimpleGeometricObjects.h"
#include "Utilities/StringUtilities.h"
#include <math.h>
#include "ObjectManagers/TableManager.h"


MeshGenerator::MeshGenerator():
  m_dim(0),
  m_min(),
  m_max()
{

/*
   for( int i=0 ; i<3 ; ++i )
   {
    m_wExtensionMin[i] = 0;
    m_wExtensionMax[i] = 0;
    m_nExtensionLayersMin[i] = 0;
    m_nExtensionLayersMax[i] = 0;
    m_commonRatioMin[i] = 1.5;
    m_commonRatioMax[i] = 1.5;
   }
 */
  m_dim = 3;
}

MeshGenerator::~MeshGenerator()
{
  // TODO Auto-generated destructor stub
}

/**
 * @author settgast, fu
 * @param hdn
 */
void MeshGenerator::ReadXML( TICPP::HierarchicalDataNode& hdn )
{

  m_elementType = hdn.GetStringVector("elementTypes");
  if (m_elementType.size() == 0)
    m_elementType = hdn.GetStringVector("elementType");
  if (m_elementType.size() == 0)
    throw GPException("MeshGenerator: No element type is specified");

  if (m_elementType[0] == "C3D8" || m_elementType[0] == "C3D4" || m_elementType[0] == "C3D6")
  {
    m_dim = 3;
  }
  else if (m_elementType[0] == "CPE4" || m_elementType[0] == "STRI" )
  {
    m_dim = 2;
  }
  else
  {
    throw GPException("MeshGenerator: incorrect element type!");
  }


  m_vertices[0] = hdn.GetAttributeVector<realT>("xcoords");
  m_nElems[0] = hdn.GetAttributeVector<int>("nx");
  m_nElemBias[0] = hdn.GetAttributeVectorOrDefault<realT>("xbias", std::vector<realT>(m_nElems[0].size(), 0.0));

  if( m_dim > 1 )
  {
    m_vertices[1] = hdn.GetAttributeVector<realT>("ycoords");
    m_nElems[1] = hdn.GetAttributeVector<int>("ny");
    m_nElemBias[1] = hdn.GetAttributeVectorOrDefault<realT>("ybias", std::vector<realT>(m_nElems[1].size(), 0.0));
  }
  else
  {
    m_vertices[1].push_back(0.0);
    m_vertices[1].push_back(0.0);
    m_nElems[1].push_back(1);
    m_nElemBias[1].push_back(0);
  }

  if( m_dim > 2 )
  {
    m_vertices[2] = hdn.GetAttributeVector<realT>("zcoords");
    m_nElems[2] = hdn.GetAttributeVector<int>("nz");
    m_nElemBias[2] = hdn.GetAttributeVectorOrDefault<realT>("zbias", std::vector<realT>(m_nElems[2].size(), 0.0));
  }
  else
  {
    m_vertices[2].push_back(0.0);
    m_vertices[2].push_back(0.0);
    m_nElems[2].push_back(1);
    m_nElemBias[2].push_back(0);
  }

  {
    bool failFlag = false;
    for( int i=0 ; i<m_dim ; ++i )
    {
      failFlag += ( m_nElems[i].size() != m_vertices[i].size()-1 );
    }
    if( failFlag )
    {
      throw GPException("vertex/element mismatch MeshGenerator::ReadXML()");
    }
  }

  m_numElePerBox.resize(m_nElems[0].size() * m_nElems[1].size() * m_nElems[2].size());

  if (m_elementType.size() != m_numElePerBox.size())
  {
    if (m_elementType.size() == 1)
    {
      m_elementType.resize(m_numElePerBox.size());
      m_elementType = m_elementType[0];
    }
    else
    {
      throw GPException("MeshGenerator: The number of element types is inconsistent with the number of total block.");
    }
  }


  for (localIndex i = 0 ; i < m_elementType.size() ; ++i)
  {
    if (m_elementType[i] == "C3D8")
    {
      m_numElePerBox[i] = 1;
      m_dim = 3;
    }
    else if (m_elementType[i] == "C3D4")
    {
      m_numElePerBox[i] = 6;
      m_dim = 3;
    }
    else if (m_elementType[i] == "C3D6")
    {
      m_numElePerBox[i] = 2;
      m_dim = 3;
      m_trianglePattern = hdn.GetAttributeOrDefault<int>("trianglePattern", 0);
    }
    else if ( m_elementType[i] == "CPE4")
    {
      m_numElePerBox[i] = 1;
      m_dim = 2;
    }
    else if (m_elementType[i] == "STRI")
    {
      m_numElePerBox[i] = 2;
      m_trianglePattern = hdn.GetAttributeOrDefault<int>("trianglePattern", 0);
      m_dim = 2;
    }
  }


  m_regionNames = hdn.GetStringVector("regionNames");
  ExpandMultipleTokens(m_regionNames);
  {
    unsigned int numBlocks = 1;
    for( int i=0 ; i<m_dim ; ++i )
    {
      numBlocks *= m_nElems[i].size();
    }
    if( numBlocks != m_regionNames.size() )
    {
      if (m_regionNames.size() == 1)
      {
        m_regionNames.resize(numBlocks);
        m_regionNames = m_regionNames[0];
      }
      else
      {
        throw GPException("Incorrect number of regionLayout entries specified in MeshGenerator::ReadXML()");
      }
    }
  }

  for( int i=0 ; i<3 ; ++i )
  {
    m_min[i] = m_vertices[i].front();
    m_max[i] = m_vertices[i].back();
  }

  for( int dir=0 ; dir<3 ; ++dir )
  {
    m_firstElemIndexForBlock[dir].resize( m_nElems[dir].size() );
    m_lastElemIndexForBlock[dir].resize( m_nElems[dir].size() );
    m_firstElemIndexForBlock[dir][0] = 0;
    m_lastElemIndexForBlock[dir][0] = m_nElems[dir][0]-1;
    for( unsigned int block=1 ; block<m_nElems[dir].size() ; ++block )
    {
      m_firstElemIndexForBlock[dir][block] = m_lastElemIndexForBlock[dir][block-1] + 1;
      m_lastElemIndexForBlock[dir][block] = m_firstElemIndexForBlock[dir][block] + m_nElems[dir][block]-1;
    }
  }


  m_fPerturb = hdn.GetAttributeOrDefault<realT>("perturbationFactor", 0.0);
  m_randSeed = hdn.GetAttributeOrDefault<int>("perturbationSeed", time(NULL));
  srand(m_randSeed);

  m_mapToRadial = hdn.GetAttributeOrDefault<int>("mapToRadial", 0);

  m_skewAngle = hdn.GetAttributeOrDefault<realT>("skewAngle", 0.0);
  m_skewAngle *= 3.14159265/180;
  R1Tensor zeroVector;
  zeroVector *= 0.0;
  m_skewCenter = hdn.GetAttributeOrDefault<R1Tensor>("skewCenter", zeroVector);


  // Mesh deformation
  m_delayMeshDeformation = hdn.GetAttributeOrDefault<int>("delayMeshDeformation", 0);
  m_meshDx = hdn.GetAttributeString("dxTable");
  m_meshDy = hdn.GetAttributeString("dyTable");
  m_meshDz = hdn.GetAttributeString("dzTable");

}

/**
 * @author settgast
 * @param domain
 */
void MeshGenerator::GenerateElementRegions( PhysicalDomainT * domain )
{
  lvector numElements;

  for( array<string>::size_type r=0 ; r<m_regionNames.size() ; ++r )
  {
    numElements.push_back( 0 );
  }

  domain->m_feElementManager.resize( numElements, m_regionNames, m_elementType );

}

/**
 * @author settgast, fu, sherman
 * @param partition
 * @param domain
 */
void MeshGenerator::GenerateMesh( SpatialPartition& partition,
                                  PhysicalDomainT * domain )
{
  // special case
  bool isRadialWithOneThetaPartition = (m_mapToRadial > 0) && (partition.GetPartitions()[1]==1);


  // partition based on even spacing to get load balance
  // Partition geometrical boundaries will be corrected in the end.
  {
    R1Tensor temp1(m_min);
    R1Tensor temp2(m_max);
    partition.setSizes( temp1, temp2 );
  }

  // find elemCenters for even uniform element sizes
  array<array<real64> > elemCenterCoords(3);
  for( int i=0 ; i<3 ; ++i )
  {
    m_numElemsTotal[i] = 0;
    for( unsigned int block=0 ; block<m_nElems[i].size() ; ++block )
    {
      m_numElemsTotal[i] += m_nElems[i][block];
    }

    elemCenterCoords[i].resize(m_numElemsTotal[i]);
    array<real64> elemCenterCoordsLocal(m_numElemsTotal[i]);
    for( int k=0 ; k<m_numElemsTotal[i] ; ++k )
    {
      elemCenterCoordsLocal[k] = m_min[i] + ( m_max[i] - m_min[i] ) * (k+0.5) / m_numElemsTotal[i];
    }
    MPI_Allreduce( elemCenterCoordsLocal.data(),
                   elemCenterCoords[i].data(),
                   m_numElemsTotal[i],
                   MPI_DOUBLE,
                   MPI_MAX,
                   MPI_COMM_WORLD );
  }

  // find starting/ending index

  // get the first and last indices in this partition each direction
  int firstElemIndexInPartition[3] = {-1,-1,-1};
  int lastElemIndexInPartition[3] = {-2,-2,-2};

  for( int i=0 ; i<3 ; ++i )
  {
//    firstElemIndexInPartition[i] = -1;
//    lastElemIndexInPartition[i] = -2;
    for( int k=0 ; k<m_numElemsTotal[i] ; ++k )
    {
      if( partition.IsCoordInPartition( elemCenterCoords[i][k], i ) )
      {
        firstElemIndexInPartition[i] = k;
        break;
      }
    }

    if( firstElemIndexInPartition[i] > -1 )
    {
      for( int k=firstElemIndexInPartition[i] ; k<m_numElemsTotal[i] ; ++k )
      {
        if( partition.IsCoordInPartition( elemCenterCoords[i][k], i ) )
        {
          lastElemIndexInPartition[i] = k;
        }
      }
    }
  }



  // calculate number of elements in this partition from each region, and the
  // total number of nodes


  std::map< std::string, int > numElemsInRegions;
  std::map< std::string, std::string > elemTypeInRegions;



  array<integer> firstElemIndexForBlockInPartition[3];
  array<integer> lastElemIndexForBlockInPartition[3];

  for( int dir=0 ; dir<3 ; ++dir )
  {
    firstElemIndexForBlockInPartition[dir] = m_firstElemIndexForBlock[dir];
    lastElemIndexForBlockInPartition[dir] = m_lastElemIndexForBlock[dir];

    for( unsigned int block=0 ; block<m_nElems[dir].size() ; ++block )
    {
      if( firstElemIndexForBlockInPartition[dir][block] > lastElemIndexInPartition[dir] ||
          lastElemIndexForBlockInPartition[dir][block] < firstElemIndexInPartition[dir] )
      {
        firstElemIndexForBlockInPartition[dir][block] = -1;
        lastElemIndexForBlockInPartition[dir][block] = -2;
      }
      else
      {
        if( firstElemIndexForBlockInPartition[dir][block] < firstElemIndexInPartition[dir] )
        {
          firstElemIndexForBlockInPartition[dir][block] = firstElemIndexInPartition[dir];
        }
        if( lastElemIndexForBlockInPartition[dir][block] > lastElemIndexInPartition[dir] )
        {
          lastElemIndexForBlockInPartition[dir][block] = lastElemIndexInPartition[dir];
        }
      }
    }
  }


  // TODO This needs to be rewritten for dimensions lower than 3.
  array<string>::const_iterator iterRegion=m_regionNames.begin();
  for( unsigned int iblock=0 ; iblock<m_nElems[0].size() ; ++iblock )
  {
    for( unsigned int jblock=0 ; jblock<m_nElems[1].size() ; ++jblock )
    {
      for( unsigned int kblock=0 ; kblock<m_nElems[2].size() ; ++kblock, ++iterRegion )
      {
        numElemsInRegions[*iterRegion] = 0;
        elemTypeInRegions[*iterRegion] = "";
      }
    }
  }

  iterRegion=m_regionNames.begin();
  {
    localIndex iR = 0;
    for( unsigned int iblock=0 ; iblock<m_nElems[0].size() ; ++iblock )
    {
      for( unsigned int jblock=0 ; jblock<m_nElems[1].size() ; ++jblock )
      {
        for( unsigned int kblock=0 ; kblock<m_nElems[2].size() ; ++kblock, ++iterRegion, ++iR )
        {
          int numElemsInRegion = 1;
          numElemsInRegion *= lastElemIndexForBlockInPartition[0][iblock] - firstElemIndexForBlockInPartition[0][iblock] + 1;

          if( m_dim > 1 )
          {
            numElemsInRegion *= lastElemIndexForBlockInPartition[1][jblock] - firstElemIndexForBlockInPartition[1][jblock] + 1;
          }
          if( m_dim > 2 )
          {
            numElemsInRegion *= lastElemIndexForBlockInPartition[2][kblock] - firstElemIndexForBlockInPartition[2][kblock] + 1;
          }

          numElemsInRegion *= m_numElePerBox[iR];
          numElemsInRegions[*iterRegion] += numElemsInRegion;
          elemTypeInRegions[*iterRegion] = m_elementType[iR];

        }
      }
    }
  }



  localIndex numNodes = 1;



//  int numElemsInDir[3] = {1,1,1};
  int numNodesInDir[3] = {1,1,1};

  for( int i=0 ; i<m_dim ; ++i )
  {
//    numElemsInDir[i] = lastElemIndexInPartition[i] -
// firstElemIndexInPartition[i] + 1;
    numNodesInDir[i] = lastElemIndexInPartition[i] - firstElemIndexInPartition[i] + 2;
    if(isRadialWithOneThetaPartition && i==1)
    {
      numNodesInDir[1] -= 1;
    }
    numNodes *= numNodesInDir[i];
  }

  domain->m_feNodeManager.resize( numNodes );
  {
    localIndex localNodeIndex = 0;
    for( int i=0 ; i<numNodesInDir[0] ; ++i )
    {
      for( int j=0 ; j<numNodesInDir[1] ; ++j )
      {
        for( int k=0 ; k<numNodesInDir[2] ; ++k )
        {
          int index[3] = {i,j,k};
          for( int a=0 ; a<m_dim ; ++a )
          {
            index[a] += firstElemIndexInPartition[a];
          }

          (*domain->m_feNodeManager.m_refposition)[localNodeIndex] = NodePosition( index, m_trianglePattern );

          // alter global node map for radial mesh
          if (m_mapToRadial > 0)
          {
            if (isEqual( (*domain->m_feNodeManager.m_refposition)[localNodeIndex][1], m_max[1], 1e-10 ))
            {
              index[1] = 0;
            }
          }

          domain->m_feNodeManager.m_localToGlobalMap[localNodeIndex] = NodeGlobalIndex( index );

          // cartesian-specific nodesets
          if (m_mapToRadial == 0)
          {
            if( isEqual( (*domain->m_feNodeManager.m_refposition)[localNodeIndex][0], m_min[0], 1e-10 ) )
            {
              domain->m_feNodeManager.m_Sets["xneg"].insert(localNodeIndex);
            }
            if( isEqual( (*domain->m_feNodeManager.m_refposition)[localNodeIndex][0], m_max[0], 1e-10 ) )
            {
              domain->m_feNodeManager.m_Sets["xpos"].insert(localNodeIndex);
            }
            if( isEqual( (*domain->m_feNodeManager.m_refposition)[localNodeIndex][1], m_min[1], 1e-10 ) )
            {
              domain->m_feNodeManager.m_Sets["yneg"].insert(localNodeIndex);
            }
            if( isEqual( (*domain->m_feNodeManager.m_refposition)[localNodeIndex][1], m_max[1], 1e-10 ) )
            {
              domain->m_feNodeManager.m_Sets["ypos"].insert(localNodeIndex);
            }
          }
          else
          {
            // radial-specific nodesets
            domain->m_feNodeManager.m_Sets["rin"]; // generate the empty sets
            domain->m_feNodeManager.m_Sets["rout"];
            if( isEqual( (*domain->m_feNodeManager.m_refposition)[localNodeIndex][0], m_min[0], 1e-10 ) )
            {
              domain->m_feNodeManager.m_Sets["rin"].insert(localNodeIndex);
            }
            if( isEqual( (*domain->m_feNodeManager.m_refposition)[localNodeIndex][0], m_max[0], 1e-10 ) )
            {
              domain->m_feNodeManager.m_Sets["rout"].insert(localNodeIndex);
            }
          }

          // general nodesets
          if( isEqual( (*domain->m_feNodeManager.m_refposition)[localNodeIndex][2], m_min[2], 1e-10 ) )
          {
            domain->m_feNodeManager.m_Sets["zneg"].insert(localNodeIndex);
          }
          if( isEqual( (*domain->m_feNodeManager.m_refposition)[localNodeIndex][2], m_max[2], 1e-10 ) )
          {
            domain->m_feNodeManager.m_Sets["zpos"].insert(localNodeIndex);
          }
          domain->m_feNodeManager.m_Sets["all"].insert(localNodeIndex);


          ++localNodeIndex;

        }
      }
    }
  }


  {
    lvector numElements;
    array<string> elementRegionNames;
    array<string> elementTypes;
    std::map<std::string,localIndex> localElemIndexInRegion;

    for( std::map< std::string, int >::iterator iterNumElemsInRegion=numElemsInRegions.begin() ;
         iterNumElemsInRegion!=numElemsInRegions.end() ; ++iterNumElemsInRegion )
    {

      numElements.push_back( iterNumElemsInRegion->second );
      elementRegionNames.push_back( iterNumElemsInRegion->first );
      elementTypes.push_back(elemTypeInRegions[iterNumElemsInRegion->first]);

      localElemIndexInRegion[iterNumElemsInRegion->first] = 0;
    }
    domain->m_feElementManager.resize( numElements, elementRegionNames, elementTypes );

    // assign global numbers to elements
    iterRegion=m_regionNames.begin();
    set<std::string> processedRegionNames;
    localIndex iR = 0;

    for( unsigned int iblock=0 ; iblock<m_nElems[0].size() ; ++iblock )
    {
      for( unsigned int jblock=0 ; jblock<m_nElems[1].size() ; ++jblock )
      {
        for( unsigned int kblock=0 ; kblock<m_nElems[2].size() ; ++kblock, ++iterRegion, ++iR )
        {
          ElementRegionT& elemRegion = domain->m_feElementManager.m_ElementRegions[*iterRegion];


          int numElemsInDirForRegion[3] = { lastElemIndexForBlockInPartition[0][iblock] - firstElemIndexForBlockInPartition[0][iblock] + 1,
                                            lastElemIndexForBlockInPartition[1][jblock] - firstElemIndexForBlockInPartition[1][jblock] + 1,
                                            lastElemIndexForBlockInPartition[2][kblock] - firstElemIndexForBlockInPartition[2][kblock] + 1 };


          for( int i=0 ; i<numElemsInDirForRegion[0] ; ++i )
          {
            for( int j=0 ; j<numElemsInDirForRegion[1] ; ++j )
            {
              for( int k=0 ; k<numElemsInDirForRegion[2] ; ++k )
              {
                int index[3] = {i + firstElemIndexForBlockInPartition[0][iblock],
                                j + firstElemIndexForBlockInPartition[1][jblock],
                                k + firstElemIndexForBlockInPartition[2][kblock]};

                const localIndex firstNodeIndex = numNodesInDir[1]*numNodesInDir[2]* ( index[0] - firstElemIndexInPartition[0] )
                                                  + numNodesInDir[2]* ( index[1] - firstElemIndexInPartition[1] )
                                                  + ( index[2] - firstElemIndexInPartition[2] );
                int nodeOfBox[8];


                if (m_elementType[iR] == "CPE4" || m_elementType[iR] == "STRI")
                {

                  nodeOfBox[0] = firstNodeIndex;
                  nodeOfBox[1] = numNodesInDir[1]*numNodesInDir[2] + firstNodeIndex;
                  nodeOfBox[2] = numNodesInDir[1]*numNodesInDir[2] + numNodesInDir[2] + firstNodeIndex;
                  nodeOfBox[3] = numNodesInDir[2] + firstNodeIndex;

                }
                else
                {
                  nodeOfBox[0] = firstNodeIndex;
                  nodeOfBox[1] = numNodesInDir[1]*numNodesInDir[2] + firstNodeIndex;
                  nodeOfBox[2] = numNodesInDir[1]*numNodesInDir[2] + numNodesInDir[2] + firstNodeIndex;
                  nodeOfBox[3] = numNodesInDir[2] + firstNodeIndex;

                  nodeOfBox[4] = firstNodeIndex+1;
                  nodeOfBox[5] = numNodesInDir[1]*numNodesInDir[2] + firstNodeIndex+1;
                  nodeOfBox[6] = numNodesInDir[1]*numNodesInDir[2] + numNodesInDir[2] + firstNodeIndex+1;
                  nodeOfBox[7] = numNodesInDir[2] + firstNodeIndex+1;

                  //               7___________________ 6
                  //               /                   /|
                  //              /                   / |f
                  //             /                   /  |
                  //           4/__________________5/   |
                  //            |                   |   |
                  //            |                   |   |
                  //            |                   |   |
                  //            |                   |   |
                  //            |                   |   |
                  //            |   3               |   /2        z
                  //            |                   |  /          |   y
                  //            |                   | /           |  /
                  //            |___________________|/            | /
                  //            0                   1             |/____ x

                }
                // fix local connectivity for single theta (y) partition (radial
                // meshes only)
                if( isRadialWithOneThetaPartition )
                {
                  if(j == numElemsInDirForRegion[1]-1 && jblock == m_nElems[1].size()-1 )  // last
                                                                                           // set
                                                                                           // of
                                                                                           // elements
                  {
                    index[1] = -1;
                    const localIndex firstNodeIndexR = numNodesInDir[1]*numNodesInDir[2]* ( index[0] - firstElemIndexInPartition[0] )
                                                       + numNodesInDir[2]* ( index[1] - firstElemIndexInPartition[1] )
                                                       + ( index[2] - firstElemIndexInPartition[2] );
                    nodeOfBox[2] = numNodesInDir[1]*numNodesInDir[2] + numNodesInDir[2] + firstNodeIndexR;
                    nodeOfBox[3] = numNodesInDir[2] + firstNodeIndexR;
                    nodeOfBox[6] = numNodesInDir[1]*numNodesInDir[2] + numNodesInDir[2] + firstNodeIndexR+1;
                    nodeOfBox[7] = numNodesInDir[2] + firstNodeIndexR+1;
                  }
                }


                for ( int iEle=0 ; iEle<m_numElePerBox[iR] ; ++iEle)
                {
                  localIndex& localElemIndex = localElemIndexInRegion[*iterRegion];
                  elemRegion.m_localToGlobalMap[localElemIndex] = ElemGlobalIndex( index ) * m_numElePerBox[iR] + iEle;

                  array<integer> nodeIDInBox(elemRegion.m_numNodesPerElem);

                  GetElemToNodesRelationInBox ( m_elementType[iR], index, iEle, nodeIDInBox.data(), elemRegion.m_numNodesPerElem);

                  for ( localIndex iN = 0 ; iN < elemRegion.m_numNodesPerElem ; ++iN)
                  {
                    elemRegion.m_toNodesRelation[localElemIndex][iN] = nodeOfBox[nodeIDInBox[iN]];
                  }
                  ++localElemIndex;


                }
              }
            }
          }
        }
      }
    }
  }

  // Correct partition geometrical boundary.
  {
    R1Tensor pMin, pMax;
    partition.getPartitionGeometricalBoundary(pMin, pMax);
    for( int i=0 ; i<m_dim ; ++i )
    {
      realT xMinByNumElems = (pMin[i] - m_min[i]) / (m_max[i] - m_min[i]) * m_numElemsTotal[i];
      realT xMaxByNumElems = (pMax[i] - m_min[i]) / (m_max[i] - m_min[i]) * m_numElemsTotal[i];

      localIndex iBlockMin(0), iBlockMax(0);
      while ( (xMinByNumElems < m_firstElemIndexForBlock[i][iBlockMin] * 1.0 || xMinByNumElems > m_lastElemIndexForBlock[i][iBlockMin] * 1.0 +1.0)
              && iBlockMin < m_nElems[i].size() - 1)
      {
        ++iBlockMin;
      }
      while ( (xMaxByNumElems < m_firstElemIndexForBlock[i][iBlockMax] * 1.0 || xMaxByNumElems > m_lastElemIndexForBlock[i][iBlockMax] * 1.0 +1.0)
              && iBlockMax < m_nElems[i].size() - 1)
      {
        ++iBlockMax;
      }
      //sanity check


      if (xMinByNumElems < m_firstElemIndexForBlock[i][iBlockMin] * 1.0 ||
          xMinByNumElems > m_lastElemIndexForBlock[i][iBlockMin] * 1.0 + 1.0 ||
          iBlockMin >= m_nElems[i].size() )
      {
        std::cout << "WARNING: Had some trouble in correcting partition geometric boundaries.  If this affects the results, contact a developer" << std::endl;
      }
      if (xMaxByNumElems < m_firstElemIndexForBlock[i][iBlockMax] * 1.0 ||
          xMaxByNumElems > m_lastElemIndexForBlock[i][iBlockMax] * 1.0 + 1.0 ||
          iBlockMax >= m_nElems[i].size())
      {
        std::cout << "WARNING: Had some trouble in correcting partition geometric boundaries.  If this affects the results, contact a developer" << std::endl;
      }

      pMin[i] = m_vertices[i][iBlockMin] + (xMinByNumElems - m_firstElemIndexForBlock[i][iBlockMin]) * (m_vertices[i][iBlockMin+1] - m_vertices[i][iBlockMin]) /
                m_nElems[i][iBlockMin];
      pMax[i] = m_vertices[i][iBlockMax] + (xMaxByNumElems - m_firstElemIndexForBlock[i][iBlockMax]) * (m_vertices[i][iBlockMax+1] - m_vertices[i][iBlockMax]) /
                m_nElems[i][iBlockMax];
    }

    partition.SetPartitionGeometricalBoundary(pMin, pMax);

  }


  /*
     {// Move nodes in the extension layers.

     for (localIndex iN = 0; iN != domain->m_feNodeManager.DataLengths(); ++iN)
     {
      for (int i=0; i<m_dim; ++i)
      {
        if ( (*domain->m_feNodeManager.m_refposition)[iN][i] < m_min[i])
        {
          int eLayer = (int) ((m_min[i]
             -(*domain->m_feNodeManager.m_refposition)[iN][i]) / ((m_max[i] -
             m_min[i]) / m_numElems[i]) + 0.5);
          (*domain->m_feNodeManager.m_refposition)[iN][i] = m_min[i] -
             ((m_max[i] - m_min[i]) / m_numElems[i]) * m_commonRatioMin[i] * (1-
             pow(m_commonRatioMin[i], eLayer)) / (1 - m_commonRatioMin[i]);
        }
        else if ((*domain->m_feNodeManager.m_refposition)[iN][i] > m_max[i])
        {
          int eLayer = (int) (((*domain->m_feNodeManager.m_refposition)[iN][i] -
             m_max[i] ) / ((m_max[i] - m_min[i]) / m_numElems[i]) + 0.5);
          (*domain->m_feNodeManager.m_refposition)[iN][i] = m_max[i] +
             ((m_max[i] - m_min[i]) / m_numElems[i]) * m_commonRatioMax[i] * (1-
             pow(m_commonRatioMax[i], eLayer)) / (1 - m_commonRatioMax[i]);
        }

      }
     }
     }
   */


  // Node perturbation
  if (m_fPerturb > 0)
  {
    for (localIndex iN = 0 ; iN != domain->m_feNodeManager.DataLengths() ; ++iN)
    {

      for (int i=0 ; i<m_dim ; ++i)
      {
        if ( (*domain->m_feNodeManager.m_refposition)[iN][i] > m_min[i] && (*domain->m_feNodeManager.m_refposition)[iN][i] < m_max[i] )
        {
          srand(domain->m_feNodeManager.m_localToGlobalMap[iN] + m_randSeed + i); // This
                                                                                  // ensures
                                                                                  // that
                                                                                  // the
                                                                                  // perturbation
                                                                                  // pattern
                                                                                  // is
                                                                                  // unaffected
                                                                                  // by
                                                                                  // domain
                                                                                  // partitioning.
          (*domain->m_feNodeManager.m_refposition)[iN][i] += ((m_max[i] - m_min[i]) / m_numElemsTotal[i]) * ( (rand()*1.0) / RAND_MAX - 0.5) * 2 * m_fPerturb;
        }
      }
    }
  }

  if (std::fabs(m_skewAngle) > 0.0)
  {
    for (localIndex iN = 0 ; iN != domain->m_feNodeManager.DataLengths() ; ++iN)
    {
      (*domain->m_feNodeManager.m_refposition)[iN][0] -= ((*domain->m_feNodeManager.m_refposition)[iN][1] - m_skewCenter[1]) * std::tan(m_skewAngle);
    }
  }

  if (m_mapToRadial > 0)
  {
    // Map to radial mesh
    for (localIndex iN = 0 ; iN != domain->m_feNodeManager.DataLengths() ; ++iN)
    {
      meshTheta = (*domain->m_feNodeManager.m_refposition)[iN][1]*3.141592654/180.0;
      meshAxis = round(meshTheta*2.0/3.141592654);
      meshPhi = fabs(meshTheta - meshAxis*3.141592654/2.0);
      meshRout = m_max[0]/cos(meshPhi);

      if (m_mapToRadial > 1)
      {
        meshRact = ((meshRout - m_min[0])/(m_max[0] - m_min[0]))*((*domain->m_feNodeManager.m_refposition)[iN][0]-m_min[0]) + m_min[0];
      }
      else
      {
        meshRact = (*domain->m_feNodeManager.m_refposition)[iN][0];
      }

      (*domain->m_feNodeManager.m_refposition)[iN][0] = meshRact * cos(meshTheta);
      (*domain->m_feNodeManager.m_refposition)[iN][1] = meshRact * sin(meshTheta);

      // add mapped values to nodesets
      if (m_mapToRadial > 1)
      {
        if( isEqual( (*domain->m_feNodeManager.m_refposition)[iN][0], -1*m_max[0], 1e-6 ) )
        {
          domain->m_feNodeManager.m_Sets["xneg"].insert(iN);
        }
        if( isEqual( (*domain->m_feNodeManager.m_refposition)[iN][0], m_max[0], 1e-6 ) )
        {
          domain->m_feNodeManager.m_Sets["xpos"].insert(iN);
        }
        if( isEqual( (*domain->m_feNodeManager.m_refposition)[iN][1], -1*m_max[0], 1e-6 ) )
        {
          domain->m_feNodeManager.m_Sets["yneg"].insert(iN);
        }
        if( isEqual( (*domain->m_feNodeManager.m_refposition)[iN][1], m_max[0], 1e-6 ) )
        {
          domain->m_feNodeManager.m_Sets["ypos"].insert(iN);
        }
      }
    }
  }

  if (m_delayMeshDeformation == 0)
  {
    RemapMesh(domain);
  }
}


/**
 * @author fu
 * @param elementType
 * @param index
 * @param iEle
 * @param nodeIDInBox
 * @param size
 */
void MeshGenerator::GetElemToNodesRelationInBox ( const std::string& elementType,
                                                  const int index[],
                                                  const int& iEle,
                                                  int nodeIDInBox[],
                                                  const int size)

{
  if (elementType == "C3D8")
  {
    nodeIDInBox[0] = 0;
    nodeIDInBox[1] = 1;
    nodeIDInBox[2] = 3;
    nodeIDInBox[3] = 2;
    nodeIDInBox[4] = 4;
    nodeIDInBox[5] = 5;
    nodeIDInBox[6] = 7;
    nodeIDInBox[7] = 6;
  }
  else if ( (elementType == "C3D6") && (m_trianglePattern == 0) )
  {
    if ( (index[0] + index[1] ) % 2  ==  1)
    {
      if  ( iEle % 2 ==0)
      {
        nodeIDInBox[0] = 6;
        nodeIDInBox[1] = 2;
        nodeIDInBox[2] = 4;
        nodeIDInBox[3] = 0;
        nodeIDInBox[4] = 7;
        nodeIDInBox[5] = 3;
        nodeIDInBox[6] = 7;
        nodeIDInBox[7] = 3;
      }
      else
      {
        nodeIDInBox[0] = 5;
        nodeIDInBox[1] = 1;
        nodeIDInBox[2] = 4;
        nodeIDInBox[3] = 0;
        nodeIDInBox[4] = 6;
        nodeIDInBox[5] = 2;
        nodeIDInBox[6] = 6;
        nodeIDInBox[7] = 2;
      }
    }
    else
    {
      if  ( iEle % 2 ==0)
      {
        nodeIDInBox[0] = 5;
        nodeIDInBox[1] = 1;
        nodeIDInBox[2] = 4;
        nodeIDInBox[3] = 0;
        nodeIDInBox[4] = 7;
        nodeIDInBox[5] = 3;
        nodeIDInBox[6] = 7;
        nodeIDInBox[7] = 3;
      }
      else
      {
        nodeIDInBox[0] = 5;
        nodeIDInBox[1] = 1;
        nodeIDInBox[2] = 7;
        nodeIDInBox[3] = 3;
        nodeIDInBox[4] = 6;
        nodeIDInBox[5] = 2;
        nodeIDInBox[6] = 6;
        nodeIDInBox[7] = 2;
      }
    }
  }
  else if ( (elementType == "C3D6") && (m_trianglePattern == 1) )
  {
    if ( index[1] % 2  ==  0)
    {
      if  ( iEle % 2 ==0)
      {
        nodeIDInBox[0] = 5;
        nodeIDInBox[1] = 1;
        nodeIDInBox[2] = 4;
        nodeIDInBox[3] = 0;
        nodeIDInBox[4] = 6;
        nodeIDInBox[5] = 2;
        nodeIDInBox[6] = 6;
        nodeIDInBox[7] = 2;
      }
      else
      {
        nodeIDInBox[0] = 6;
        nodeIDInBox[1] = 2;
        nodeIDInBox[2] = 4;
        nodeIDInBox[3] = 0;
        nodeIDInBox[4] = 7;
        nodeIDInBox[5] = 3;
        nodeIDInBox[6] = 7;
        nodeIDInBox[7] = 3;
      }
    }
    else
    {
      if  ( iEle % 2 ==0)
      {
        nodeIDInBox[0] = 5;
        nodeIDInBox[1] = 1;
        nodeIDInBox[2] = 4;
        nodeIDInBox[3] = 0;
        nodeIDInBox[4] = 7;
        nodeIDInBox[5] = 3;
        nodeIDInBox[6] = 7;
        nodeIDInBox[7] = 3;
      }
      else
      {
        nodeIDInBox[0] = 5;
        nodeIDInBox[1] = 1;
        nodeIDInBox[2] = 7;
        nodeIDInBox[3] = 3;
        nodeIDInBox[4] = 6;
        nodeIDInBox[5] = 2;
        nodeIDInBox[6] = 6;
        nodeIDInBox[7] = 2;
      }

    }
  }
  else if (elementType == "CPE4")
  {
    nodeIDInBox[0] = 0;
    nodeIDInBox[1] = 1;
    nodeIDInBox[2] = 3;
    nodeIDInBox[3] = 2;
  }
  else if ( (elementType == "STRI") && (m_trianglePattern == 0) )
  {
    if ( (index[0] + index[1] ) % 2  ==  1)
    {
      if  ( iEle % 2 ==0)
      {
        nodeIDInBox[0] = 0;
        nodeIDInBox[1] = 2;
        nodeIDInBox[2] = 3;
      }
      else
      {
        nodeIDInBox[0] = 0;
        nodeIDInBox[1] = 1;
        nodeIDInBox[2] = 2;
      }
    }
    else
    {
      if  ( iEle % 2 ==0)
      {
        nodeIDInBox[0] = 0;
        nodeIDInBox[1] = 1;
        nodeIDInBox[2] = 3;
      }
      else
      {
        nodeIDInBox[0] = 1;
        nodeIDInBox[1] = 2;
        nodeIDInBox[2] = 3;
      }

    }
  }
  else if ( (elementType == "STRI") && (m_trianglePattern == 1) )
  {
    if ( index[1] % 2  ==  0)
    {
      if  ( iEle % 2 ==0)
      {
        nodeIDInBox[0] = 0;
        nodeIDInBox[1] = 1;
        nodeIDInBox[2] = 2;
      }
      else
      {
        nodeIDInBox[0] = 2;
        nodeIDInBox[1] = 3;
        nodeIDInBox[2] = 0;
      }
    }
    else
    {
      if  ( iEle % 2 ==0)
      {
        nodeIDInBox[0] = 0;
        nodeIDInBox[1] = 1;
        nodeIDInBox[2] = 3;
      }
      else
      {
        nodeIDInBox[0] = 1;
        nodeIDInBox[1] = 2;
        nodeIDInBox[2] = 3;
      }

    }
  }
  else if (elementType == "C3D4")
  {
    int mapBoxTet[8][6][4] =
    {
      {
        { 0, 3, 7, 6 }, { 0, 7, 4, 6 }, { 0, 5, 1, 6 }, { 0, 4, 5, 6 }, { 0, 1, 2, 6 }, { 0, 2, 3, 6 }
      },
      {
        { 1, 5, 6, 7 }, { 1, 6, 2, 7 }, { 0, 4, 1, 7 }, { 1, 4, 5, 7 }, { 0, 1, 3, 7 }, { 1, 2, 3, 7 }
      },
      {
        { 0, 3, 4, 2 }, { 3, 7, 4, 2 }, { 0, 4, 1, 2 }, { 1, 4, 5, 2 }, { 4, 7, 6, 2 }, { 4, 6, 5, 2 }
      },
      {
        { 1, 5, 2, 3 }, { 2, 5, 6, 3 }, { 0, 5, 1, 3 }, { 0, 4, 5, 3 }, { 4, 7, 5, 3 }, { 5, 7, 6, 3 }
      },
      {
        { 0, 3, 4, 5 }, { 3, 7, 4, 5 }, { 3, 6, 7, 5 }, { 3, 2, 6, 5 }, { 3, 0, 1, 5 }, { 1, 2, 3, 5 }
      },
      {
        { 1, 5, 2, 4 }, { 2, 5, 6, 4 }, { 2, 6, 7, 4 }, { 2, 7, 3, 4 }, { 0, 2, 3, 4 }, { 0, 1, 2, 4 }
      },
      {
        { 0, 7, 4, 1 }, { 0, 3, 7, 1 }, { 2, 7, 3, 1 }, { 2, 6, 7, 1 }, { 4, 7, 5, 1 }, { 5, 7, 6, 1 }
      },
      {
        { 1, 5, 6, 0 }, { 1, 6, 2, 0 }, { 2, 6, 3, 0 }, { 3, 6, 7, 0 }, { 4, 6, 5, 0 }, { 4, 7, 6, 0 }
      }
    };


    int mapBoxType[2][2][2];
    mapBoxType[0][0][0]=0;
    mapBoxType[1][0][0]=1;
    mapBoxType[0][0][1]=2;
    mapBoxType[1][0][1]=3;
    mapBoxType[0][1][0]=4;
    mapBoxType[1][1][0]=5;
    mapBoxType[0][1][1]=6;
    mapBoxType[1][1][1]=7;

    int boxType = mapBoxType[index[0]%2][index[1]%2][index[2]%2];
    for (int i=0 ; i<size ; ++i)
    {
      nodeIDInBox[i] = mapBoxTet[boxType][iEle][i];
    }

  }
}


void MeshGenerator::RemapMesh ( PhysicalDomainT * domain )
{
  // Node mapping
  if (!m_meshDx.empty())
  {
    const Table3D* tableDx = stlMapLookupPointer(TableManager::Instance().Tables<3>(), m_meshDx);

    for (localIndex iN=0 ; iN!=domain->m_feNodeManager.DataLengths() ; ++iN)
    {
      realT dx=tableDx->Lookup((*domain->m_feNodeManager.m_refposition)[iN]);
      (*domain->m_feNodeManager.m_refposition)[iN][0] += dx;
    }
  }

  if (!m_meshDy.empty())
  {
    const Table3D* tableDy = stlMapLookupPointer(TableManager::Instance().Tables<3>(), m_meshDy);

    for (localIndex iN=0 ; iN!=domain->m_feNodeManager.DataLengths() ; ++iN)
    {
      realT dy=tableDy->Lookup((*domain->m_feNodeManager.m_refposition)[iN]);
      (*domain->m_feNodeManager.m_refposition)[iN][1] += dy;
    }
  }

  if (!m_meshDz.empty())
  {
    const Table3D* tableDz = stlMapLookupPointer(TableManager::Instance().Tables<3>(), m_meshDz);

    for (localIndex iN=0 ; iN!=domain->m_feNodeManager.DataLengths() ; ++iN)
    {
      realT dz=tableDz->Lookup((*domain->m_feNodeManager.m_refposition)[iN]);
      (*domain->m_feNodeManager.m_refposition)[iN][2] += dz;
    }
  }

}
