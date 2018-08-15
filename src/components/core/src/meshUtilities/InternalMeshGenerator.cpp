/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
 * InternalMeshGenerator.cpp
 *
 *  Created on: Nov 19, 2012
 *      Author: settgast1
 */

#include "InternalMeshGenerator.hpp"

#include "managers/DomainPartition.hpp"

#include "codingUtilities/StringUtilities.hpp"
#include <math.h>
#include <algorithm>
//#include "managers/TableManager.hpp"
//#include "SimpleGeometricObjects.hpp"

#ifdef USE_ATK
#include "slic/slic.hpp"
#endif

#include "MPI_Communications/PartitionBase.hpp"
#include "MPI_Communications/SpatialPartition.hpp"

#include "mesh/MeshBody.hpp"

namespace geosx
{
using namespace dataRepository;

InternalMeshGenerator::InternalMeshGenerator( string const & name, ManagedGroup * const parent ):
  MeshGeneratorBase( name, parent ),
//    m_vertices({this->RegisterViewWrapper<real64_array>(keys::xCoords).reference(),
//                this->RegisterViewWrapper<real64_array>(keys::yCoords).reference(),
//                this->RegisterViewWrapper<real64_array>(keys::zCoords).reference()
// }),
  m_dim( 0 ),
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

InternalMeshGenerator::~InternalMeshGenerator()
{
  // TODO Auto-generated destructor stub
}

void InternalMeshGenerator::FillDocumentationNode()
{
  //MeshLevel * const mesh =
  // domain->group_cast<DomainPartition*>()->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  //NodeManager * const nodes    = mesh->getNodeManager();
  // CellBlockManager * elems =
  // domain->GetGroup<CellBlockManager>(keys::cellManager);

  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName( "InternalMesh" );
  docNode->setSchemaType( "Node" );
  docNode->setShortDescription( "a mesh generator" );


//  nodes->getDocumentationNode()->AllocateChildNode( keys::ReferencePosition,
//                                                   keys::ReferencePosition,
//                                                   -1,
//                                                   "r1_array",
//                                                   "r1_array",
//                                                   "Reference position of mesh
// vertex points",
//                                                   "Reference position of mesh
// vertex points",
//                                                   "1",
//                                                   "",
//                                                   1,
//                                                   0,
//                                                   0 );

  docNode->AllocateChildNode( keys::xCoords,
                              keys::xCoords,
                              -1,
                              "real64_array",
                              "real64_array",
                              "x-coordinates of mesh vertex points",
                              "x-coordinates of mesh vertex points",
                              "1",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::yCoords,
                              keys::yCoords,
                              -1,
                              "real64_array",
                              "real64_array",
                              "y-coordinates of mesh vertex points",
                              "y-coordinates of mesh vertex points",
                              "1",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::zCoords,
                              keys::zCoords,
                              -1,
                              "real64_array",
                              "real64_array",
                              "z-coordinates of mesh vertex points",
                              "z-coordinates of mesh vertex points",
                              "1",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::xElems,
                              keys::xElems,
                              -1,
                              "integer_array",
                              "integer_array",
                              "number of elements in x-direction",
                              "number of elements in x-direction",
                              "1",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::yElems,
                              keys::yElems,
                              -1,
                              "integer_array",
                              "integer_array",
                              "number of elements in y-direction",
                              "number of elements in y-direction",
                              "1",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::zElems,
                              keys::zElems,
                              -1,
                              "integer_array",
                              "integer_array",
                              "number of elements in z-direction",
                              "number of elements in z-direction",
                              "1",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::xBias,
                              keys::xBias,
                              -1,
                              "real64_array",
                              "real64_array",
                              "spacing bias in x-direction",
                              "spacing bias in x-direction",
                              "0",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::yBias,
                              keys::yBias,
                              -1,
                              "real64_array",
                              "real64_array",
                              "spacing bias in y-direction",
                              "spacing bias in y-direction",
                              "0",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::zBias,
                              keys::zBias,
                              -1,
                              "real64_array",
                              "real64_array",
                              "spacing bias in z-direction",
                              "spacing bias in z-direction",
                              "0",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::cellBlockNames,
                              keys::cellBlockNames,
                              -1,
                              "string_array",
                              "string_array",
                              "names of the regions",
                              "names of the regions",
                              "Region",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::elementTypes,
                              keys::elementTypes,
                              -1,
                              "string_array",
                              "string_array",
                              "topology of discrete volumes",
                              "topology of discrete volumes",
                              "C3D8",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::trianglePattern,
                              keys::trianglePattern,
                              -1,
                              "integer",
                              "integer",
                              "",
                              "",
                              "0",
                              "",
                              0,
                              1,
                              0 );

}

//}
/**
 * @author settgast
 * @param domain
 */
void InternalMeshGenerator::GenerateElementRegions( DomainPartition& domain )
{
  //  lvector numElements;
  //
  //  for( array1d<string>::size_type r=0 ; r<m_regionNames.size() ; ++r )
  //  {
  //    numElements.push_back( 0 );
  //  }
  //
  //  domain.m_feElementManager->resize( numElements, m_regionNames,
  // m_elementType );

}

void InternalMeshGenerator::ReadXML_PostProcess()
{

  real64_array const &  xCoords = this->getReference<real64_array>(keys::xCoords);
  real64_array const &  yCoords = this->getReference<real64_array>(keys::yCoords);
  real64_array const &  zCoords = this->getReference<real64_array>(keys::zCoords);
  m_vertices[0] = xCoords;
  m_vertices[1] = yCoords;
  m_vertices[2] = zCoords;

  integer_array const &  xElems = this->getReference<integer_array>(keys::xElems);
  integer_array const &  yElems = this->getReference<integer_array>(keys::yElems);
  integer_array const &  zElems = this->getReference<integer_array>(keys::zElems);
  m_nElems[0] = xElems;
  m_nElems[1] = yElems;
  m_nElems[2] = zElems;

  real64_array const &  xBias = this->getReference<real64_array>(keys::xBias);
  real64_array const &  yBias = this->getReference<real64_array>(keys::yBias);
  real64_array const &  zBias = this->getReference<real64_array>(keys::zBias);
  m_nElemBias[0] = xBias;
  m_nElemBias[1] = yBias;
  m_nElemBias[2] = zBias;

  m_regionNames = this->getReference<string_array>(keys::cellBlockNames);
  m_elementType = this->getReference<string_array>(keys::elementTypes);
  m_trianglePattern = *(this->getData<integer>(keys::trianglePattern));



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
    GEOS_ERROR("InternalMeshGenerator: incorrect element type!");
  }

  {
    bool failFlag = false;
    for( int i=0 ; i<m_dim ; ++i )
    {
      failFlag += ( m_nElems[i].size() != m_vertices[i].size()-1 );
    }
    if( failFlag )
    {
      GEOS_ERROR("vertex/element mismatch InternalMeshGenerator::ReadXMLPost()");
    }
  }

  m_numElePerBox.resize(m_nElems[0].size() * m_nElems[1].size() * m_nElems[2].size());

  if (static_cast<long>(m_elementType.size()) != m_numElePerBox.size())
  {
    if (m_elementType.size() == 1)
    {
      m_elementType.resize(m_numElePerBox.size());
      std::fill(m_elementType.begin(), m_elementType.end(), m_elementType[0]);
    }
    else
    {
#ifdef USE_ATK
      SLIC_ERROR("InternalMeshGenerator: The number of element types is inconsistent with the number of total block.");
#endif
    }
  }


  for (localIndex i = 0 ; i < static_cast<localIndex>( m_elementType.size() ) ; ++i)
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
    }
    else if ( m_elementType[i] == "CPE4")
    {
      m_numElePerBox[i] = 1;
      m_dim = 2;
    }
    else if (m_elementType[i] == "STRI")
    {
      m_numElePerBox[i] = 2;
      m_dim = 2;
    }
  }


//    ExpandMultipleTokens(m_regionNames);
  {
    int numBlocks = 1;
    for( int i=0 ; i<m_dim ; ++i )
    {
      numBlocks *= m_nElems[i].size();
    }
    if( numBlocks != static_cast<int>(m_regionNames.size()) )
    {
      if (m_regionNames.size() == 1)
      {
        m_regionNames.resize(numBlocks);
        std::fill(m_regionNames.begin(), m_regionNames.end(), m_regionNames[0]);
      }
      else
      {
#ifdef USE_ATK
        SLIC_ERROR("Incorrect number of regionLayout entries specified in InternalMeshGenerator::ReadXML()");
#endif
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
    for( int block=1 ; block<m_nElems[dir].size() ; ++block )
    {
      m_firstElemIndexForBlock[dir][block] = m_lastElemIndexForBlock[dir][block-1] + 1;
      m_lastElemIndexForBlock[dir][block] = m_firstElemIndexForBlock[dir][block] + m_nElems[dir][block]-1;
    }
  }

  m_fPerturb = 0.0;

//    m_fPerturb = hdn.GetAttributeOrDefault<realT>("perturbationFactor", 0.0);
//    m_randSeed = hdn.GetAttributeOrDefault<int>("perturbationSeed",
// time(NULL));
//    srand(m_randSeed);
//
//    m_mapToRadial = hdn.GetAttributeOrDefault<int>("mapToRadial", 0);
//
//    m_skewAngle = hdn.GetAttributeOrDefault<realT>("skewAngle", 0.0);
//    m_skewAngle *= 3.14159265/180;
//    R1Tensor zeroVector;
//    zeroVector *= 0.0;
//    m_skewCenter = hdn.GetAttributeOrDefault<R1Tensor>("skewCenter",
// zeroVector);
//
//
//    // Mesh deformation
//    m_delayMeshDeformation =
// hdn.GetAttributeOrDefault<int>("delayMeshDeformation", 0);
//    m_meshDx = hdn.GetAttributeString("dxTable");
//    m_meshDy = hdn.GetAttributeString("dyTable");
//    m_meshDz = hdn.GetAttributeString("dzTable");

}



void InternalMeshGenerator::CreateChild( string const & childKey, string const & childName )
{
}


/**
 * @author settgast, fu, sherman
 * @param partition
 * @param domain
 */
void InternalMeshGenerator::GenerateMesh( dataRepository::ManagedGroup * const domain )
{
  // This cannot find groupkeys:
  // ManagedGroup * const meshBodies = domain->GetGroup(domain->groupKeys.meshBodies);
  ManagedGroup * const meshBodies = domain->GetGroup(std::string("MeshBodies"));
  MeshBody * const meshBody = meshBodies->RegisterGroup<MeshBody>( this->getName() );
  MeshLevel * const meshLevel0 = meshBody->RegisterGroup<MeshLevel>(std::string("Level0"));

  // special case
  //  bool isRadialWithOneThetaPartition = (m_mapToRadial > 0) &&
  // (partition.GetPartitions()[1]==1);

  NodeManager * nodeManager = meshLevel0->getNodeManager();

  // Make sure that the node manager fields are initialized
  nodeManager->SetDocumentationNodes();

  CellBlockManager * elementManager = domain->GetGroup<CellBlockManager>( keys::cellManager );
  ManagedGroup * nodeSets = nodeManager->GetGroup( std::string( "Sets" ) );

  PartitionBase & partition = domain->getReference<PartitionBase>(keys::partitionManager);

  bool isRadialWithOneThetaPartition = false;


  // This should probably handled elsewhere:
  for( auto & cellBlockName : m_regionNames )
  {
    CellBlock * cellBlock = elementManager->GetGroup(keys::cellBlocks)->RegisterGroup<CellBlock>(cellBlockName);
    cellBlock->SetDocumentationNodes();
    cellBlock->RegisterDocumentationNodes();
    cellBlock->ReadXML_PostProcess();
  }


  localIndex_set & xnegNodes = nodeSets->RegisterViewWrapper<localIndex_set>( std::string("xneg") )->reference();
  localIndex_set & xposNodes = nodeSets->RegisterViewWrapper<localIndex_set>( std::string("xpos") )->reference();
  localIndex_set & ynegNodes = nodeSets->RegisterViewWrapper<localIndex_set>( std::string("yneg") )->reference();
  localIndex_set & yposNodes = nodeSets->RegisterViewWrapper<localIndex_set>( std::string("ypos") )->reference();
  localIndex_set & znegNodes = nodeSets->RegisterViewWrapper<localIndex_set>( std::string("zneg") )->reference();
  localIndex_set & zposNodes = nodeSets->RegisterViewWrapper<localIndex_set>( std::string("zpos") )->reference();
  localIndex_set & allNodes  = nodeSets->RegisterViewWrapper<localIndex_set>( std::string("all") )->reference();


  // partition based on even spacing to get load balance
  // Partition geometrical boundaries will be corrected in the end.
  {
    m_min[0] = m_vertices[0].front();
    m_min[1] = m_vertices[1].front();
    m_min[2] = m_vertices[2].front();

    m_max[0] = m_vertices[0].back();
    m_max[1] = m_vertices[1].back();
    m_max[2] = m_vertices[2].back();

    R1Tensor temp1( m_min );
    R1Tensor temp2( m_max );
    partition.setSizes( temp1, temp2 );
  }

  // find elemCenters for even uniform element sizes
  array1d<array1d<real64> > elemCenterCoords( 3 );
  for( int i = 0 ; i < 3 ; ++i )
  {
    m_numElemsTotal[i] = 0;
    for( int block = 0 ; block < m_nElems[i].size() ; ++block )
    {
      m_numElemsTotal[i] += m_nElems[i][block];
    }

    elemCenterCoords[i].resize( m_numElemsTotal[i] );
    array1d<real64> elemCenterCoordsLocal( m_numElemsTotal[i] );
    for( int k = 0 ; k < m_numElemsTotal[i] ; ++k )
    {
      elemCenterCoordsLocal[k] = m_min[i] + ( m_max[i] - m_min[i] ) * ( k + 0.5 ) / m_numElemsTotal[i];
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
  int firstElemIndexInPartition[3] =
  { -1, -1, -1 };
  int lastElemIndexInPartition[3] =
  { -2, -2, -2 };

  for( int i = 0 ; i < 3 ; ++i )
  {
    //    firstElemIndexInPartition[i] = -1;
    //    lastElemIndexInPartition[i] = -2;
    for( int k = 0 ; k < m_numElemsTotal[i] ; ++k )
    {
      if( partition.IsCoordInPartition( elemCenterCoords[i][k], i ) )
      {
        firstElemIndexInPartition[i] = k;
        break;
      }
    }

    if( firstElemIndexInPartition[i] > -1 )
    {
      for( int k = firstElemIndexInPartition[i] ; k < m_numElemsTotal[i] ; ++k )
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

  std::map<std::string, int> numElemsInRegions;
  std::map<std::string, std::string> elemTypeInRegions;

  integer_array firstElemIndexForBlockInPartition[3];
  integer_array lastElemIndexForBlockInPartition[3];

  for( int dir = 0 ; dir < 3 ; ++dir )
  {
    firstElemIndexForBlockInPartition[dir] = m_firstElemIndexForBlock[dir];
    lastElemIndexForBlockInPartition[dir] = m_lastElemIndexForBlock[dir];

    for( int block = 0 ; block < m_nElems[dir].size() ; ++block )
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
  string_array::const_iterator iterRegion = m_regionNames.begin();
  for( int iblock = 0 ; iblock < m_nElems[0].size() ; ++iblock )
  {
    for( int jblock = 0 ; jblock < m_nElems[1].size() ; ++jblock )
    {
      for( int kblock = 0 ; kblock < m_nElems[2].size() ; ++kblock, ++iterRegion )
      {
        numElemsInRegions[*iterRegion] = 0;
        elemTypeInRegions[*iterRegion] = "";
      }
    }
  }

  iterRegion = m_regionNames.begin();
  {
    localIndex iR = 0;
    for( int iblock = 0 ; iblock < m_nElems[0].size() ; ++iblock )
    {
      for( int jblock = 0 ; jblock < m_nElems[1].size() ; ++jblock )
      {
        for( int kblock = 0 ; kblock < m_nElems[2].size() ; ++kblock, ++iterRegion, ++iR )
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
  localIndex numNodesInDir[3] =
  { 1, 1, 1 };

  for( int i = 0 ; i < m_dim ; ++i )
  {
    //    numElemsInDir[i] = lastElemIndexInPartition[i] -
    // firstElemIndexInPartition[i] + 1;
    numNodesInDir[i] = lastElemIndexInPartition[i] - firstElemIndexInPartition[i] + 2;
    if( isRadialWithOneThetaPartition && i == 1 )
    {
      numNodesInDir[1] -= 1;
    }
    numNodes *= numNodesInDir[i];
  }

  nodeManager->resize( numNodes );
  view_rtype<r1_array> X = nodeManager->getData<r1_array>( keys::referencePositionString );

  {
    localIndex localNodeIndex = 0;
    for( int i = 0 ; i < numNodesInDir[0] ; ++i )
    {
      for( int j = 0 ; j < numNodesInDir[1] ; ++j )
      {
        for( int k = 0 ; k < numNodesInDir[2] ; ++k )
        {
          int index[3] =
          { i, j, k };
          for( int a = 0 ; a < m_dim ; ++a )
          {
            index[a] += firstElemIndexInPartition[a];
          }

          X[localNodeIndex] = NodePosition( index, m_trianglePattern );

          // alter global node map for radial mesh
          if( m_mapToRadial > 0 )
          {
            if( isEqual( X[localNodeIndex][1], m_max[1], 1e-10 ) )
            {
              index[1] = 0;
            }
          }

          nodeManager->m_localToGlobalMap[localNodeIndex] = NodeGlobalIndex( index );

          // cartesian-specific nodesets
          if( m_mapToRadial == 0 )
          {
            if( isEqual( X[localNodeIndex][0], m_min[0], 1e-10 ) )
            {
              xnegNodes.insert( localNodeIndex );
            }
            if( isEqual( X[localNodeIndex][0], m_max[0], 1e-10 ) )
            {
              xposNodes.insert( localNodeIndex );
            }
            if( isEqual( X[localNodeIndex][1], m_min[1], 1e-10 ) )
            {
              ynegNodes.insert( localNodeIndex );
            }
            if( isEqual( X[localNodeIndex][1], m_max[1], 1e-10 ) )
            {
              yposNodes.insert( localNodeIndex );
            }
          }
          else
          {
            // radial-specific nodesets
            if( isEqual( X[localNodeIndex][0], m_min[0], 1e-10 ) )
            {
              xnegNodes.insert( localNodeIndex );
            }
            if( isEqual( X[localNodeIndex][0], m_max[0], 1e-10 ) )
            {
              xposNodes.insert( localNodeIndex );
            }
          }

          // general nodesets
          if( isEqual( X[localNodeIndex][2], m_min[2], 1e-10 ) )
          {
            znegNodes.insert( localNodeIndex );
          }
          if( isEqual( X[localNodeIndex][2], m_max[2], 1e-10 ) )
          {
            zposNodes.insert( localNodeIndex );
          }
          allNodes.insert( localNodeIndex );

          ++localNodeIndex;

        }
      }
    }
  }

  {
    integer_array numElements;
    string_array elementRegionNames;
    string_array elementTypes;
    std::map<std::string, localIndex> localElemIndexInRegion;

    for( std::map<std::string, int>::iterator iterNumElemsInRegion = numElemsInRegions.begin() ;
         iterNumElemsInRegion != numElemsInRegions.end() ; ++iterNumElemsInRegion )
    {

      numElements.push_back( iterNumElemsInRegion->second );
      elementRegionNames.push_back( iterNumElemsInRegion->first );
      elementTypes.push_back( elemTypeInRegions[iterNumElemsInRegion->first] );

      localElemIndexInRegion[iterNumElemsInRegion->first] = 0;
    }

    elementManager->resize( numElements, elementRegionNames, elementTypes );

    // assign global numbers to elements
    iterRegion = m_regionNames.begin();
    set<std::string> processedRegionNames;
    localIndex iR = 0;

    for( int iblock = 0 ; iblock < m_nElems[0].size() ; ++iblock )
    {
      for( int jblock = 0 ; jblock < m_nElems[1].size() ; ++jblock )
      {
        for( int kblock = 0 ; kblock < m_nElems[2].size() ; ++kblock, ++iterRegion, ++iR )
        {
//          ElementRegionT& elemRegion =
// domain->m_feElementManager->m_ElementRegions[*iterRegion];

          CellBlock * elemRegion =  elementManager->GetRegion(*iterRegion);
          int const numNodesPerElem = integer_conversion<int>(elemRegion->numNodesPerElement());
          FixedOneToManyRelation & elemsToNodes = elemRegion->nodeList();

          int numElemsInDirForRegion[3] =
          { lastElemIndexForBlockInPartition[0][iblock] - firstElemIndexForBlockInPartition[0][iblock] + 1,
            lastElemIndexForBlockInPartition[1][jblock] - firstElemIndexForBlockInPartition[1][jblock] + 1,
            lastElemIndexForBlockInPartition[2][kblock] - firstElemIndexForBlockInPartition[2][kblock] + 1 };

          for( int i = 0 ; i < numElemsInDirForRegion[0] ; ++i )
          {
            for( int j = 0 ; j < numElemsInDirForRegion[1] ; ++j )
            {
              for( int k = 0 ; k < numElemsInDirForRegion[2] ; ++k )
              {
                int index[3] =
                { i + firstElemIndexForBlockInPartition[0][iblock],
                  j + firstElemIndexForBlockInPartition[1][jblock],
                  k + firstElemIndexForBlockInPartition[2][kblock] };

                const localIndex firstNodeIndex = numNodesInDir[1] * numNodesInDir[2] * ( index[0] - firstElemIndexInPartition[0] )
                                                  + numNodesInDir[2] * ( index[1] - firstElemIndexInPartition[1] )
                                                  + ( index[2] - firstElemIndexInPartition[2] );
                localIndex nodeOfBox[8];

                if( m_elementType[iR] == "CPE4" || m_elementType[iR] == "STRI" )
                {

                  nodeOfBox[0] = firstNodeIndex;
                  nodeOfBox[1] = numNodesInDir[1] * numNodesInDir[2] + firstNodeIndex;
                  nodeOfBox[2] = numNodesInDir[1] * numNodesInDir[2] + numNodesInDir[2] + firstNodeIndex;
                  nodeOfBox[3] = numNodesInDir[2] + firstNodeIndex;

                }
                else
                {
                  nodeOfBox[0] = firstNodeIndex;
                  nodeOfBox[1] = numNodesInDir[1] * numNodesInDir[2] + firstNodeIndex;
                  nodeOfBox[2] = numNodesInDir[1] * numNodesInDir[2] + numNodesInDir[2] + firstNodeIndex;
                  nodeOfBox[3] = numNodesInDir[2] + firstNodeIndex;

                  nodeOfBox[4] = firstNodeIndex + 1;
                  nodeOfBox[5] = numNodesInDir[1] * numNodesInDir[2] + firstNodeIndex + 1;
                  nodeOfBox[6] = numNodesInDir[1] * numNodesInDir[2] + numNodesInDir[2] + firstNodeIndex + 1;
                  nodeOfBox[7] = numNodesInDir[2] + firstNodeIndex + 1;

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
                  if( j == numElemsInDirForRegion[1] - 1 && jblock == m_nElems[1].size() - 1 )
                  { // last set of elements
                    index[1] = -1;
                    const localIndex firstNodeIndexR = numNodesInDir[1] * numNodesInDir[2] * ( index[0] - firstElemIndexInPartition[0] )
                                                       + numNodesInDir[2] * ( index[1] - firstElemIndexInPartition[1] )
                                                       + ( index[2] - firstElemIndexInPartition[2] );
                    nodeOfBox[2] = numNodesInDir[1] * numNodesInDir[2] + numNodesInDir[2] + firstNodeIndexR;
                    nodeOfBox[3] = numNodesInDir[2] + firstNodeIndexR;
                    nodeOfBox[6] = numNodesInDir[1] * numNodesInDir[2] + numNodesInDir[2] + firstNodeIndexR + 1;
                    nodeOfBox[7] = numNodesInDir[2] + firstNodeIndexR + 1;
                  }
                }

                for( int iEle = 0 ; iEle < m_numElePerBox[iR] ; ++iEle )
                {
                  localIndex& localElemIndex = localElemIndexInRegion[*iterRegion];
                  elemRegion->m_localToGlobalMap[localElemIndex] = ElemGlobalIndex( index ) * m_numElePerBox[iR] + iEle;

                  integer_array nodeIDInBox( numNodesPerElem );

                  GetElemToNodesRelationInBox( m_elementType[iR], index, iEle, nodeIDInBox.data(),
                                               numNodesPerElem );

                  for( localIndex iN = 0 ; iN < numNodesPerElem ; ++iN )
                  {
// #ifdef USE_ATK
//                    SLIC_ERROR("not implemented");
// #endif
                    elemsToNodes[localElemIndex][iN] = nodeOfBox[nodeIDInBox[iN]];
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

#if 0
  // Correct partition geometrical boundary.
  {
    R1Tensor pMin, pMax;
    partition.getPartitionGeometricalBoundary( pMin, pMax );
    for( int i = 0 ; i < m_dim ; ++i )
    {
      realT xMinByNumElems = ( pMin[i] - m_min[i] ) / ( m_max[i] - m_min[i] ) * m_numElemsTotal[i];
      realT xMaxByNumElems = ( pMax[i] - m_min[i] ) / ( m_max[i] - m_min[i] ) * m_numElemsTotal[i];

      localIndex iBlockMin( 0 ), iBlockMax( 0 );
      while( ( xMinByNumElems < m_firstElemIndexForBlock[i][iBlockMin] * 1.0 || xMinByNumElems > m_lastElemIndexForBlock[i][iBlockMin] * 1.0 + 1.0 )
             && iBlockMin < m_nElems[i].size() - 1 )
      {
        ++iBlockMin;
      }
      while( ( xMaxByNumElems < m_firstElemIndexForBlock[i][iBlockMax] * 1.0 || xMaxByNumElems > m_lastElemIndexForBlock[i][iBlockMax] * 1.0 + 1.0 )
             && iBlockMax < m_nElems[i].size() - 1 )
      {
        ++iBlockMax;
      }
      //sanity check

      if( xMinByNumElems < m_firstElemIndexForBlock[i][iBlockMin] * 1.0 ||
          xMinByNumElems > m_lastElemIndexForBlock[i][iBlockMin] * 1.0 + 1.0 ||
          iBlockMin >= m_nElems[i].size() )
      {
        std::cout
          << "WARNING: Had some trouble in correcting partition geometric boundaries.  If this affects the results, contact a developer"
          << std::endl;
      }
      if( xMaxByNumElems < m_firstElemIndexForBlock[i][iBlockMax] * 1.0 ||
          xMaxByNumElems > m_lastElemIndexForBlock[i][iBlockMax] * 1.0 + 1.0 ||
          iBlockMax >= m_nElems[i].size() )
      {
        std::cout
          << "WARNING: Had some trouble in correcting partition geometric boundaries.  If this affects the results, contact a developer"
          << std::endl;
      }

      pMin[i] = m_vertices[i][iBlockMin] + ( xMinByNumElems - m_firstElemIndexForBlock[i][iBlockMin] ) *
                ( m_vertices[i][iBlockMin + 1] - m_vertices[i][iBlockMin] ) / m_nElems[i][iBlockMin];
      pMax[i] = m_vertices[i][iBlockMax] + ( xMaxByNumElems - m_firstElemIndexForBlock[i][iBlockMax] ) *
                ( m_vertices[i][iBlockMax + 1] - m_vertices[i][iBlockMax] ) / m_nElems[i][iBlockMax];
    }

    partition.SetPartitionGeometricalBoundary( pMin, pMax );

  }
#endif

  /*
     {// Move nodes in the extension layers.

     for (localIndex iN = 0; iN != nodeManager->DataLengths(); ++iN)
     {
     for (int i=0; i<m_dim; ++i)
     {
     if ( X[iN][i] < m_min[i])
     {
     int eLayer = (int) ((m_min[i] -X[iN][i]) / ((m_max[i] - m_min[i]) /
        m_numElems[i]) + 0.5);
     X[iN][i] = m_min[i] - ((m_max[i] - m_min[i]) / m_numElems[i]) *
        m_commonRatioMin[i] * (1- pow(m_commonRatioMin[i], eLayer)) / (1 -
        m_commonRatioMin[i]);
     }
     else if (X[iN][i] > m_max[i])
     {
     int eLayer = (int) ((X[iN][i] - m_max[i] ) / ((m_max[i] - m_min[i]) /
        m_numElems[i]) + 0.5);
     X[iN][i] = m_max[i] + ((m_max[i] - m_min[i]) / m_numElems[i]) *
        m_commonRatioMax[i] * (1- pow(m_commonRatioMax[i], eLayer)) / (1 -
        m_commonRatioMax[i]);
     }

     }
     }
     }
   */

  // Node perturbation
  if( m_fPerturb > 0 )
  {
    for( localIndex iN = 0 ; iN != nodeManager->size() ; ++iN )
    {

      for( int i = 0 ; i < m_dim ; ++i )
      {
        if( X[iN][i] > m_min[i] && X[iN][i] < m_max[i] )
        {
          srand( integer_conversion<int>(nodeManager->m_localToGlobalMap[iN]) + m_randSeed + i ); // This
                                                                                                  // ensures
                                                                                                  // that
                                                                                                  // the
                                                                                                  // perturbation
                                                                                                  // pattern
                                                                                                  // is
                                                                                                  // unaffected
                                                                                                  // by
                                                                                                  // domain
          X[iN][i] += ( ( m_max[i] - m_min[i] ) / m_numElemsTotal[i] ) * ( ( rand() * 1.0 ) / RAND_MAX - 0.5 ) * 2 * m_fPerturb;
        }
      }
    }
  }

  if( std::fabs( m_skewAngle ) > 0.0 )
  {
    for( localIndex iN = 0 ; iN != nodeManager->size() ; ++iN )
    {
      X[iN][0] -= ( X[iN][1] - m_skewCenter[1] ) * std::tan( m_skewAngle );
    }
  }

  if( m_mapToRadial > 0 )
  {
    // Map to radial mesh
    for( localIndex iN = 0 ; iN != nodeManager->size() ; ++iN )
    {
      m_meshTheta = X[iN][1] * M_PI / 180.0;
      m_meshAxis = static_cast<int>(round( m_meshTheta * 2.0 / M_PI ));
      m_meshPhi = fabs( m_meshTheta - m_meshAxis * M_PI / 2.0 );
      m_meshRout = m_max[0] / cos( m_meshPhi );

      if( m_mapToRadial > 1 )
      {
        m_meshRact = ( ( m_meshRout - m_min[0] ) / ( m_max[0] - m_min[0] ) ) * ( X[iN][0] - m_min[0] ) + m_min[0];
      }
      else
      {
        m_meshRact = X[iN][0];
      }

      X[iN][0] = m_meshRact * cos( m_meshTheta );
      X[iN][1] = m_meshRact * sin( m_meshTheta );

      // add mapped values to nodesets
      if( m_mapToRadial > 1 )
      {
        if( isEqual( X[iN][0], -1 * m_max[0], 1e-6 ) )
        {
          xnegNodes.insert( iN );
        }
        if( isEqual( X[iN][0], m_max[0], 1e-6 ) )
        {
          xposNodes.insert( iN );
        }
        if( isEqual( X[iN][1], -1 * m_max[0], 1e-6 ) )
        {
          ynegNodes.insert( iN );
        }
        if( isEqual( X[iN][1], m_max[0], 1e-6 ) )
        {
          yposNodes.insert( iN );
        }
      }
    }
  }

  if( m_delayMeshDeformation == 0 )
  {
    RemapMesh( domain );
  }
}

/**
 * @author fu
 * @param elementType
 * @param index
 * @param iEle
 * @param nodeIDInBox
 * @param node_size
 */
void InternalMeshGenerator::GetElemToNodesRelationInBox( const std::string& elementType,
                                                         const int index[],
                                                         const int& iEle,
                                                         int nodeIDInBox[],
                                                         const int node_size )

{
  if( elementType == "C3D8" )
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
  else if( ( elementType == "C3D6" ) && ( m_trianglePattern == 0 ) )
  {
    if( ( index[0] + index[1] ) % 2 == 1 )
    {
      if( iEle % 2 == 0 )
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
      if( iEle % 2 == 0 )
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
  else if( ( elementType == "C3D6" ) && ( m_trianglePattern == 1 ) )
  {
    if( index[1] % 2 == 0 )
    {
      if( iEle % 2 == 0 )
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
      if( iEle % 2 == 0 )
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
  else if( elementType == "CPE4" )
  {
    nodeIDInBox[0] = 0;
    nodeIDInBox[1] = 1;
    nodeIDInBox[2] = 3;
    nodeIDInBox[3] = 2;
  }
  else if( ( elementType == "STRI" ) && ( m_trianglePattern == 0 ) )
  {
    if( ( index[0] + index[1] ) % 2 == 1 )
    {
      if( iEle % 2 == 0 )
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
      if( iEle % 2 == 0 )
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
  else if( ( elementType == "STRI" ) && ( m_trianglePattern == 1 ) )
  {
    if( index[1] % 2 == 0 )
    {
      if( iEle % 2 == 0 )
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
      if( iEle % 2 == 0 )
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
  else if( elementType == "C3D4" )
  {
    int mapBoxTet[8][6][4] =
    {
      {
        { 0, 3, 7, 6 },
        { 0, 7, 4, 6 },
        { 0, 5, 1, 6 },
        { 0, 4, 5, 6 },
        { 0, 1, 2, 6 },
        { 0, 2, 3, 6 }
      },
      {
        { 1, 5, 6, 7 },
        { 1, 6, 2, 7 },
        { 0, 4, 1, 7 },
        { 1, 4, 5, 7 },
        { 0, 1, 3, 7 },
        { 1, 2, 3, 7 }
      },
      {
        { 0, 3, 4, 2 },
        { 3, 7, 4, 2 },
        { 0, 4, 1, 2 },
        { 1, 4, 5, 2 },
        { 4, 7, 6, 2 },
        { 4, 6, 5, 2 }
      },
      {
        { 1, 5, 2, 3 },
        { 2, 5, 6, 3 },
        { 0, 5, 1, 3 },
        { 0, 4, 5, 3 },
        { 4, 7, 5, 3 },
        { 5, 7, 6, 3 }
      },
      {
        { 0, 3, 4, 5 },
        { 3, 7, 4, 5 },
        { 3, 6, 7, 5 },
        { 3, 2, 6, 5 },
        { 3, 0, 1, 5 },
        { 1, 2, 3, 5 }
      },
      {
        { 1, 5, 2, 4 },
        { 2, 5, 6, 4 },
        { 2, 6, 7, 4 },
        { 2, 7, 3, 4 },
        { 0, 2, 3, 4 },
        { 0, 1, 2, 4 }
      },
      {
        { 0, 7, 4, 1 },
        { 0, 3, 7, 1 },
        { 2, 7, 3, 1 },
        { 2, 6, 7, 1 },
        { 4, 7, 5, 1 },
        { 5, 7, 6, 1 }
      },
      {
        { 1, 5, 6, 0 },
        { 1, 6, 2, 0 },
        { 2, 6, 3, 0 },
        { 3, 6, 7, 0 },
        { 4, 6, 5, 0 },
        { 4, 7, 6, 0 }
      }
    };

    int mapBoxType[2][2][2];
    mapBoxType[0][0][0] = 0;
    mapBoxType[1][0][0] = 1;
    mapBoxType[0][0][1] = 2;
    mapBoxType[1][0][1] = 3;
    mapBoxType[0][1][0] = 4;
    mapBoxType[1][1][0] = 5;
    mapBoxType[0][1][1] = 6;
    mapBoxType[1][1][1] = 7;

    int boxType = mapBoxType[index[0] % 2][index[1] % 2][index[2] % 2];
    for( int i = 0 ; i < node_size ; ++i )
    {
      nodeIDInBox[i] = mapBoxTet[boxType][iEle][i];
    }

  }
}

void InternalMeshGenerator::RemapMesh( dataRepository::ManagedGroup * const domain )
{
  //  // Node mapping
  //  if (!m_meshDx.empty())
  //  {
  //    const Table3D* tableDx =
  // stlMapLookupPointer(TableManager::Instance().Tables<3>(), m_meshDx);
  //
  //    for (localIndex iN=0; iN!=nodeManager->DataLengths(); ++iN)
  //    {
  //      realT dx=tableDx->Lookup(X[iN]);
  //      X[iN][0] += dx;
  //    }
  //  }
  //
  //  if (!m_meshDy.empty())
  //  {
  //    const Table3D* tableDy =
  // stlMapLookupPointer(TableManager::Instance().Tables<3>(), m_meshDy);
  //
  //    for (localIndex iN=0; iN!=nodeManager->DataLengths(); ++iN)
  //    {
  //      realT dy=tableDy->Lookup(X[iN]);
  //      X[iN][1] += dy;
  //    }
  //  }
  //
  //  if (!m_meshDz.empty())
  //  {
  //    const Table3D* tableDz =
  // stlMapLookupPointer(TableManager::Instance().Tables<3>(), m_meshDz);
  //
  //    for (localIndex iN=0; iN!=nodeManager->DataLengths(); ++iN)
  //    {
  //      realT dz=tableDz->Lookup(X[iN]);
  //      X[iN][2] += dz;
  //    }
  //  }

}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, InternalMeshGenerator, std::string const &, ManagedGroup * const )
}
