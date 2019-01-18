/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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
 * PAMELAMeshGenerator.cpp
 *
 *  Created on: Oct 08, 2018
 *      Author: Antoine Mazuyer
 */

#include "PAMELAMeshGenerator.hpp"

#include "managers/DomainPartition.hpp"

#include "codingUtilities/StringUtilities.hpp"
#include <math.h>

#ifdef USE_ATK
#include "slic/slic.hpp"
#endif

#include "Mesh/MeshFactory.hpp"

#include "MeshDataWriters/MeshParts.hpp"
#include "MeshDataWriters/VTKWriter.hpp"

#include "MPI_Communications/PartitionBase.hpp"
#include "MPI_Communications/SpatialPartition.hpp"

#include "mesh/MeshBody.hpp"

namespace geosx
{
using namespace dataRepository;

PAMELAMeshGenerator::PAMELAMeshGenerator( string const & name, ManagedGroup * const parent ):
  MeshGeneratorBase( name, parent )
{
  this->RegisterViewWrapper<string>(keys::filePath)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("path to the mesh file");
}

PAMELAMeshGenerator::~PAMELAMeshGenerator()
{}

void PAMELAMeshGenerator::GenerateElementRegions( DomainPartition& domain )
{}

void PAMELAMeshGenerator::PostProcessInput()
{
  m_pamelaMesh =
    std::unique_ptr< PAMELA::Mesh >
      ( PAMELA::MeshFactory::makeMesh( this->getReference<string>( keys::filePath )));
  m_pamelaMesh->CreateFacesFromCells();
  m_pamelaMesh->PerformPolyhedronPartitioning( PAMELA::ELEMENTS::FAMILY::POLYGON,
                                               PAMELA::ELEMENTS::FAMILY::POLYGON );
  m_pamelaMesh->CreateLineGroupWithAdjacency(
    "TopologicalC2C",
    m_pamelaMesh->getAdjacencySet()->get_TopologicalAdjacency( PAMELA::ELEMENTS::FAMILY::POLYHEDRON, PAMELA::ELEMENTS::FAMILY::POLYHEDRON,
                                                               PAMELA::ELEMENTS::FAMILY::POLYGON ));

}

void PAMELAMeshGenerator::RemapMesh( dataRepository::ManagedGroup * const domain )
{
  return;
}

ManagedGroup * PAMELAMeshGenerator::CreateChild( string const & childKey, string const & childName )
{
  return nullptr;
}

void PAMELAMeshGenerator::GenerateMesh( dataRepository::ManagedGroup * const domain )
{
  int nbRanks = 1;
  int thisRank;
#ifdef GEOSX_USE_MPI
  MPI_Comm_size( MPI_COMM_GEOSX, &nbRanks );
  MPI_Comm_rank( MPI_COMM_GEOSX, &thisRank );
#endif
  ManagedGroup * const meshBodies = domain->GetGroup( std::string( "MeshBodies" ));
  MeshBody * const meshBody = meshBodies->RegisterGroup<MeshBody>( this->getName() );

  //TODO for the moment we only consider on mesh level "Level0"
  MeshLevel * const meshLevel0 = meshBody->RegisterGroup<MeshLevel>( std::string( "Level0" ));
  NodeManager * nodeManager = meshLevel0->getNodeManager();
  CellBlockManager * cellBlockManager = domain->GetGroup<CellBlockManager>( keys::cellManager );


  // Use the PartMap of PAMELA to get the mesh
  auto polyhedronPartMap = std::get<0>( PAMELA::getPolyhedronPartMap( m_pamelaMesh.get()));

  // Vertices are written first
  r1_array const & X = nodeManager->referencePosition();
  nodeManager->resize(m_pamelaMesh->get_PointCollection()->size_all());
  localIndex nbVertices = m_pamelaMesh->get_PointCollection()->size_all();
  for( auto verticesIterator : *m_pamelaMesh->get_PointCollection()) {
    localIndex vertexLocalIndex = verticesIterator->get_localIndex();
    globalIndex vertexGlobalIndex = verticesIterator->get_globalIndex();
    real64 * const pointData = X[verticesIterator->get_localIndex()].Data();
    pointData[0] = verticesIterator->get_coordinates().x;
    pointData[1] = verticesIterator->get_coordinates().y;
    pointData[2] = verticesIterator->get_coordinates().z;
    nodeManager->m_localToGlobalMap[vertexLocalIndex] = vertexGlobalIndex;
  }
  
  // First loop which iterate on the regions
  for( auto regionItr = polyhedronPartMap.begin() ; regionItr != polyhedronPartMap.end() ; ++regionItr )
  {
    auto regionPtr = regionItr->second;
    auto regionIndex = regionPtr->Index;
    auto regionIndexStr = std::to_string(regionIndex);

    // Iterate on cell types
    for( auto cellBlockIterator = regionPtr->SubParts.begin() ;
        cellBlockIterator != regionPtr->SubParts.end() ; cellBlockIterator++ )
    {
      auto cellBlockPAMELA = cellBlockIterator->second;
      auto cellBlockType = cellBlockPAMELA->ElementType;
      auto cellBlockName = ElementToLabel.at( cellBlockType );
      CellBlock * cellBlock = nullptr;
      if( cellBlockName == "HEX" )
      {
        cellBlock =
          cellBlockManager->GetGroup( keys::cellBlocks )->RegisterGroup<CellBlock>( regionIndexStr + "_" + cellBlockName);
        cellBlock -> SetElementType("C3D8");
        auto nbCells = cellBlockPAMELA->SubCollection.size_owned();
        auto & cellToVertex = cellBlock->nodeList();
        cellBlock->resize( nbCells );
        cellToVertex.resize( nbCells, 8 );

        // Iterate on cells
        for( auto cellItr = cellBlockPAMELA->SubCollection.begin_owned() ;
            cellItr != cellBlockPAMELA->SubCollection.end_owned() ;
            cellItr++ )
        {
          localIndex cellLocalIndex = (*cellItr)->get_localIndex();
          globalIndex cellGlobalIndex = (*cellItr)->get_globalIndex();
          auto cornerList = (*cellItr)->get_vertexList();

          cellToVertex[cellLocalIndex][0] =
            cornerList[0]->get_localIndex();
          cellToVertex[cellLocalIndex][1] =
            cornerList[1]->get_localIndex();
          cellToVertex[cellLocalIndex][2] =
            cornerList[3]->get_localIndex();
          cellToVertex[cellLocalIndex][3] =
            cornerList[2]->get_localIndex();
          cellToVertex[cellLocalIndex][4] =
            cornerList[4]->get_localIndex();
          cellToVertex[cellLocalIndex][5] =
            cornerList[5]->get_localIndex();
          cellToVertex[cellLocalIndex][6] =
            cornerList[7]->get_localIndex();
          cellToVertex[cellLocalIndex][7] =
            cornerList[6]->get_localIndex();

          cellBlock->m_localToGlobalMap[cellLocalIndex] = cellGlobalIndex;
        }
      }
      else if( cellBlockName == "TETRA" )
      {
        cellBlock =
          cellBlockManager->GetGroup( keys::cellBlocks )->RegisterGroup<CellBlock>( regionIndexStr + "_" + cellBlockName);
        cellBlock -> SetElementType("C3D4");
        auto nbCells = cellBlockPAMELA->SubCollection.size_owned();
        auto & cellToVertex = cellBlock->nodeList();
        cellBlock->resize( nbCells );
        cellToVertex.resize( nbCells, 4 );

        // Iterate on cells
        for( auto cellItr = cellBlockPAMELA->SubCollection.begin_owned() ;
            cellItr != cellBlockPAMELA->SubCollection.end_owned() ;
            cellItr++ )
        {
          localIndex cellLocalIndex = (*cellItr)->get_localIndex();
          globalIndex cellGlobalIndex = (*cellItr)->get_globalIndex();
          auto cornerList = (*cellItr)->get_vertexList();

          cellToVertex[cellLocalIndex][0] =
            cornerList[0]->get_localIndex();
          cellToVertex[cellLocalIndex][1] =
            cornerList[1]->get_localIndex();
          cellToVertex[cellLocalIndex][2] =
            cornerList[2]->get_localIndex();
          cellToVertex[cellLocalIndex][3] =
            cornerList[3]->get_localIndex();

          cellBlock->m_localToGlobalMap[cellLocalIndex] = cellGlobalIndex;
        }
      }
      else if( cellBlockName == "WEDGE" )
      {
        cellBlock =
          cellBlockManager->GetGroup( keys::cellBlocks )->RegisterGroup<CellBlock>( regionIndexStr + "_" + cellBlockName);
        cellBlock -> SetElementType("C3D6");
        auto nbCells = cellBlockPAMELA->SubCollection.size_owned();
        auto & cellToVertex = cellBlock->nodeList();
        cellBlock->resize( nbCells );
        cellToVertex.resize( nbCells, 6 );

        // Iterate on cells
        for( auto cellItr = cellBlockPAMELA->SubCollection.begin_owned() ;
            cellItr != cellBlockPAMELA->SubCollection.end_owned() ;
            cellItr++ )
        {
          localIndex cellLocalIndex = (*cellItr)->get_localIndex();
          globalIndex cellGlobalIndex = (*cellItr)->get_globalIndex();
          auto cornerList = (*cellItr)->get_vertexList();

          cellToVertex[cellLocalIndex][0] =
            cornerList[0]->get_localIndex();
          cellToVertex[cellLocalIndex][1] =
            cornerList[1]->get_localIndex();
          cellToVertex[cellLocalIndex][2] =
            cornerList[2]->get_localIndex();
          cellToVertex[cellLocalIndex][3] =
            cornerList[3]->get_localIndex();
          cellToVertex[cellLocalIndex][4] =
            cornerList[4]->get_localIndex();
          cellToVertex[cellLocalIndex][5] =
            cornerList[5]->get_localIndex();

          cellBlock->m_localToGlobalMap[cellLocalIndex] = cellGlobalIndex;
        }
      }
      else if( cellBlockName == "PYRAMID" )
      {
        cellBlock =
          cellBlockManager->GetGroup( keys::cellBlocks )->RegisterGroup<CellBlock>( regionIndexStr + "_" + cellBlockName);
        cellBlock -> SetElementType("C3D5");
        auto nbCells = cellBlockPAMELA->SubCollection.size_owned();
        auto & cellToVertex = cellBlock->nodeList();
        cellBlock->resize( nbCells );
        cellToVertex.resize( nbCells, 5 );

        // Iterate on cells
        for( auto cellItr = cellBlockPAMELA->SubCollection.begin_owned() ;
            cellItr != cellBlockPAMELA->SubCollection.end_owned() ;
            cellItr++ )
        {
          localIndex cellLocalIndex = (*cellItr)->get_localIndex();
          globalIndex cellGlobalIndex = (*cellItr)->get_globalIndex();
          auto cornerList = (*cellItr)->get_vertexList();

          cellToVertex[cellLocalIndex][0] =
            cornerList[0]->get_localIndex();
          cellToVertex[cellLocalIndex][1] =
            cornerList[1]->get_localIndex();
          cellToVertex[cellLocalIndex][2] =
            cornerList[2]->get_localIndex();
          cellToVertex[cellLocalIndex][3] =
            cornerList[3]->get_localIndex();
          cellToVertex[cellLocalIndex][4] =
            cornerList[4]->get_localIndex();

          cellBlock->m_localToGlobalMap[cellLocalIndex] = cellGlobalIndex;
        }
      }
    }
  }

}

void PAMELAMeshGenerator::GetElemToNodesRelationInBox( const std::string& elementType,
                                                       const int index[],
                                                       const int& iEle,
                                                       int nodeIDInBox[],
                                                       const int node_size )
{}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, PAMELAMeshGenerator, std::string const &, ManagedGroup * const )
}
