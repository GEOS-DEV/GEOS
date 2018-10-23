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
    // TODO : To be written --> Next PR !
}

PAMELAMeshGenerator::~PAMELAMeshGenerator()
{
    // TODO : To be written --> Next PR !
}

void PAMELAMeshGenerator::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName( "PAMELAMeshGenerator" );
  docNode->setSchemaType( "Node" );
  docNode->setShortDescription( "Generate a mesh using PAMELA library" );

  docNode->AllocateChildNode( keys::filePath,
                              keys::filePath,
                              -1,
                              "string",
                              "string",
                              "path to the mesh file",
                              "path to the mesh file",
                              "filePath",
                              "",
                              0,
                              1,
                              0 );
}

void PAMELAMeshGenerator::GenerateElementRegions( DomainPartition& domain )
{
    // TODO : To be written --> Next PR !
}

void PAMELAMeshGenerator::ReadXML_PostProcess()
{
    m_pamelaMesh = 
      std::unique_ptr< PAMELA::Mesh >
      (PAMELA::MeshFactory::makeMesh(this->getReference<string>(keys::filePath)));
}

void PAMELAMeshGenerator::RemapMesh(dataRepository::ManagedGroup * const domain)
{
  return;
}

void PAMELAMeshGenerator::CreateChild( string const & childKey, string const & childName )
{
  return;
}

void PAMELAMeshGenerator::GenerateMesh( dataRepository::ManagedGroup * const domain )
{
  int nranks = 1;
#ifdef GEOSX_USE_MPI
  MPI_Comm_size(MPI_COMM_GEOSX, &nranks);
#endif

  // We throw an error if GEOSX in launched with MPI. Multirank will be handle in another PR
  GEOS_ERROR_IF(nranks > 1, "PAMELA Mesh Generator not yet compatible with multidomains");

  ManagedGroup * const meshBodies = domain->GetGroup(std::string("MeshBodies"));
  MeshBody * const meshBody = meshBodies->RegisterGroup<MeshBody>( this->getName() );

  //TODO for the moment we only consider on mesh level "Level0"
  MeshLevel * const meshLevel0 = meshBody->RegisterGroup<MeshLevel>(std::string("Level0"));
  NodeManager * nodeManager = meshLevel0->getNodeManager();
  arrayView1d<R1Tensor> X = nodeManager->referencePosition();
  //CellBlockManager * elementManager = domain->GetGroup<CellBlockManager>( keys::cellManager );
  ElementRegionManager * elementManager = meshBody->getMeshLevel(0)->getElemManager();

  //TODO for the moment we only write the polyhedron and the associated vertices
  auto polyhedronCollection = m_pamelaMesh->get_PolyhedronCollection();

  // Use the PartMap of PAMELA to get the mesh
  auto polyhedronPartMap = std::get<0>(PAMELA::getPolyhedronPartMap(m_pamelaMesh.get()));
  
  // First loop which iterate on the regions
  for(auto regionItr = polyhedronPartMap.begin() ; regionItr != polyhedronPartMap.end();++regionItr) {
    auto regionPtr = regionItr->second;
    auto regionName = regionPtr->Label;

    // Check if the region is declared.
    GEOS_ERROR_IF(!elementManager->GetRegion(regionName), "Region " + regionName +
        "found in the mesh file but not declared in the XML command file");
    auto region = elementManager->GetRegion(regionName);

    // Iterate on vertices
    for(auto verticesIterator = regionPtr->Points.begin() ;
        verticesIterator != regionPtr->Points.end(); verticesIterator++) {
      real64 * const pointData = X[(*verticesIterator)->get_globalIndex()].Data();
      pointData[0] = (*verticesIterator)->get_coordinates().x;
      pointData[1] = (*verticesIterator)->get_coordinates().y;
      pointData[2] = (*verticesIterator)->get_coordinates().z;
    }

    // Iterate on cell types
    for(auto cellBlockIterator = regionPtr->SubParts.begin() ;
        cellBlockIterator != regionPtr->SubParts.end() ; cellBlockIterator++) {
      // Check if there is cell of this type
      if(cellBlockIterator->second->SubCollection.size_owned() > 0) {
        auto cellBlockPAMELA = cellBlockIterator->second;
        auto cellBlockType = cellBlockPAMELA->ElementType;
        auto cellBlockName = ElementToLabel.at(cellBlockType);
        CellBlock * cellBlock =
          elementManager->GetGroup(keys::cellBlocks)->RegisterGroup<CellBlock>(cellBlockName);

        if( cellBlockName == "HEX") {
          auto nbCells = cellBlockPAMELA->SubCollection.size_owned();
          auto & cellToVertex = cellBlock->nodeList();
          cellBlock->resize(nbCells);
          cellToVertex.resize(nbCells, 8 );
          
          // Iterate on cells
          for(auto cellItr = cellBlockPAMELA->SubCollection.begin_owned() ;
              cellItr != cellBlockPAMELA->SubCollection.end_owned() ;
              cellItr++ ) {
            auto cellIndex = (*cellItr)->get_globalIndex();
            auto cornerList = (*cellItr)->get_vertexList();
            cellToVertex[cellIndex][0] =
              regionPtr->GlobalToLocalPointMapping.at(cornerList[0]->get_globalIndex());
            cellToVertex[cellIndex][1] =
              regionPtr->GlobalToLocalPointMapping.at(cornerList[1]->get_globalIndex());
            cellToVertex[cellIndex][2] =
              regionPtr->GlobalToLocalPointMapping.at(cornerList[3]->get_globalIndex());
            cellToVertex[cellIndex][3] =
              regionPtr->GlobalToLocalPointMapping.at(cornerList[2]->get_globalIndex());
            cellToVertex[cellIndex][4] =
              regionPtr->GlobalToLocalPointMapping.at(cornerList[4]->get_globalIndex());
            cellToVertex[cellIndex][5] =
              regionPtr->GlobalToLocalPointMapping.at(cornerList[5]->get_globalIndex());
            cellToVertex[cellIndex][6] =
              regionPtr->GlobalToLocalPointMapping.at(cornerList[7]->get_globalIndex());
            cellToVertex[cellIndex][7] =
              regionPtr->GlobalToLocalPointMapping.at(cornerList[6]->get_globalIndex());
          }
        }
      }
    }
  }
  // Writting nodes in GEOSX data structure
  //for( int vertexIndex = 0 ; vertexIndex < polyhedronCollection-
}

void PAMELAMeshGenerator::GetElemToNodesRelationInBox( const std::string& elementType,
    const int index[],
    const int& iEle,
    int nodeIDInBox[],
    const int node_size )
{

}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, PAMELAMeshGenerator, std::string const &, ManagedGroup * const )
}
