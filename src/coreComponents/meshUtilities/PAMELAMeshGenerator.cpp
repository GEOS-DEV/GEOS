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

#include "Mesh/MeshFactory.hpp"

#include "MeshDataWriters/MeshParts.hpp"

#include "MPI_Communications/PartitionBase.hpp"
#include "MPI_Communications/SpatialPartition.hpp"

#include "mesh/MeshBody.hpp"

namespace geosx
{
using namespace dataRepository;

PAMELAMeshGenerator::PAMELAMeshGenerator( string const & name, ManagedGroup * const parent ):
  MeshGeneratorBase( name, parent )
{

  RegisterViewWrapper(viewKeyStruct::filePathString, &m_filePath, false)->
    setInputFlag(InputFlags::REQUIRED)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("path to the mesh file");
  RegisterViewWrapper(viewKeyStruct::fieldsToImportString, &m_fieldsToImport, false)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Fields to be imported from the external mesh file");
  RegisterViewWrapper(viewKeyStruct::fieldNamesInGEOSXString, &m_fieldNamesInGEOSX, false)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Name of the fields within GEOSX");
  RegisterViewWrapper(viewKeyStruct::scaleString, &m_scale, false)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDefaultValue(1.)->setDescription("Scale the coordinates of the vertices");
}

PAMELAMeshGenerator::~PAMELAMeshGenerator()
{}

void PAMELAMeshGenerator::GenerateElementRegions( DomainPartition& domain )
{}

void PAMELAMeshGenerator::PostProcessInput()
{
  m_pamelaMesh =
    std::unique_ptr< PAMELA::Mesh >
      ( PAMELA::MeshFactory::makeMesh( m_filePath ) );
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

void PAMELAMeshGenerator::GenerateMesh( DomainPartition * const domain )
{
  GEOS_LOG_RANK_0("Writing into the GEOSX mesh data structure");
  domain->getMetisNeighborList() = m_pamelaMesh->getNeighborList();
  ManagedGroup * const meshBodies = domain->GetGroup( std::string( "MeshBodies" ));
  MeshBody * const meshBody = meshBodies->RegisterGroup<MeshBody>( this->getName() );

  //TODO for the moment we only consider on mesh level "Level0"
  MeshLevel * const meshLevel0 = meshBody->RegisterGroup<MeshLevel>( std::string( "Level0" ));
  NodeManager * nodeManager = meshLevel0->getNodeManager();
  CellBlockManager * cellBlockManager = domain->GetGroup<CellBlockManager>( keys::cellManager );


  // Use the PartMap of PAMELA to get the mesh
  auto polyhedronPartMap = std::get<0>( PAMELA::getPolyhedronPartMap( m_pamelaMesh.get(), 0));

  // Vertices are written first
  r1_array const & X = nodeManager->referencePosition();
  nodeManager->resize(m_pamelaMesh->get_PointCollection()->size_all());
    R1Tensor xMax( std::numeric_limits< real64 >::min(),
                   std::numeric_limits< real64 >::min(),
                   std::numeric_limits< real64 >::min());

    R1Tensor xMin( std::numeric_limits< real64 >::max(),
                   std::numeric_limits< real64 >::max(),
                   std::numeric_limits< real64 >::max());
  for( auto verticesIterator : *m_pamelaMesh->get_PointCollection()) {
    localIndex vertexLocalIndex = verticesIterator->get_localIndex();
    globalIndex vertexGlobalIndex = verticesIterator->get_globalIndex();
    real64 * const pointData = X[verticesIterator->get_localIndex()].Data();
    pointData[0] = verticesIterator->get_coordinates().x;
    pointData[1] = verticesIterator->get_coordinates().y;
    pointData[2] = verticesIterator->get_coordinates().z;
    nodeManager->m_localToGlobalMap[vertexLocalIndex] = vertexGlobalIndex;
    for( int i = 0; i < 3 ; i++ )
    {
      if( pointData[i] > xMax[i] )
      {
        xMax[i] = pointData[i] ;
      }
      if( pointData[i] < xMin[i] )
      {
        xMin[i] = pointData[i] ;
      }
    }
  }
  xMax -= xMin;
  meshBody->setGlobalLengthScale( std::fabs( xMax.L2_Norm() ) );
  
  // First loop which iterate on the regions
  array1d< globalIndex > globalIndexRegionOffset( polyhedronPartMap.size() +1 );
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
        auto nbCells = cellBlockPAMELA->SubCollection.size_owned();
        if ( nbCells == 0 ) continue;
        cellBlock =
          cellBlockManager->GetGroup( keys::cellBlocks )->RegisterGroup<CellBlock>( regionIndexStr + "_" + cellBlockName);
        cellBlock -> SetElementType("C3D8");
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
        auto nbCells = cellBlockPAMELA->SubCollection.size_owned();
        if ( nbCells == 0 ) continue;
        cellBlock =
          cellBlockManager->GetGroup( keys::cellBlocks )->RegisterGroup<CellBlock>( regionIndexStr + "_" + cellBlockName);
        cellBlock -> SetElementType("C3D4");
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
        auto nbCells = cellBlockPAMELA->SubCollection.size_owned();
        if ( nbCells == 0 ) continue;
        cellBlock =
          cellBlockManager->GetGroup( keys::cellBlocks )->RegisterGroup<CellBlock>( regionIndexStr + "_" + cellBlockName);
        cellBlock -> SetElementType("C3D6");
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
        auto nbCells = cellBlockPAMELA->SubCollection.size_owned();
        if ( nbCells == 0 ) continue;
        cellBlock =
          cellBlockManager->GetGroup( keys::cellBlocks )->RegisterGroup<CellBlock>( regionIndexStr + "_" + cellBlockName);
        cellBlock -> SetElementType("C3D5");
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
      /// Import ppt
      if( cellBlock != nullptr )
      {
        for( localIndex fieldIndex = 0; fieldIndex < m_fieldNamesInGEOSX.size(); fieldIndex++ )
        {
          auto meshProperty = regionPtr->FindVariableByName( m_fieldsToImport[fieldIndex] );
          auto dimension = meshProperty->Dimension;
          if( dimension == PAMELA::VARIABLE_DIMENSION::SCALAR )
          {
            real64_array & property = cellBlock->AddProperty< real64_array >( m_fieldNamesInGEOSX[fieldIndex] );
            GEOS_ERROR_IF(property.size() != integer_conversion< localIndex >( meshProperty->size() ),
                "Viewer size (" << property.size() << ") mismatch with property size in PAMELA ("
                << meshProperty->size() << ") on " <<cellBlock->getName() );
            for( int cellIndex = 0; cellIndex < property.size(); cellIndex++ )
            {
              property[cellIndex] = meshProperty->get_data( cellIndex )[0];
            }
          }
          else if( dimension == PAMELA::VARIABLE_DIMENSION::VECTOR )
          {
            array2d< real64 > & property = cellBlock->AddProperty< array2d< real64 > >( m_fieldNamesInGEOSX[fieldIndex] );
            property.resizeDimension<1>( static_cast< localIndex >( dimension ) );
            GEOS_ERROR_IF(property.size() != integer_conversion< localIndex >( meshProperty->size() ),
                "Viewer size (" << property.size() << ") mismatch with property size in PAMELA ("
                << meshProperty->size() << ") on " <<cellBlock->getName() );
            for( int cellIndex = 0; cellIndex < cellBlock->size(); cellIndex++ )
            {
              for( int dim = 0; dim < 3; dim++ )
              {
                property[cellIndex][dim] = meshProperty->get_data( cellIndex )[dim];
              }
            }
          }
          else
          {
            GEOS_ERROR("Dimension of " <<  m_fieldNamesInGEOSX[fieldIndex] << " is not supported for import in GEOSX");
          }
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
