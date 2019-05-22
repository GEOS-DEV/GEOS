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
 * VTMMeshGenerator.cpp
 *
 *  Created on: Aug 16, 2018
 *      Author: Antoine Mazuyer
 */

#include "VTMMeshGenerator.hpp"

#include "managers/DomainPartition.hpp"

#include "codingUtilities/StringUtilities.hpp"
#include <math.h>
//#include "managers/TableManager.hpp"
//#include "SimpleGeometricObjects.hpp"

#include "MPI_Communications/PartitionBase.hpp"
#include "MPI_Communications/SpatialPartition.hpp"

#include "mesh/MeshBody.hpp"

namespace geosx
{
using namespace dataRepository;

VTMMeshGenerator::VTMMeshGenerator( string const & name, ManagedGroup * const parent ):
  MeshGeneratorBase( name, parent )
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

  RegisterViewWrapper<string>(keys::filePath)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("path to the vtm file");
}

VTMMeshGenerator::~VTMMeshGenerator()
{
  // TODO Auto-generated destructor stub
}



//}
/**
 * @param domain
 */
void VTMMeshGenerator::GenerateElementRegions( DomainPartition& domain )
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

void VTMMeshGenerator::PostProcessInput()
{
  m_fileName = this->getReference<string>(keys::filePath);
  m_vtmFile.Load(m_fileName, true, false);
}



void VTMMeshGenerator::RemapMesh(dataRepository::ManagedGroup * const domain)
{

}

ManagedGroup * VTMMeshGenerator::CreateChild( string const & childKey, string const & childName )
{
  return nullptr;
}

void VTMMeshGenerator::GenerateMesh( DomainPartition * const domain )
{
    /// Basic mesh registration
    ManagedGroup * const meshBodies = domain->GetGroup(std::string("MeshBodies"));
    MeshBody * const meshBody = meshBodies->RegisterGroup<MeshBody>( this->getName() );
    MeshLevel * const meshLevel0 = meshBody->RegisterGroup<MeshLevel>(std::string("Level0"));
    NodeManager * nodeManager = meshLevel0->getNodeManager();
    CellBlockManager * elementManager = domain->GetGroup<CellBlockManager>( keys::cellManager );

    /// Region registrations
    for( localIndex rankBlockIndex = 0 ; rankBlockIndex < m_vtmFile.NumRankBlocks(); rankBlockIndex++) {
        const auto & rankBlock = m_vtmFile.GetRankBlock(rankBlockIndex); 
        for( localIndex meshBlockIndex = 0 ; meshBlockIndex < rankBlock.NumMeshBlocks(); meshBlockIndex++) {
            const auto & meshBlock = rankBlock.GetMeshBlock(meshBlockIndex);
            if( meshBlock.IsARegionBlock() ) {
                const auto & mesh = meshBlock.mesh();
                /// Write nodes
                nodeManager->resize(mesh.NumVertices());
                arrayView1d<R1Tensor> & X = nodeManager->referencePosition();
                for( localIndex a=0 ; a< mesh.NumVertices() ; ++a )
                {
                    real64 * const tensorData = X[a].Data();
                    tensorData[0] = mesh.Vertex(a)[0];
                    tensorData[1] = mesh.Vertex(a)[1];
                    tensorData[2] = mesh.Vertex(a)[2];
                }
                /// Cell blocks registrations
                if( mesh.NumHex() > 0) {
                    CellBlock * cellBlock = elementManager->GetGroup(keys::cellBlocks)->RegisterGroup<CellBlock>("HEX");
                    cellBlock -> SetElementType("C3D8");
                    auto & cellToVertex = cellBlock->nodeList();
                    cellBlock->resize( mesh.NumCells());
                    cellToVertex.resize( mesh.NumCells(), mesh.NumVerticesInCell(0));

                    for( localIndex k=0 ; k<mesh.NumCells() ; ++k )
                    {
                      cellToVertex[k][0] = mesh.CellVertexIndex(k,0);
                      cellToVertex[k][1] = mesh.CellVertexIndex(k,1);
                      cellToVertex[k][2] = mesh.CellVertexIndex(k,3);
                      cellToVertex[k][3] = mesh.CellVertexIndex(k,2);
                      cellToVertex[k][4] = mesh.CellVertexIndex(k,4);
                      cellToVertex[k][5] = mesh.CellVertexIndex(k,5);
                      cellToVertex[k][6] = mesh.CellVertexIndex(k,7);
                      cellToVertex[k][7] = mesh.CellVertexIndex(k,6);
                    }
                }
                if( mesh.NumTetra() > 0) {
                    CellBlock * cellBlock = elementManager->GetGroup(keys::cellBlocks)->RegisterGroup<CellBlock>("TETRA");
                    cellBlock -> SetElementType("C3D4");
                    auto & cellToVertex = cellBlock->nodeList();
                    cellBlock->resize( mesh.NumCells());
                    cellToVertex.resize( mesh.NumCells(), mesh.NumVerticesInCell(0));

                    for( localIndex k=0 ; k<mesh.NumCells() ; ++k )
                    {
                      cellToVertex[k][0] = mesh.CellVertexIndex(k,0);
                      cellToVertex[k][1] = mesh.CellVertexIndex(k,1);
                      cellToVertex[k][2] = mesh.CellVertexIndex(k,2);
                      cellToVertex[k][3] = mesh.CellVertexIndex(k,3);
                    }
                }
                if( mesh.NumPrism() > 0) {
                    CellBlock * cellBlock = elementManager->GetGroup(keys::cellBlocks)->RegisterGroup<CellBlock>("WEDGE");
                    cellBlock -> SetElementType("C3D6");
                    auto & cellToVertex = cellBlock->nodeList();
                    cellBlock->resize( mesh.NumCells());
                    cellToVertex.resize(mesh.NumCells(), mesh.NumVerticesInCell(0));

                    for( localIndex k=0 ; k<mesh.NumCells() ; ++k )
                    {
                      cellToVertex[k][0] = mesh.CellVertexIndex(k,0);
                      cellToVertex[k][1] = mesh.CellVertexIndex(k,1);
                      cellToVertex[k][2] = mesh.CellVertexIndex(k,2);
                      cellToVertex[k][3] = mesh.CellVertexIndex(k,3);
                      cellToVertex[k][4] = mesh.CellVertexIndex(k,4);
                      cellToVertex[k][5] = mesh.CellVertexIndex(k,5);
                    }
                }
                if( mesh.NumPyr() > 0) {
                    CellBlock * cellBlock = elementManager->GetGroup(keys::cellBlocks)->RegisterGroup<CellBlock>("PYR");
                    cellBlock -> SetElementType("C3D5");
                    auto & cellToVertex = cellBlock->nodeList();
                    cellBlock->resize(mesh.NumCells());
                    cellToVertex.resize(mesh.NumCells(), mesh.NumVerticesInCell(0));

                    for( localIndex k=0 ; k<mesh.NumCells() ; ++k )
                    {
                      cellToVertex[k][0] = mesh.CellVertexIndex(k,0);
                      cellToVertex[k][1] = mesh.CellVertexIndex(k,1);
                      cellToVertex[k][2] = mesh.CellVertexIndex(k,2);
                      cellToVertex[k][3] = mesh.CellVertexIndex(k,3);
                      cellToVertex[k][4] = mesh.CellVertexIndex(k,4);
                    }
                }
            }
        }
    }
}

void VTMMeshGenerator::GetElemToNodesRelationInBox( const std::string& elementType,
                                                         const int index[],
                                                         const int& iEle,
                                                         int nodeIDInBox[],
                                                         const int node_size )

{

}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, VTMMeshGenerator, std::string const &, ManagedGroup * const )
}
