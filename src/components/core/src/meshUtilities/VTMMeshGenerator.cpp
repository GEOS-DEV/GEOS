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

#ifdef USE_ATK
#include "slic/slic.hpp"
#endif

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
}

VTMMeshGenerator::~VTMMeshGenerator()
{
  // TODO Auto-generated destructor stub
}

void VTMMeshGenerator::FillDocumentationNode()
{
  //MeshLevel * const mesh =
  // domain->group_cast<DomainPartition*>()->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  //NodeManager * const nodes    = mesh->getNodeManager();
  // CellBlockManager * elems =
  // domain->GetGroup<CellBlockManager>(keys::cellManager);

  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName( "MeshFile" );
  docNode->setSchemaType( "Node" );
  docNode->setShortDescription( "a mesh generator" );

  docNode->AllocateChildNode( keys::filePath,
                              keys::filePath,
                              -1,
                              "string",
                              "string",
                              "path to the vtm file",
                              "path to the vtm file",
                              "filePath",
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

void VTMMeshGenerator::ReadXML_PostProcess()
{
  m_fileName = this->getReference<string>(keys::filePath);
  std::cout << "Load m_fileName : " << m_fileName << std::endl;
  m_vtmFile.Load(m_fileName);
  std::cout << "VTM loaded ! " << m_fileName << std::endl;

}



void VTMMeshGenerator::RemapMesh(dataRepository::ManagedGroup * const domain)
{

}

void VTMMeshGenerator::CreateChild( string const & childKey, string const & childName )
{
}

void VTMMeshGenerator::GenerateMesh( dataRepository::ManagedGroup * const domain )
{
  ManagedGroup * const meshBodies = domain->GetGroup(std::string("MeshBodies"));
  MeshBody * const meshBody = meshBodies->RegisterGroup<MeshBody>( this->getName() );
  MeshLevel * const meshLevel0 = meshBody->RegisterGroup<MeshLevel>(std::string("Level0"));
  NodeManager * nodeManager = meshLevel0->getNodeManager();
  nodeManager->SetDocumentationNodes();
  CellBlockManager * elementManager = domain->GetGroup<CellBlockManager>( keys::cellManager );
  ManagedGroup * nodeSets = nodeManager->GetGroup( std::string( "Sets" ) );
  PartitionBase & partition = domain->getReference<PartitionBase>(keys::partitionManager);
  //TODO for the moment we handle only one block
    CellBlock * cellBlock = elementManager->GetGroup(keys::cellBlocks)->RegisterGroup<CellBlock>("PART00001_POLYHEDRON_POLYHEDRON_GROUP_1");
    cellBlock->SetDocumentationNodes();
    cellBlock->RegisterDocumentationNodes();
    cellBlock->ReadXML_PostProcess();
  m_vtmFile.FromVtmToGEOS(meshLevel0);
  std::cout << "Mesh Generated !" << std::endl;
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
