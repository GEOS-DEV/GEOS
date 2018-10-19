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
    // TODO : To be written --> Next PR !
}

void PAMELAMeshGenerator::GenerateElementRegions( DomainPartition& domain )
{
    // TODO : To be written --> Next PR !
}

void PAMELAMeshGenerator::ReadXML_PostProcess()
{
    // TODO : To be written --> Next PR !
}

void PAMELAMeshGenerator::RemapMesh(dataRepository::ManagedGroup * const domain)
{
    // TODO : To be written --> Next PR !
}

void PAMELAMeshGenerator::CreateChild( string const & childKey, string const & childName )
{
    // TODO : To be written --> Next PR !
}

void PAMELAMeshGenerator::GenerateMesh( dataRepository::ManagedGroup * const domain )
{
    // TODO : To be written --> Next PR !
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
