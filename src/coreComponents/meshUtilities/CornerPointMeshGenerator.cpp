/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CornerPointMeshGenerator.cpp
 */

#include "CornerPointMeshGenerator.hpp"

#include "mesh/MeshBody.hpp"
#include "mpiCommunications/PartitionBase.hpp"
#include "mpiCommunications/SpatialPartition.hpp"

namespace geosx
{
using namespace dataRepository;


CornerPointMeshGenerator::CornerPointMeshGenerator( string const & name, Group * const parent ):
  MeshGeneratorBase( name, parent )
{
  registerWrapper( viewKeyStruct::filePathString(), &m_filePath ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Path to the mesh file" );

  string const builderName = "cpMeshBuilder";
  m_cPBuilder = std::make_unique< CPMesh::CPMeshBuilder >( builderName );
}

CornerPointMeshGenerator::~CornerPointMeshGenerator()
{}

void CornerPointMeshGenerator::generateElementRegions( DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

void CornerPointMeshGenerator::postProcessInput()
{
  m_cPBuilder->buildMesh( m_filePath );
  GEOSX_ERROR( "For now, stop the simulation here" );
}

void CornerPointMeshGenerator::remapMesh( dataRepository::Group & GEOSX_UNUSED_PARAM( domain ) )
{
  return;
}

Group * CornerPointMeshGenerator::createChild( string const & GEOSX_UNUSED_PARAM( childKey ),
                                               string const & GEOSX_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

void CornerPointMeshGenerator::generateMesh( DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

void CornerPointMeshGenerator::getElemToNodesRelationInBox( const string & GEOSX_UNUSED_PARAM( elementType ),
                                                            const int GEOSX_UNUSED_PARAM( index )[],
                                                            const int & GEOSX_UNUSED_PARAM( iEle ),
                                                            int GEOSX_UNUSED_PARAM( nodeIDInBox )[],
                                                            const int GEOSX_UNUSED_PARAM( node_size ) )
{}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, CornerPointMeshGenerator, string const &, Group * const )
}
