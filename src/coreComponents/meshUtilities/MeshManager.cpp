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


#include "MeshManager.hpp"

#include "mpiCommunications/SpatialPartition.hpp"
#include "MeshGeneratorBase.hpp"
#include "common/TimingMacros.hpp"

namespace geosx
{

using namespace dataRepository;

MeshManager::MeshManager( std::string const & name,
                          Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::REQUIRED );
}

MeshManager::~MeshManager()
{}

Group * MeshManager::createChild( string const & childKey, string const & childName )
{
  GEOSX_LOG_RANK_0( "Adding Mesh: " << childKey << ", " << childName );
  std::unique_ptr< MeshGeneratorBase > solver = MeshGeneratorBase::CatalogInterface::Factory( childKey, childName, this );
  return this->registerGroup< MeshGeneratorBase >( childName, std::move( solver ) );
}


void MeshManager::ExpandObjectCatalogs()
{
  // During schema generation, register one of each type derived from MeshGeneratorBase here
  for( auto & catalogIter: MeshGeneratorBase::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}


void MeshManager::GenerateMeshes( DomainPartition * const domain )
{
  forSubGroups< MeshGeneratorBase >( [&]( MeshGeneratorBase & meshGen )
  {
    meshGen.GenerateMesh( domain );
  } );
}


void MeshManager::GenerateMeshLevels( DomainPartition * const domain )
{
  this->forSubGroups< MeshGeneratorBase >( [&]( MeshGeneratorBase & meshGen )
  {
    string meshName = meshGen.getName();

    // THIS IS A HACK
    if( meshName.find( "well" ) == std::string::npos )
    {
      domain->getMeshBodies()->registerGroup< MeshBody >( meshName )->CreateMeshLevel( 0 );
    }
  } );
}


} /* namespace geosx */
