/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


#include "MeshManager.hpp"

#include "mesh/mpiCommunications/SpatialPartition.hpp"
#include "generators/MeshGeneratorBase.hpp"
#include "common/TimingMacros.hpp"

namespace geosx
{

using namespace dataRepository;

MeshManager::MeshManager( string const & name,
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
  std::unique_ptr< MeshGeneratorBase > solver = MeshGeneratorBase::CatalogInterface::factory( childKey, childName, this );
  return &this->registerGroup< MeshGeneratorBase >( childName, std::move( solver ) );
}


void MeshManager::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from MeshGeneratorBase here
  for( auto & catalogIter: MeshGeneratorBase::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}


void MeshManager::generateMeshes( DomainPartition & domain )
{
  forSubGroups< MeshGeneratorBase >( [&]( MeshGeneratorBase & meshGen )
  {
    meshGen.generateMesh( domain );
  } );
}


void MeshManager::generateMeshLevels( DomainPartition & domain )
{
  this->forSubGroups< MeshGeneratorBase >( [&]( MeshGeneratorBase & meshGen )
  {
    string const & meshName = meshGen.getName();

    // THIS IS A HACK
    if( meshName.find( "well" ) == string::npos )
    {
      domain.getMeshBodies().registerGroup< MeshBody >( meshName ).createMeshLevel( MeshLevel::groupStructKeys::baseDiscretizationString() );
    }
  } );
}

void MeshManager::importFields( DomainPartition & domain )
{
  forSubGroups< MeshGeneratorBase >( [&]( MeshGeneratorBase & meshGen )
  {
    meshGen.importFields( domain );
  } );
}


} /* namespace geosx */
