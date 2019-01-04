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


#include "MeshManager.hpp"

#include "MeshGeneratorBase.hpp"
#include "MPI_Communications/SpatialPartition.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

MeshManager::MeshManager( std::string const & name,
                          ManagedGroup * const parent ):
  ManagedGroup( name, parent)
{}

MeshManager::~MeshManager()
{}

ManagedGroup * MeshManager::CreateChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0("Adding Mesh: " << childKey << ", " << childName);
  std::unique_ptr<MeshGeneratorBase> solver = MeshGeneratorBase::CatalogInterface::Factory( childKey, childName, this );
  return this->RegisterGroup<MeshGeneratorBase>( childName, std::move(solver) );
}


void MeshManager::ExpandObjectCatalogs()
{
  // During schema generation, register one of each type derived from MeshGeneratorBase here
  for (auto& catalogIter: MeshGeneratorBase::GetCatalog())
  {
    CreateChild( catalogIter.first, catalogIter.first );
  }
}


void MeshManager::GenerateMeshes( DomainPartition * const domain )
{
  forSubGroups<MeshGeneratorBase>([&]( MeshGeneratorBase * meshGen ) -> void
  {
    meshGen->GenerateMesh( domain );
  });

}


void MeshManager::GenerateMeshLevels( DomainPartition * const domain )
{
  this->forSubGroups<MeshGeneratorBase>([&]( MeshGeneratorBase * meshGen ) -> void
  {
    string meshName = meshGen->getName();
    domain->getMeshBodies()->RegisterGroup<MeshBody>(meshName)->CreateMeshLevel(0);
  });
}


} /* namespace geosx */
