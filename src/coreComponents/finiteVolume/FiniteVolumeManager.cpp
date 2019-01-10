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
 * @file FiniteVolumeManager.cpp
 *
 */

#include "FiniteVolumeManager.hpp"

#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"

namespace geosx
{

using namespace dataRepository;


FiniteVolumeManager::FiniteVolumeManager(string const &name, ManagedGroup *const parent)
  : ManagedGroup(name, parent)
{
  setSchemaFlags(SchemaFlags::UNIQUE_NODE);
}

FiniteVolumeManager::~FiniteVolumeManager()
{

}

ManagedGroup * FiniteVolumeManager::CreateChild(string const &childKey, string const &childName)
{
  std::unique_ptr<FluxApproximationBase> approx = FluxApproximationBase::CatalogInterface::Factory(childKey, childName, this);
  return this->RegisterGroup<FluxApproximationBase>(childName, std::move(approx));
}


void FiniteVolumeManager::ExpandObjectCatalogs()
{
  // During schema generation, register one of each type derived from FluxApproximationBase here
  for (auto& catalogIter: FluxApproximationBase::GetCatalog())
  {
    CreateChild( catalogIter.first, catalogIter.first );
  }
}


FluxApproximationBase const * FiniteVolumeManager::getFluxApproximation(std::string const &name) const
{
  return this->GetGroup<FluxApproximationBase>(name);
}


} // namespace geosx
