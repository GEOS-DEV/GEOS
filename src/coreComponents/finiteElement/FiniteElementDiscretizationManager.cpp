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
 * FiniteElementSpaceManager.cpp
 *
 *  Created on: Dec 5, 2017
 *      Author: sherman
 */

#include "FiniteElementDiscretizationManager.hpp"

#include "FiniteElementDiscretization.hpp"

namespace geosx
{
using namespace dataRepository;

FiniteElementDiscretizationManager::FiniteElementDiscretizationManager( string const & name, ManagedGroup * const parent ):
  ManagedGroup(name,parent)
{
  setInputFlags(InputFlags::OPTIONAL);
}

FiniteElementDiscretizationManager::~FiniteElementDiscretizationManager()
{
  // TODO Auto-generated destructor stub
}


ManagedGroup * FiniteElementDiscretizationManager::CreateChild( string const & childKey, string const & childName )
{
  // These objects should probably not be registered on managed group...
  std::unique_ptr<ManagedGroup> fem = ManagedGroup::CatalogInterface::Factory( childKey, childName, this );
  return this->RegisterGroup( childName, std::move(fem) );
}


void FiniteElementDiscretizationManager::ExpandObjectCatalogs()
{
  // During schema generation, register one of each type derived from ManagedGroup here
  for (auto& catalogIter: ManagedGroup::GetCatalog())
  {
    CreateChild( catalogIter.first, catalogIter.first );
  }
}


} /* namespace geosx */
