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

/**
 * @file FiniteElementSpaceManager.cpp
 */

#include "FiniteElementDiscretizationManager.hpp"

#include "FiniteElementDiscretization.hpp"

namespace geosx
{
using namespace dataRepository;

FiniteElementDiscretizationManager::FiniteElementDiscretizationManager( string const & name, Group * const parent ):
  Group(name,parent)
{
  setInputFlags(InputFlags::OPTIONAL);
}

FiniteElementDiscretizationManager::~FiniteElementDiscretizationManager()
{
  // TODO Auto-generated destructor stub
}


Group * FiniteElementDiscretizationManager::CreateChild( string const & childKey, string const & childName )
{
  // These objects should probably not be registered on managed group...
  std::unique_ptr<Group> fem = Group::CatalogInterface::Factory( childKey, childName, this );
  return this->RegisterGroup( childName, std::move(fem) );
}


void FiniteElementDiscretizationManager::ExpandObjectCatalogs()
{
  // During schema generation, register one of each type derived from Group here
  for (auto& catalogIter: Group::GetCatalog())
  {
    CreateChild( catalogIter.first, catalogIter.first );
  }
}


} /* namespace geosx */
