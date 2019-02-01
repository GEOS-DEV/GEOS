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

#include "BasisFunctionManager.hpp"
#include "BasisBase.hpp"
#include "dataRepository/RestartFlags.hpp"

namespace geosx
{
using namespace dataRepository;

BasisFunctionManager::BasisFunctionManager( string const & name, ManagedGroup * const parent ):
  ManagedGroup(name,parent)
{
  setInputFlags(InputFlags::OPTIONAL);
}

BasisFunctionManager::~BasisFunctionManager()
{
  // TODO Auto-generated destructor stub
}


ManagedGroup * BasisFunctionManager::CreateChild( string const & childKey, string const & childName )
{
  std::unique_ptr<BasisBase> basis = BasisBase::CatalogInterface::Factory( childKey, childName, this );
  return this->RegisterGroup<BasisBase>( childName, std::move(basis) );
}


void BasisFunctionManager::ExpandObjectCatalogs()
{
  // During schema generation, register one of each type derived from BasisBase here
  for (auto& catalogIter: BasisBase::GetCatalog())
  {
    CreateChild( catalogIter.first, catalogIter.first );
  }
}


} /* namespace geosx */
