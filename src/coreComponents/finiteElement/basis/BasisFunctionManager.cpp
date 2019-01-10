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
  setSchemaFlags(SchemaFlags::UNIQUE_NODE);
}

BasisFunctionManager::~BasisFunctionManager()
{
  // TODO Auto-generated destructor stub
}


ManagedGroup * BasisFunctionManager::CreateChild( string const & childKey, string const & childName )
{
  std::unique_ptr<BasisBase> basis = BasisBase::CatalogInterface::Factory( childKey );
  this->RegisterViewWrapper( childName, std::move(basis) )->setRestartFlags(RestartFlags::NO_WRITE);
  return nullptr;
}


void BasisFunctionManager::ExpandObjectCatalogs()
{
  // During schema generation, register one of each type derived from BasisBase here
  for (auto& catalogIter: BasisBase::GetCatalog())
  {
    CreateChild( catalogIter.first, catalogIter.first );
  }
}


// Basis Base is not derived from ManagedGroup, so we need to do this manually:
void BasisFunctionManager::ProcessInputFile( xmlWrapper::xmlNode const & targetNode )
{
  for (xmlWrapper::xmlNode childNode=targetNode.first_child() ; childNode ; childNode=childNode.next_sibling())
  {
    std::string childName = childNode.attribute("name").value();
    BasisBase * basis = this->getPointer<BasisBase>(childName);

    if (basis != nullptr)
    {
      basis->ReadXML(childNode);
    }
  }
}



} /* namespace geosx */
