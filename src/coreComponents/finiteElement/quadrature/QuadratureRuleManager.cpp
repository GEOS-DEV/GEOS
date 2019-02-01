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

#include "QuadratureRuleManager.hpp"
#include "QuadratureBase.hpp"

namespace geosx
{
using namespace dataRepository;

QuadratureRuleManager::QuadratureRuleManager( string const & name, ManagedGroup * const parent ):
  ManagedGroup(name,parent)
{
  setInputFlags(InputFlags::OPTIONAL);
}

QuadratureRuleManager::~QuadratureRuleManager()
{
  // TODO Auto-generated destructor stub
}


ManagedGroup * QuadratureRuleManager::CreateChild( string const & childKey, string const & childName )
{
  std::unique_ptr<QuadratureBase> quadrature = QuadratureBase::CatalogInterface::Factory( childKey, childName, this );
  return this->RegisterGroup<QuadratureBase>( childName, std::move(quadrature) );
}


void QuadratureRuleManager::ExpandObjectCatalogs()
{
  // During schema generation, register one of each type derived from QuadratureBase here
  for (auto& catalogIter: QuadratureBase::GetCatalog())
  {
    CreateChild( catalogIter.first, catalogIter.first );
  }
}


} /* namespace geosx */
