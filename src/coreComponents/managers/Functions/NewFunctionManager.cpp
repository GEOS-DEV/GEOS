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
 * NewFunctionManager.cpp
 *
 *  Created on: June 16, 2017
 *      Author: sherman
 */

#include "NewFunctionManager.hpp"
#include "dataRepository/ManagedGroup.hpp"
#include "FunctionBase.hpp"
#include "common/DataTypes.hpp"

namespace geosx
{

using namespace dataRepository;


NewFunctionManager::NewFunctionManager( const std::string& name,
                                        ManagedGroup * const parent ):
  ManagedGroup( name, parent )
{
  setSchemaFlags(SchemaFlags::UNIQUE_NODE);
}

NewFunctionManager::~NewFunctionManager()
{
  // TODO Auto-generated destructor stub
}


ManagedGroup * NewFunctionManager::CreateChild( string const & functionCatalogKey,
                                      string const & functionName )
{
  GEOS_LOG_RANK_0("   " << functionCatalogKey << ": " << functionName);
  std::unique_ptr<FunctionBase> function = FunctionBase::CatalogInterface::Factory( functionCatalogKey, functionName, this );
  return this->RegisterGroup<FunctionBase>( functionName, std::move(function) );
}


void NewFunctionManager::ExpandObjectCatalogs()
{
  // During schema generation, register one of each type derived from FunctionBase here
  for (auto& catalogIter: FunctionBase::GetCatalog())
  {
    CreateChild( catalogIter.first, catalogIter.first );
  }
}

} /* namespace ANST */
