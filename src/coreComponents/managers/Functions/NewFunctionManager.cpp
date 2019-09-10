/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file NewFunctionManager.cpp
 */

#include "NewFunctionManager.hpp"

#include "dataRepository/Group.hpp"
#include "FunctionBase.hpp"
#include "common/DataTypes.hpp"

namespace geosx
{

using namespace dataRepository;


NewFunctionManager::NewFunctionManager( const std::string& name,
                                        Group * const parent ):
  Group( name, parent )
{
  setInputFlags(InputFlags::OPTIONAL);
}

NewFunctionManager::~NewFunctionManager()
{
  // TODO Auto-generated destructor stub
}


Group * NewFunctionManager::CreateChild( string const & functionCatalogKey,
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
