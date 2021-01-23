/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FunctionManager.cpp
 */

#include "FunctionManager.hpp"

#include "dataRepository/Group.hpp"
#include "FunctionBase.hpp"
#include "common/DataTypes.hpp"

namespace geosx
{

using namespace dataRepository;


FunctionManager::FunctionManager( const std::string & name,
                                  Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
}

FunctionManager::~FunctionManager()
{
  // TODO Auto-generated destructor stub
}


Group * FunctionManager::createChild( string const & functionCatalogKey,
                                      string const & functionName )
{
  GEOSX_LOG_RANK_0( "   " << functionCatalogKey << ": " << functionName );
  std::unique_ptr< FunctionBase > function = FunctionBase::CatalogInterface::factory( functionCatalogKey, functionName, this );
  return this->registerGroup< FunctionBase >( functionName, std::move( function ) );
}


void FunctionManager::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from FunctionBase here
  for( auto & catalogIter: FunctionBase::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}

} /* namespace ANST */
