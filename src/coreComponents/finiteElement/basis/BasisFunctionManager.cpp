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

#include "BasisFunctionManager.hpp"
#include "BasisBase.hpp"
#include "dataRepository/RestartFlags.hpp"

namespace geosx
{
using namespace dataRepository;

BasisFunctionManager::BasisFunctionManager( string const & name, Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
}

BasisFunctionManager::~BasisFunctionManager()
{
  // TODO Auto-generated destructor stub
}


Group * BasisFunctionManager::CreateChild( string const & childKey, string const & childName )
{
  std::unique_ptr< BasisBase > basis = BasisBase::CatalogInterface::Factory( childKey, childName, this );
  return this->RegisterGroup< BasisBase >( childName, std::move( basis ) );
}


void BasisFunctionManager::ExpandObjectCatalogs()
{
  // During schema generation, register one of each type derived from BasisBase here
  for( auto & catalogIter: BasisBase::GetCatalog())
  {
    CreateChild( catalogIter.first, catalogIter.first );
  }
}


} /* namespace geosx */
