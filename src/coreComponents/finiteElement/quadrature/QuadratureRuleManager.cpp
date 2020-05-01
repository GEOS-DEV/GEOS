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

#include "QuadratureRuleManager.hpp"
#include "QuadratureBase.hpp"

namespace geosx
{
using namespace dataRepository;

QuadratureRuleManager::QuadratureRuleManager( string const & name, Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
}

QuadratureRuleManager::~QuadratureRuleManager()
{
  // TODO Auto-generated destructor stub
}


Group * QuadratureRuleManager::CreateChild( string const & childKey, string const & childName )
{
  std::unique_ptr< QuadratureBase > quadrature = QuadratureBase::CatalogInterface::Factory( childKey, childName, this );
  return this->RegisterGroup< QuadratureBase >( childName, std::move( quadrature ) );
}


void QuadratureRuleManager::ExpandObjectCatalogs()
{
  // During schema generation, register one of each type derived from QuadratureBase here
  for( auto & catalogIter: QuadratureBase::GetCatalog())
  {
    CreateChild( catalogIter.first, catalogIter.first );
  }
}


} /* namespace geosx */
