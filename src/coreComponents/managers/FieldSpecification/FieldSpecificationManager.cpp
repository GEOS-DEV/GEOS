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

#include "FieldSpecificationManager.hpp"

#include "codingUtilities/StringUtilities.hpp"
#include "constitutive/ConstitutiveManager.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;
FieldSpecificationManager::FieldSpecificationManager( string const & name, Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
}

FieldSpecificationManager::~FieldSpecificationManager()
{
  // TODO Auto-generated destructor stub
}

Group * FieldSpecificationManager::createChild( string const & childKey, string const & childName )
{
  std::unique_ptr< FieldSpecificationBase > bc = FieldSpecificationBase::CatalogInterface::factory( childKey, childName, this );
  return &this->registerGroup( childName, std::move( bc ) );
}


void FieldSpecificationManager::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from BoundaryConditionBase here
  for( auto & catalogIter: FieldSpecificationBase::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}


void FieldSpecificationManager::applyInitialConditions( DomainPartition & domain ) const
{
  apply( 0.0, domain, "", "",
         [&]( FieldSpecificationBase const & bc,
              string const &,
              SortedArrayView< localIndex const > const & targetSet,
              Group & targetGroup,
              string const fieldName )
  {
    bc.applyFieldValue< FieldSpecificationEqual >( targetSet, 0.0, targetGroup, fieldName );
  } );
}

} /* namespace geosx */
