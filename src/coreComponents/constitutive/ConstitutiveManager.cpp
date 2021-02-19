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
 * @file ConstitutiveManager.cpp
 */

#include "ConstitutiveManager.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


ConstitutiveManager::ConstitutiveManager( string const & name,
                                          Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
}

ConstitutiveManager::~ConstitutiveManager()
{}


Group * ConstitutiveManager::createChild( string const & childKey, string const & childName )
{
  std::unique_ptr< ConstitutiveBase > material = ConstitutiveBase::CatalogInterface::factory( childKey, childName, this );
  return &registerGroup< ConstitutiveBase >( childName, std::move( material ) );
}


void ConstitutiveManager::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from ConstitutiveBase here
  for( auto & catalogIter: ConstitutiveBase::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}


ConstitutiveBase &
ConstitutiveManager::hangConstitutiveRelation( string const & constitutiveRelationInstanceName,
                                               dataRepository::Group * const parent,
                                               localIndex const numConstitutivePointsPerParentIndex ) const
{
  ConstitutiveBase const & constitutiveRelation = getConstitutiveRelation( constitutiveRelationInstanceName );

  std::unique_ptr< ConstitutiveBase > material = constitutiveRelation.deliverClone( constitutiveRelationInstanceName, parent );

  material->allocateConstitutiveData( *parent,
                                      numConstitutivePointsPerParentIndex );

  dataRepository::Group * constitutiveGroup = parent->getGroupPointer( groupKeyStruct::constitutiveModelsString() );
  if( constitutiveGroup == nullptr )
  {
    constitutiveGroup = &parent->registerGroup( groupKeyStruct::constitutiveModelsString() ).setSizedFromParent( 1 );
    constitutiveGroup->resize( parent->size() );
  }


  ConstitutiveBase &
  rval = constitutiveGroup->registerGroup< ConstitutiveBase >( constitutiveRelationInstanceName, std::move( material ) );
  rval.setSizedFromParent( 1 );
  rval.resize( constitutiveGroup->size() );
  return rval;
}

}

} /* namespace geosx */
