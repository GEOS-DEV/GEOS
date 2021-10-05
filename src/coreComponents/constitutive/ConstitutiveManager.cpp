/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
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


void
ConstitutiveManager::hangConstitutiveRelation( string const & constitutiveRelationInstanceName,
                                               dataRepository::Group * const parent,
                                               localIndex const numConstitutivePointsPerParentIndex ) const
{
  dataRepository::Group * constitutiveGroup = parent->getGroupPointer( groupKeyStruct::constitutiveModelsString() );
  if( constitutiveGroup == nullptr )
  {
    constitutiveGroup = &parent->registerGroup( groupKeyStruct::constitutiveModelsString() ).setSizedFromParent( 1 );
    constitutiveGroup->resize( parent->size() );
  }

  // 1. Allocate constitutive relation
  // we only register the constitutive relation if it has not been registered yet.
  GEOSX_ERROR_IF( constitutiveGroup->hasGroup( constitutiveRelationInstanceName ),
                  GEOSX_FMT( "Error! The constitutive relation {} has already been registered on the subRegion {}. "
                             "Make sure that the same constitutive model is not listed as a material on a"
                             " region both as a stand-alone one and as part of a compound constitutive model.",
                             constitutiveRelationInstanceName, parent->getName() ) );

  ConstitutiveBase const & constitutiveRelation = getConstitutiveRelation( constitutiveRelationInstanceName );

  std::unique_ptr< ConstitutiveBase > material = constitutiveRelation.deliverClone( constitutiveRelationInstanceName, parent );

  material->allocateConstitutiveData( *parent,
                                      numConstitutivePointsPerParentIndex );

  ConstitutiveBase &
  materialGroup = constitutiveGroup->registerGroup< ConstitutiveBase >( constitutiveRelationInstanceName, std::move( material ) );
  materialGroup.setSizedFromParent( 1 );
  materialGroup.resize( constitutiveGroup->size() );

  // 2. Allocate subrelations (for compound models)
  std::vector< string > const subRelationNames = constitutiveRelation.getSubRelationNames();
  for( string const & subRelationName : subRelationNames )
  {
    // we only want to register the subRelation if it has not been registered yet.
    GEOSX_ERROR_IF( constitutiveGroup->hasGroup( subRelationName ),
                    GEOSX_FMT( "Error! The constitutive relation {} has already been registered on the subRegion {}. "
                               "Make sure that the same constitutive model is not listed as a material on a"
                               " region both as a stand-alone one and as part of a compound constitutive model.",
                               subRelationName, parent->getName() ) );

    ConstitutiveBase const & subRelation = getConstitutiveRelation( subRelationName );

    std::unique_ptr< ConstitutiveBase > constitutiveModel = subRelation.deliverClone( subRelationName, parent );

    constitutiveModel->allocateConstitutiveData( *parent,
                                                 numConstitutivePointsPerParentIndex );


    ConstitutiveBase &
    group = constitutiveGroup->registerGroup< ConstitutiveBase >( subRelationName, std::move( constitutiveModel ) );
    group.setSizedFromParent( 1 );
    group.resize( constitutiveGroup->size() );
  }
}

} /* namespace constitutive */

} /* namespace geosx */
