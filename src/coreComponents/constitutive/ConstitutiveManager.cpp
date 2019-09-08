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

/*
 * ConstitutiveManager.cpp
 *
 *  Created on: Aug 4, 2016
 *      Author: rrsettgast
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
  setInputFlags(InputFlags::OPTIONAL);
}

ConstitutiveManager::~ConstitutiveManager()
{}


Group * ConstitutiveManager::CreateChild( string const & childKey, string const & childName )
{
  std::unique_ptr<ConstitutiveBase> material = ConstitutiveBase::CatalogInterface::Factory( childKey, childName, this );
  return RegisterGroup<ConstitutiveBase>( childName, std::move( material ) );
}


void ConstitutiveManager::ExpandObjectCatalogs()
{
  // During schema generation, register one of each type derived from ConstitutiveBase here
  for (auto& catalogIter: ConstitutiveBase::GetCatalog())
  {
    CreateChild( catalogIter.first, catalogIter.first );
  }
}


ConstitutiveBase *
ConstitutiveManager::HangConstitutiveRelation( string const & constitutiveRelationInstanceName,
                                               dataRepository::Group * const parent,
                                               localIndex const numConstitutivePointsPerParentIndex ) const
{
  ConstitutiveBase const * const
  constitutiveRelation = GetConstitutiveRelation( constitutiveRelationInstanceName );

  std::unique_ptr<ConstitutiveBase> material;
  constitutiveRelation->DeliverClone( constitutiveRelationInstanceName, parent, material );

  material->AllocateConstitutiveData( parent,
                                      numConstitutivePointsPerParentIndex );

  dataRepository::Group * constitutiveGroup = parent->GetGroup( groupKeyStruct::constitutiveModelsString );
  if( constitutiveGroup == nullptr )
  {
    constitutiveGroup = parent->RegisterGroup( groupKeyStruct::constitutiveModelsString );
  }

  return constitutiveGroup->RegisterGroup<ConstitutiveBase>( constitutiveRelationInstanceName,
                                                             std::move( material ) );


}

}

} /* namespace geosx */
