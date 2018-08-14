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
                                          ManagedGroup * const parent ):
  ManagedGroup( name, parent )
{}

ConstitutiveManager::~ConstitutiveManager()
{}

void ConstitutiveManager::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName( "Constitutive" );
  docNode->setSchemaType( "Node" );
}

void ConstitutiveManager::CreateChild( string const & childKey, string const & childName )
{
  std::unique_ptr<ConstitutiveBase> material = ConstitutiveBase::CatalogInterface::Factory( childKey, childName, this );
  ConstitutiveBase * newMaterial = this->RegisterGroup<ConstitutiveBase>( childName, std::move( material ) );
}

void ConstitutiveManager::HangConstitutiveRelation( string const & constitutiveRelationInstanceName,
                                                    dataRepository::ManagedGroup * const parent,
                                                    localIndex const numConstitutivePointsPerParentIndex ) const
{
  ConstitutiveBase const * const
  constitutiveRelation = GetConstitituveRelation( constitutiveRelationInstanceName );

  std::unique_ptr<ConstitutiveBase>
  material = constitutiveRelation->DeliverClone( constitutiveRelationInstanceName, parent );

  material->AllocateConstitutiveData( parent,
                                      numConstitutivePointsPerParentIndex );

  dataRepository::ManagedGroup * constitutiveGroup = parent->GetGroup( groupKeyStruct::constitutiveModelsString );
  if( constitutiveGroup == nullptr )
  {
    constitutiveGroup = parent->RegisterGroup( groupKeyStruct::constitutiveModelsString );
  }

  constitutiveGroup->RegisterGroup<ConstitutiveBase>( constitutiveRelationInstanceName,
                                                      std::move( material ) );


}

}

} /* namespace geosx */
