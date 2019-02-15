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

#include "FieldSpecificationManager.hpp"

#include "FieldSpecificationBase.hpp"
#include "constitutive/ConstitutiveManager.hpp"

#include "mesh/MeshBody.hpp"

#include "managers/NumericalMethodsManager.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteElement/ElementLibrary/FiniteElement.h"
#include "codingUtilities/StringUtilities.hpp"

#include "managers/DomainPartition.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;
FieldSpecificationManager::FieldSpecificationManager( string const & name, ManagedGroup * const parent ):
  ManagedGroup( name, parent )
{
  setInputFlags(InputFlags::OPTIONAL);
}


FieldSpecificationManager * FieldSpecificationManager::get()
{
  static FieldSpecificationManager bcman( "FieldSpecifications", nullptr );
  return &bcman;
}


FieldSpecificationManager::~FieldSpecificationManager()
{
  // TODO Auto-generated destructor stub
}

ManagedGroup * FieldSpecificationManager::CreateChild( string const & childKey, string const & childName )
{
  std::unique_ptr<FieldSpecificationBase> bc = FieldSpecificationBase::CatalogInterface::Factory( childKey, childName, this );
  return this->RegisterGroup( childName, std::move( bc ) );
}


void FieldSpecificationManager::ExpandObjectCatalogs()
{
  // During schema generation, register one of each type derived from BoundaryConditionBase here
  for (auto& catalogIter: FieldSpecificationBase::GetCatalog())
  {
    CreateChild( catalogIter.first, catalogIter.first );
  }
}


void FieldSpecificationManager::ApplyInitialConditions( ManagedGroup * domain ) const
{

  Apply( 0.0, domain, "", "",
         [&]( FieldSpecificationBase const * const bc,
         string const &,
         set<localIndex> const & targetSet,
         ManagedGroup * const targetGroup,
         string const fieldName )
    {
      bc->ApplyFieldValue<FieldSpecificationEqual>( targetSet, 0.0, targetGroup, fieldName );
    } );
}
void FieldSpecificationManager::ApplyInitialConditionsFromMesh( dataRepository::ManagedGroup * domain,
                                                                MeshManager * meshManager) const
{
  std::cout << "This mesh mananager has nb subgroups " << meshManager->numSubGroups() << std::endl;
  Apply2( 0.0, domain, "", "", meshManager,
         [&]( FieldSpecificationBase const * const bc,
         string const &,
         set<localIndex> const & targetSet,
         ManagedGroup * const targetGroup,
         string const fieldName,
         const real64_array& fieldArray)
    {
      bc->ApplyFieldValue<FieldSpecificationEqual>( targetSet, 0.0, targetGroup, fieldName, fieldArray );
    } );
}

} /* namespace geosx */
