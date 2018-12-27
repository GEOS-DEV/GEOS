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
 * BoundaryConditionManager.cpp
 *
 *  Created on: May 26, 2017
 *      Author: rrsettgast
 */

#include "BoundaryConditionManager.hpp"
#include "BoundaryConditionBase.hpp"
#include "constitutive/ConstitutiveManager.hpp"

#include "mesh/MeshBody.hpp"

#include "managers/NumericalMethodsManager.hpp"
#include "finiteElement/FiniteElementSpaceManager.hpp"
#include "finiteElement/ElementLibrary/FiniteElement.h"
#include "codingUtilities/StringUtilities.hpp"

#include "managers/DomainPartition.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;
BoundaryConditionManager::BoundaryConditionManager( string const & name, ManagedGroup * const parent ):
  ManagedGroup( name, parent )
{}


BoundaryConditionManager * BoundaryConditionManager::get()
{
  static BoundaryConditionManager bcman( "bcMan", nullptr );
  return &bcman;
}


BoundaryConditionManager::~BoundaryConditionManager()
{
  // TODO Auto-generated destructor stub
}

ManagedGroup * BoundaryConditionManager::CreateChild( string const & childKey, string const & childName )
{
  std::unique_ptr<BoundaryConditionBase> bc = BoundaryConditionBase::CatalogInterface::Factory( childKey, childName, this );
  return this->RegisterGroup( childName, std::move( bc ) );
}



void BoundaryConditionManager::ApplyInitialConditions( ManagedGroup * domain ) const
{

  ApplyBoundaryCondition( 0.0, domain, "", "",
                          [&]( BoundaryConditionBase const * const bc,
                               string const &,
                               set<localIndex> const & targetSet,
                               ManagedGroup * const targetGroup,
                               string const fieldName )
    {
      bc->ApplyBoundaryConditionToField<BcEqual>( targetSet, 0.0, targetGroup, fieldName );
    } );
}

} /* namespace geosx */
