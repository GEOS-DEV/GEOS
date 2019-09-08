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
 * NumericalMethodsManager.cpp
 *
 *  Created on: Apr 18, 2017
 *      Author: rrsettgast
 */

#include "NumericalMethodsManager.hpp"

#include "finiteElement/basis/BasisFunctionManager.hpp"
#include "finiteElement/quadrature/QuadratureRuleManager.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"

namespace geosx
{
using namespace dataRepository;

NumericalMethodsManager::NumericalMethodsManager( string const & name, Group * const parent ):
  Group(name,parent)
{
  setInputFlags(InputFlags::OPTIONAL);

  this->RegisterGroup<BasisFunctionManager>(keys::basisFunctions);
  this->RegisterGroup<QuadratureRuleManager>(keys::quadratureRules);
  this->RegisterGroup<FiniteElementDiscretizationManager>(keys::finiteElementDiscretizations);
  this->RegisterGroup<FiniteVolumeManager>(keys::finiteVolumeManager);
}

NumericalMethodsManager::~NumericalMethodsManager()
{
  // TODO Auto-generated destructor stub
}

Group * NumericalMethodsManager::CreateChild( string const & childKey, string const & childName )
{
  return nullptr;
}

dataRepository::Group const * NumericalMethodsManager::FindNumericalMethodByName(string const & name) const
{
  for( auto & iterNumericalMethod : this->GetSubGroups() )
  {
    if( iterNumericalMethod.second->getName() == name)
    {
      return iterNumericalMethod.second;
    }
  }
  GEOS_ERROR("Can't find subgroup named " + name + " in " + this->getName());
  return nullptr;
}



} /* namespace geosx */
