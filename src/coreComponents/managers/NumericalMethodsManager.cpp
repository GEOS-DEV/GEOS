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

/**
 * @file NumericalMethodsManager.cpp
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

Group * NumericalMethodsManager::CreateChild( string const & GEOSX_UNUSED_PARAM( childKey ), string const & GEOSX_UNUSED_PARAM( childName ) )
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
  GEOSX_ERROR("Can't find subgroup named " + name + " in " + this->getName());
  return nullptr;
}



} /* namespace geosx */
