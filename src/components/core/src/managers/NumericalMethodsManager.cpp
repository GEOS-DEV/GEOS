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
 * NumericalMethodsManager.cpp
 *
 *  Created on: Apr 18, 2017
 *      Author: rrsettgast
 */

#include "NumericalMethodsManager.hpp"
#include "finiteElement/basis/BasisFunctionManager.hpp"
#include "finiteElement/quadrature/QuadratureRuleManager.hpp"
#include "finiteElement/FiniteElementSpaceManager.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"

namespace geosx
{
using namespace dataRepository;

NumericalMethodsManager::NumericalMethodsManager( string const & name, ManagedGroup * const parent ):
  ManagedGroup(name,parent)
{
  this->RegisterGroup<BasisFunctionManager>(keys::basisFunctions);
  this->RegisterGroup<QuadratureRuleManager>(keys::quadratureRules);
  this->RegisterGroup<FiniteElementSpaceManager>(keys::finiteElementSpaces);
  this->RegisterGroup<FiniteVolumeManager>(keys::finiteVolumeManager);
}

NumericalMethodsManager::~NumericalMethodsManager()
{
  // TODO Auto-generated destructor stub
}

void NumericalMethodsManager::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName("NumericalMethods");
  docNode->setSchemaType("Node");
}

void NumericalMethodsManager::CreateChild( string const & childKey, string const & childName )
{
}



} /* namespace geosx */
