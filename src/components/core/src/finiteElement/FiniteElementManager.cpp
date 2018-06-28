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
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
 * FiniteElementManager.cpp
 *
 *  Created on: Apr 18, 2017
 *      Author: rrsettgast
 */

#include "FiniteElementManager.hpp"
#include "basis/BasisFunctionManager.hpp"
#include "quadrature/QuadratureRuleManager.hpp"
#include "FiniteElementSpaceManager.hpp"


namespace geosx
{
using namespace dataRepository;

FiniteElementManager::FiniteElementManager( string const & name, ManagedGroup * const parent ):
  ManagedGroup(name,parent)
{
  this->RegisterGroup<BasisFunctionManager>(keys::basisFunctions);
  this->RegisterGroup<QuadratureRuleManager>(keys::quadratureRules);
  this->RegisterGroup<FiniteElementSpaceManager>(keys::finiteElementSpaces);
}

FiniteElementManager::~FiniteElementManager()
{
  // TODO Auto-generated destructor stub
}

void FiniteElementManager::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName("NumericalMethods");
  docNode->setSchemaType("Node");
}

void FiniteElementManager::CreateChild( string const & childKey, string const & childName )
{
}



} /* namespace geosx */
