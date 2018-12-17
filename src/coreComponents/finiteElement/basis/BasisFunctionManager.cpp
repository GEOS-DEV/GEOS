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
 * BasisFunctionManager.cpp
 *
 *  Created on: Dec 5, 2017
 *      Author: sherman
 */

#include "BasisFunctionManager.hpp"
#include "BasisBase.hpp"
#include "dataRepository/RestartFlags.hpp"

namespace geosx
{
using namespace dataRepository;

BasisFunctionManager::BasisFunctionManager( string const & name, ManagedGroup * const parent ):
  ManagedGroup(name,parent)
{}

BasisFunctionManager::~BasisFunctionManager()
{
  // TODO Auto-generated destructor stub
}

void BasisFunctionManager::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName("BasisFunctions");
  docNode->setSchemaType("UniqueNode");
}

void BasisFunctionManager::CreateChild( string const & childKey, string const & childName )
{
  std::unique_ptr<BasisBase> basis = BasisBase::CatalogInterface::Factory( childKey );
  this->RegisterViewWrapper( childName, std::move(basis) )->setRestartFlags(RestartFlags::NO_WRITE);
}

// Basis Base is not derived from ManagedGroup, so we need to do this manually:
void BasisFunctionManager::ReadXMLsub( xmlWrapper::xmlNode const & targetNode )
{
  for (xmlWrapper::xmlNode childNode=targetNode.first_child() ; childNode ; childNode=childNode.next_sibling())
  {
    std::string childName = childNode.attribute("name").value();
    BasisBase * basis = this->getPointer<BasisBase>(childName);

    if (basis != nullptr)
    {
      basis->ReadXML(childNode);
    }
  }
}

void BasisFunctionManager::ProcessInputFile( xmlWrapper::xmlNode const & targetNode )
{
  ReadXMLsub( targetNode );
}


} /* namespace geosx */
