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
 * @file PerforationManager.cpp
 *
 */

#include "PerforationManager.hpp"
#include "Perforation.hpp"

namespace geosx
{

using namespace dataRepository;

class WellConnectionManager
{

};

PerforationManager::PerforationManager(string const & name, ManagedGroup * const parent)
  : ManagedGroup(name, parent)
{

}

PerforationManager::~PerforationManager()
{

}

void PerforationManager::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName("Perforations");
  docNode->setSchemaType("UniqueNode");
}

void PerforationManager::CreateChild(string const & childKey, string const & childName)
{
  if ( childKey == groupKeyStruct::perforationString )
  {
    this->RegisterGroup<Perforation>(childName);
  }
  else
  {
    GEOS_ERROR("Unrecognized node: " << childKey);
  }
}


} //namespace geosx