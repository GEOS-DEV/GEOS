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
 * GeometricObjectManager.cpp
 *
 *  Created on: Oct 18, 2017
 *      Author: sherman
 */

#include "GeometricObjectManager.hpp"

#include "../../../../cxx-utilities/src/src/DocumentationNode.hpp"
#include "SimpleGeometricObjectBase.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

GeometricObjectManager::GeometricObjectManager( std::string const & name,
                                                ManagedGroup * const parent ):
  ManagedGroup( name, parent)
{}

GeometricObjectManager::~GeometricObjectManager()
{}


void GeometricObjectManager::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName("Geometry");
  docNode->setSchemaType("UniqueNode");
  docNode->setShortDescription("Geometric object manager");
}


void GeometricObjectManager::CreateChild( string const & childKey, string const & childName )
{
  std::cout << "Adding Geometric Object: " << childKey << ", " << childName << std::endl;
  std::unique_ptr<SimpleGeometricObjectBase> geometriObject = SimpleGeometricObjectBase::CatalogInterface::Factory( childKey, childName, this );
  this->RegisterGroup<SimpleGeometricObjectBase>( childName, std::move(geometriObject) );
}



} /* namespace geosx */
