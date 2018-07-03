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
 * NewFunctionManager.cpp
 *
 *  Created on: June 16, 2017
 *      Author: sherman
 */

#include "NewFunctionManager.hpp"
#include "dataRepository/ManagedGroup.hpp"
#include "FunctionBase.hpp"
#include "common/DataTypes.hpp"

namespace geosx
{

using namespace dataRepository;


NewFunctionManager::NewFunctionManager( const std::string& name,
                                        ManagedGroup * const parent ):
  ManagedGroup( name, parent )
{}

NewFunctionManager::~NewFunctionManager()
{
  // TODO Auto-generated destructor stub
}

void NewFunctionManager::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName("Functions");
  docNode->setSchemaType("Node");
  docNode->setShortDescription("Function manager");
}


void NewFunctionManager::CreateChild( string const & functionCatalogKey,
                                      string const & functionName )
{
  std::cout << "   " << functionCatalogKey << ": " << functionName << std::endl;
  std::unique_ptr<FunctionBase> function = FunctionBase::CatalogInterface::Factory( functionCatalogKey, functionName, this );
  this->RegisterGroup<FunctionBase>( functionName, std::move(function) );
}

} /* namespace ANST */
