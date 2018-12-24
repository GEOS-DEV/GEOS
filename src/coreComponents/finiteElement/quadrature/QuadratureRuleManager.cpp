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
 * QuadratureRuleManager.cpp
 *
 *  Created on: Dec 5, 2017
 *      Author: sherman
 */

#include "QuadratureRuleManager.hpp"
#include "QuadratureBase.hpp"

namespace geosx
{
using namespace dataRepository;

QuadratureRuleManager::QuadratureRuleManager( string const & name, ManagedGroup * const parent ):
  ManagedGroup(name,parent)
{}

QuadratureRuleManager::~QuadratureRuleManager()
{
  // TODO Auto-generated destructor stub
}


void QuadratureRuleManager::CreateChild( string const & childKey, string const & childName )
{
  std::unique_ptr<QuadratureBase> quadrature = QuadratureBase::CatalogInterface::Factory( childKey );
  this->RegisterViewWrapper( childName, std::move(quadrature) )->setRestartFlags(RestartFlags::NO_WRITE);
}

// Basis Base is not derived from ManagedGroup, so we need to do this manually:
void QuadratureRuleManager::ReadXMLsub( xmlWrapper::xmlNode const & targetNode )
{
  for (xmlWrapper::xmlNode childNode=targetNode.first_child() ; childNode ; childNode=childNode.next_sibling())
  {
    std::string childName = childNode.attribute("name").value();
    QuadratureBase * quadrature = this->getPointer<QuadratureBase>(childName);

    if (quadrature != nullptr)
    {
      quadrature->ReadXML(childNode);
    }
  }
}

void QuadratureRuleManager::ProcessInputFile( xmlWrapper::xmlNode const & targetNode )
{
  ReadXMLsub( targetNode );
}



} /* namespace geosx */
