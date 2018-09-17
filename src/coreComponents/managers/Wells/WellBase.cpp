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
 * @file WellBase.cpp
 *
 */

#include "WellBase.hpp"

#include "Perforation.hpp"
#include "WellManager.hpp"

namespace geosx
{


WellBase::WellBase(string const & name, dataRepository::ManagedGroup * const parent)
  : ObjectManagerBase(name, parent),
    m_referenceDepth(0),
    m_type()
{
  RegisterViewWrapper( viewKeys.referenceDepth.Key(), &m_referenceDepth, false );
  RegisterViewWrapper( viewKeys.type.Key(), &m_type, false );
}

WellBase::~WellBase()
{

}

void WellBase::FillDocumentationNode()
{
  ObjectManagerBase::FillDocumentationNode();

  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->AllocateChildNode( viewKeys.referenceDepth.Key(),
                              viewKeys.referenceDepth.Key(),
                              -1,
                              "real64",
                              "real64",
                              "Reference bottom hole depth",
                              "Reference bottom hole depth",
                              "0",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeys.type.Key(),
                              viewKeys.type.Key(),
                              -1,
                              "string",
                              "string",
                              "Well type (producer/injector)",
                              "Well type (producer/injector)",
                              "REQUIRED",
                              "",
                              0,
                              1,
                              0 );
}

void WellBase::CreateChild(string const & childKey, string const & childName)
{
  //std::cout << "Adding child: " << childKey << ", " << childName << std::endl;
  if (childKey == "Perforation")
  {
    RegisterGroup<Perforation>(childName);
  }
  else
  {
    GEOS_ERROR("Unrecognized child: " << childKey);
  }
}

void WellBase::InitializePostSubGroups(dataRepository::ManagedGroup * const group)
{
  ObjectManagerBase::InitializePostSubGroups(group);
}

R1Tensor const & WellBase::getGravityVector() const
{
  return dynamic_cast<WellManager const *>(getParent())->getGravityVector();
}

void WellBase::FinalInitialization(dataRepository::ManagedGroup * const group)
{
  // nothing yet
}


} //namespace geosx