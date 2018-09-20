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

#include "PerforationManager.hpp"
#include "WellManager.hpp"

namespace geosx
{

WellBase::WellBase(string const & name, dataRepository::ManagedGroup * const parent)
  : ObjectManagerBase(name, parent),
    m_perfManager( groupKeyStruct::perforationsString, this ),
    m_numConnections( 0 ),
    m_referenceDepth( 0.0 ),
    m_typeString( "producer" ),
    m_type( Type::PRODUCER )
{
  RegisterViewWrapper( viewKeys.referenceDepth.Key(), &m_referenceDepth, false );
  RegisterViewWrapper( viewKeys.type.Key(), &m_typeString, false );
  RegisterViewWrapper( viewKeys.connectionElementRegion.Key(), &m_connectionElementRegion, false );
  RegisterViewWrapper( viewKeys.connectionElementSubregion.Key(), &m_connectionElementSubregion, false );
  RegisterViewWrapper( viewKeys.connectionElementIndex.Key(), &m_connectionElementIndex, false );
  RegisterViewWrapper( viewKeys.connectionPerforationIndex.Key(), &m_connectionPerforationIndex, false );

  RegisterGroup( groupKeys.perforations.Key(), &m_perfManager, false );
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
}

void WellBase::InitializePostSubGroups(dataRepository::ManagedGroup * const group)
{
  ObjectManagerBase::InitializePostSubGroups(group);
}

R1Tensor const & WellBase::getGravityVector() const
{
  return getParent()->group_cast<WellManager const *>()->getGravityVector();
}

void WellBase::FinalInitialization(dataRepository::ManagedGroup * const group)
{
  // nothing yet
}

void WellBase::ReadXML_PostProcess()
{
  if (m_typeString == "producer")
  {
    m_type = Type::PRODUCER;
  }
  else if (m_typeString == "injector")
  {
    m_type = Type::INJECTOR;
  }
  else
  {
    GEOS_ERROR("Invalid well type: " << m_typeString);
  }
}


} //namespace geosx
