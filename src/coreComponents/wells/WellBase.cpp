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

using namespace dataRepository;

WellBase::WellBase(string const & name, dataRepository::ManagedGroup * const parent)
  : ObjectManagerBase( name, parent ),
    m_segManager( groupKeyStruct::segmentsString, this ),
    m_perfManager( groupKeyStruct::perforationsString, this ),
    m_referenceDepth( 0.0 ),
    m_typeString( "producer" ),
    m_type( Type::PRODUCER )
{
  RegisterViewWrapper( viewKeyStruct::referenceDepthString, &m_referenceDepth, false )->
    setDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Reference depth for well bottom hole pressure");

  RegisterViewWrapper( viewKeyStruct::typeString, &m_typeString, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Well type (producer/injector)");

  RegisterGroup( groupKeyStruct::segmentsString, &m_segManager, false );
  
  RegisterGroup( groupKeyStruct::perforationsString, &m_perfManager, false );
}

WellBase::~WellBase()
{

}

dataRepository::ManagedGroup * WellBase::CreateChild(string const & childKey, string const & childName)
{
  // child groups will be registered explicitly, rather than via input file
  return nullptr;
}

R1Tensor const & WellBase::getGravityVector() const
{
  return getParent()->group_cast<WellManager const *>()->getGravityVector();
}

void WellBase::PostProcessInput()
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
