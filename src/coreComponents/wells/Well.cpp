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
 * @file Well.cpp
 *
 */

#include "Well.hpp"

#include "WellManager.hpp"

namespace geosx
{

using namespace dataRepository;

Well::Well(string const & name, dataRepository::ManagedGroup * const parent)
  : ObjectManagerBase( name, parent ),
    m_wellElementSubRegion( groupKeyStruct::wellElementDataString, this ),
    m_wellElementManager( groupKeyStruct::wellElementsString, this ),
    m_connectionData( groupKeyStruct::connectionDataString, this ),
    m_perforationData( groupKeyStruct::perforationDataString, this ),
    m_perforationManager( groupKeyStruct::perforationsString, this ),
    m_referenceDepth( 0.0 ),
    m_typeString( "producer" ),
    m_type( Type::PRODUCER )
{
  std::cout << "Well::Well: " << name << std::endl;
  RegisterViewWrapper( viewKeyStruct::referenceDepthString, &m_referenceDepth, false )->
    setDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Reference depth for well bottom hole pressure");

  RegisterViewWrapper( viewKeyStruct::typeString, &m_typeString, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Well type (producer/injector)");

  RegisterGroup( groupKeyStruct::wellElementDataString, &m_wellElementSubRegion, false );
  RegisterGroup( groupKeyStruct::wellElementsString,    &m_wellElementManager,   false );
  
  RegisterGroup( groupKeyStruct::connectionDataString,  &m_connectionData,  false );

  RegisterGroup( groupKeyStruct::perforationDataString, &m_perforationData,    false );
  RegisterGroup( groupKeyStruct::perforationsString,    &m_perforationManager, false );
}

Well::~Well()
{

}


  
dataRepository::ManagedGroup * Well::CreateChild(string const & childKey, string const & childName)
{
  return nullptr;
}

void Well::InitializePostSubGroups( ManagedGroup * const rootGroup )
{
  resize(1);
}
  
void Well::PostProcessInput()
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


R1Tensor const & Well::getGravityVector() const
{
  return getParent()->group_cast<WellManager const *>()->getGravityVector();
}

} //namespace geosx
