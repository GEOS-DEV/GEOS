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
    m_wellElementManager(groupKeyStruct::wellElementsString, this ),
    m_connectionData( groupKeyStruct::connectionDataString, this ),
    m_perforationData( groupKeyStruct::perforationDataString, this ),
    m_perforationManager( groupKeyStruct::perforationsString, this ),
    m_referenceDepth( 0.0 ),
    m_typeString( "producer" ),
    m_type( Type::PRODUCER ),
    m_controlString( "BHP" ),
    m_currentControl( Control::BHP ),
    m_targetBHP( 0.0 ),
    m_targetRate( 0.0 )
{
  RegisterViewWrapper( viewKeyStruct::typeString, &m_typeString, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Well type (producer/injector)");

  RegisterViewWrapper( viewKeyStruct::controlString, &m_controlString, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Well control (BHP/gasRate/oilRate/waterRate)");

  // TODO: here, double-check what should be REQUIRED/OPTIONAL
  
  RegisterViewWrapper( viewKeyStruct::referenceDepthString, &m_referenceDepth, false )->
    setDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Reference depth for well bottom hole pressure");

  RegisterViewWrapper( viewKeyStruct::targetBHPString, &m_targetBHP, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Target bottom-hole pressure");

  RegisterViewWrapper( viewKeyStruct::targetRateString, &m_targetRate, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Well control (BHP/gasRate/oilRate/waterRate)");

  RegisterViewWrapper( viewKeyStruct::injectionStreamString, &m_injectionStream, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Global component densities for the injection stream");
  
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
  // set well type
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
  
  // set initial control type
  if (m_controlString == "BHP")
  {
    m_currentControl = Control::BHP;
  }
  else if (m_controlString == "gasRate")
  {
    m_currentControl = Control::GASRATE;
  }
  else if (m_controlString == "oilRate")
  {
    m_currentControl = Control::OILRATE;
  }
  else if (m_controlString == "waterRate")
  {
    m_currentControl = Control::WATERRATE;
  }
  else if (m_controlString == "liquidRate")
  {
    m_currentControl = Control::LIQUIDRATE;
  }
  else
  {
    GEOS_ERROR("Invalid initial well control: " << m_controlString);
  }

}


R1Tensor const & Well::getGravityVector() const
{
  return getParent()->group_cast<WellManager const *>()->getGravityVector();
}

} //namespace geosx
