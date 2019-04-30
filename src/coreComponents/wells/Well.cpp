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
    m_perforationData( groupKeyStruct::perforationDataString, this ),
    m_perforationManager( groupKeyStruct::perforationsString, this ),
    m_refWellElemDepth( 0.0 ),
    m_refWellElemIndex( -1 ),
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

  RegisterViewWrapper( viewKeyStruct::targetBHPString, &m_targetBHP, false )->
    setDefaultValue(-1)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Target bottom-hole pressure");

  RegisterViewWrapper( viewKeyStruct::targetRateString, &m_targetRate, false )->
    setDefaultValue(-1)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Target rate");

  RegisterViewWrapper( viewKeyStruct::refWellElemIndexString, &m_refWellElemIndex, false );

  RegisterViewWrapper( viewKeyStruct::refWellElemDepthString, &m_refWellElemDepth, false )->
    setDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Reference depth for well bottom hole pressure");

  RegisterViewWrapper( viewKeyStruct::injectionStreamString, &m_injectionStream, false )->
    setDefaultValue(-1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Global component densities for the injection stream");
  
  RegisterGroup( groupKeyStruct::wellElementDataString, &m_wellElementSubRegion, false );
  RegisterGroup( groupKeyStruct::wellElementsString,    &m_wellElementManager,   false );

  RegisterGroup( groupKeyStruct::perforationDataString, &m_perforationData,    false );
  RegisterGroup( groupKeyStruct::perforationsString,    &m_perforationManager, false );

}

Well::~Well()
{
}
 
dataRepository::ManagedGroup * Well::CreateChild( string const & childKey, 
                                                  string const & childName )
{
  return nullptr;
}

real64 Well::getInjectionStream( localIndex ic ) const
{ 
  real64 compFrac = -1;
  if (ic < m_injectionStream.size()) 
  {
    compFrac = m_injectionStream[ic];
  }
  return compFrac;
}
  
void Well::PostProcessInput()
{  
  // 1) set well type

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
  
  // 2) set initial control type

  if (m_controlString == "BHP")
  {
    m_currentControl = Control::BHP;
  }
  else if (m_controlString == "gasRate")
  {
    m_currentControl = Control::GASRATE;
    GEOS_ERROR("This control is not implemented yet");
  }
  else if (m_controlString == "oilRate")
  {
    m_currentControl = Control::OILRATE;
    GEOS_ERROR("This control is not implemented yet");
  }
  else if (m_controlString == "waterRate")
  {
    m_currentControl = Control::WATERRATE;
    GEOS_ERROR("This control is not implemented yet");
  }
  else if (m_controlString == "liquidRate")
  {
    m_currentControl = Control::LIQUIDRATE;
  }
  else
  {
    GEOS_ERROR("Invalid initial well control: " << m_controlString);
  }

  // 3.a) check target BHP
  if (m_targetBHP < 0)
  {
    GEOS_ERROR("Target bottom-hole pressure for well "<< getName() 
            << " is negative");
  }

  // 3.b) check target rate 
  if (m_targetRate < 0) // choose a default value if negative
  {
    GEOS_ERROR("Target rate for well "<< getName() 
            << " is negative");
  }

  // currently there is only one MPI rank per well
  m_wellElementSubRegion.resize( m_wellElementManager.numWellElementsGlobal() );  
}

void Well::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const rootGroup )
{
  // for a producer, the solvers compute negative rates, so we adjust the input here
  if (getType() == Type::PRODUCER && m_targetRate > 0.0)
  {
    m_targetRate *= -1;
  }

}

void Well::setControl( Control control, real64 const & val )
{ 
  m_currentControl = control; 
  if (control == Control::BHP)
  {
    m_targetBHP = val;
  }
  else 
  {
    m_targetRate = val;
  }
}

} //namespace geosx
