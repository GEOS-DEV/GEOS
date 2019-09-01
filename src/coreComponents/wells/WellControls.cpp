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
 * @file WellControls.cpp
 *
 */

#include "WellControls.hpp"

#include "dataRepository/InputFlags.hpp"

namespace geosx
{


using namespace dataRepository;

WellControls::WellControls(string const & name, Group * const parent)
  : Group(name, parent),
    m_typeString(""),
    m_type( Type::PRODUCER ),
    m_refWellElemDepth( 0.0 ),
    m_refWellElemIndex( -1 ),
    m_inputControlString(""),
    m_currentControl( Control::BHP ),
    m_targetBHP( 0.0 ),
    m_targetRate( 0.0 )
{

  registerWrapper( viewKeyStruct::typeString, &m_typeString, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Well type (producer/injector)");

  registerWrapper( viewKeyStruct::controlString, &m_inputControlString, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Well control (BHP/gasRate/oilRate/waterRate)");

  registerWrapper( viewKeyStruct::targetBHPString, &m_targetBHP, false )->
    setDefaultValue(-1)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Target bottom-hole pressure");

  registerWrapper( viewKeyStruct::targetRateString, &m_targetRate, false )->
    setDefaultValue(-1)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Target rate");

  registerWrapper( viewKeyStruct::refWellElemDepthString, &m_refWellElemDepth, false )->
    setDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Reference depth for well bottom hole pressure");

  registerWrapper( viewKeyStruct::injectionStreamString, &m_injectionStream, false )->
    setDefaultValue(-1)->
    setSizedFromParent(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Global component densities for the injection stream");

}

WellControls::~WellControls()
{

}

real64 WellControls::GetInjectionStream( localIndex ic ) const
{ 
  real64 compFrac = -1;
  if (ic < m_injectionStream.size()) 
  {
    compFrac = m_injectionStream[ic];
  }
  return compFrac;
}

void WellControls::SetControl( Control control, 
                                    real64 const & val )
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


void WellControls::PostProcessInput()
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

  if (m_inputControlString == "BHP")
  {
    m_currentControl = Control::BHP;
  }
  else if (m_inputControlString == "gasRate")
  {
    m_currentControl = Control::GASRATE;
    GEOS_ERROR("This control is not implemented yet");
  }
  else if (m_inputControlString == "oilRate")
  {
    m_currentControl = Control::OILRATE;
    GEOS_ERROR("This control is not implemented yet");
  }
  else if (m_inputControlString == "waterRate")
  {
    m_currentControl = Control::WATERRATE;
    GEOS_ERROR("This control is not implemented yet");
  }
  else if (m_inputControlString == "liquidRate")
  {
    m_currentControl = Control::LIQUIDRATE;
  }
  else
  {
    GEOS_ERROR("Invalid initial well control: " << m_inputControlString);
  }

  // 3.a) check target BHP
  if (m_targetBHP < 0)
  {
    GEOS_ERROR("Target bottom-hole pressure for well "<< getName() << " is negative");
  }

  // 3.b) check target rate
  if (m_targetRate < 0) // choose a default value if negative
  {
    GEOS_ERROR("Target rate for well "<< getName() << " is negative");
  }

  // 4) check injection stream
  if (!m_injectionStream.empty())
  {
    real64 sum = 0.0;
    for (localIndex ic = 0; ic < m_injectionStream.size(); ++ic)
    {
      GEOS_ERROR_IF( m_injectionStream[ic] < 0.0 || m_injectionStream[ic] > 1.0,
                     "Invalid injection stream for well " << getName() );
      sum += m_injectionStream[ic];
    }
    GEOS_ERROR_IF( std::abs(1.0 - sum) > std::numeric_limits<real64>::epsilon(),
                   "Invalid injection stream for well " << getName() );
  }
}

void WellControls::InitializePostInitialConditions_PreSubGroups( Group * const rootGroup )
{
  // for a producer, the solvers compute negative rates, so we adjust the input here
  if (GetType() == Type::PRODUCER && m_targetRate > 0.0)
  {
    m_targetRate *= -1;
  }

}

void WellControls::Debug() const 
{
  std::cout << "Name = " << getName() << std::endl;
  std::cout << "Type = " << static_cast<integer>(m_type) << std::endl;
  std::cout << "Current control = " << static_cast<integer>(m_currentControl) << std::endl;
  std::cout << "Target BHP = " << m_targetBHP << std::endl;
  std::cout << "Target rate = " << m_targetRate << std::endl;
  for (localIndex ic = 0; ic < m_injectionStream.size(); ++ic)
  {
    std::cout << "Injection composition: injectionStream[" << ic << "] = " << m_injectionStream[ic] << std::endl;
  }
}

} //namespace geosx
