/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/*
 * @file WellControls.cpp
 */

#include "WellControls.hpp"
#include "dataRepository/InputFlags.hpp"

namespace geosx
{

using namespace dataRepository;

WellControls::WellControls( string const & name, Group * const parent )
  : Group( name, parent ),
  m_type( Type::PRODUCER ),
  m_refWellElemIndex( -1 ),
  m_currentControl( Control::BHP ),
  m_targetBHP( 0.0 ),
  m_targetRate( 0.0 )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::typeString, &m_type )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Well type. Valid options:\n* " + EnumStrings< Type >::concat( "\n* " ) );

  registerWrapper( viewKeyStruct::controlString, &m_currentControl )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Well control. Valid options:\n* " + EnumStrings< Control >::concat( "\n* " ) );

  registerWrapper( viewKeyStruct::targetBHPString, &m_targetBHP )->
    setDefaultValue( -1 )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Target bottom-hole pressure" );

  registerWrapper( viewKeyStruct::targetRateString, &m_targetRate )->
    setDefaultValue( -1 )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Target rate" );

  registerWrapper( viewKeyStruct::injectionStreamString, &m_injectionStream )->
    setDefaultValue( -1 )->
    setSizedFromParent( 0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Global component densities for the injection stream" );

}


WellControls::~WellControls()
{}


void WellControls::setControl( Control control,
                               real64 const & val )
{
  m_currentControl = control;
  if( control == Control::BHP )
  {
    m_targetBHP = val;
  }
  else
  {
    m_targetRate = val;
  }
}


void WellControls::postProcessInput()
{
  // 3.a) check target BHP
  if( m_targetBHP < 0 )
  {
    GEOSX_ERROR( "Target bottom-hole pressure for well "<< getName() << " is negative" );
  }

  // 3.b) check target rate
  if( m_targetRate < 0 ) // choose a default value if negative
  {
    GEOSX_ERROR( "Target rate for well "<< getName() << " is negative" );
  }

  // 4) check injection stream
  if( !m_injectionStream.empty())
  {
    real64 sum = 0.0;
    for( localIndex ic = 0; ic < m_injectionStream.size(); ++ic )
    {
      GEOSX_ERROR_IF( m_injectionStream[ic] < 0.0 || m_injectionStream[ic] > 1.0,
                      "Invalid injection stream for well " << getName() );
      sum += m_injectionStream[ic];
    }
    GEOSX_ERROR_IF( std::abs( 1.0 - sum ) > std::numeric_limits< real64 >::epsilon(),
                    "Invalid injection stream for well " << getName() );
  }
}


void WellControls::initializePostInitialConditionsPreSubGroups( Group * const GEOSX_UNUSED_PARAM( rootGroup ) )
{
  // for a producer, the solvers compute negative rates, so we adjust the input here
  if( getType() == Type::PRODUCER && m_targetRate > 0.0 )
  {
    m_targetRate *= -1;
  }
}


void WellControls::debug() const
{
  std::cout << "Name = " << getName() << std::endl;
  std::cout << "Type = " << static_cast< integer >(m_type) << std::endl;
  std::cout << "Current control = " << static_cast< integer >(m_currentControl) << std::endl;
  std::cout << "Target BHP = " << m_targetBHP << std::endl;
  std::cout << "Target rate = " << m_targetRate << std::endl;
  for( localIndex ic = 0; ic < m_injectionStream.size(); ++ic )
  {
    std::cout << "Injection composition: injectionStream[" << ic << "] = " << m_injectionStream[ic] << std::endl;
  }
}

} //namespace geosx
