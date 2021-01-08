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
  m_refElevation( 0.0 ),
  m_refGravCoef( 0.0 ),
  m_currentControl( Control::BHP ),
  m_targetBHP( 0.0 ),
  m_targetTotalRate( 0.0 ),
  m_targetPhaseRate( 0.0 ),
  m_targetPhaseName( "" ),
  m_useSurfaceConditions( 0 )
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

  registerWrapper( viewKeyStruct::targetTotalRateString, &m_targetTotalRate )->
    setDefaultValue( 0.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Target total volumetric rate" );

  registerWrapper( viewKeyStruct::targetPhaseRateString, &m_targetPhaseRate )->
    setDefaultValue( 0.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Target phase volumetric rate" );

  registerWrapper( viewKeyStruct::targetPhaseNameString, &m_targetPhaseName )->
    setDefaultValue( "" )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Name of the target phase" );


  registerWrapper( viewKeyStruct::refElevString, &m_refElevation )->
    setDefaultValue( -1 )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Reference elevation where BHP control is enforced" );

  registerWrapper( viewKeyStruct::useSurfaceConditionsString, &m_useSurfaceConditions )->
    setDefaultValue( 0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Flag to specify whether rates are checked at surface or reservoir conditions.\n"
                    "Equal to 1 for surface conditions, and to 0 for reservoir conditions" );

  registerWrapper( viewKeyStruct::injectionStreamString, &m_injectionStream )->
    setDefaultValue( -1 )->
    setSizedFromParent( 0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Global component densities for the injection stream" );

}


WellControls::~WellControls()
{}

void WellControls::SwitchToBHPControl( real64 const & val )
{
  m_currentControl = Control::BHP;
  m_targetBHP = val;
}

void WellControls::SwitchToTotalRateControl( real64 const & val )
{
  m_currentControl = Control::TOTALVOLRATE;
  m_targetTotalRate = val;
}

void WellControls::SwitchToPhaseRateControl( real64 const & val )
{
  m_currentControl = Control::PHASEVOLRATE;
  m_targetPhaseRate = val;
}

void WellControls::PostProcessInput()
{
  // 3.a) check target BHP
  if( m_targetBHP < 0 )
  {
    GEOSX_ERROR( "Target bottom-hole pressure for well "<< getName() << " is negative" );
  }

  // 3.b) check target rates
  if( m_targetTotalRate < 0 )
  {
    GEOSX_ERROR( "Target rate for well "<< getName() << " is negative" );
  }
  if( m_targetPhaseRate < 0 )
  {
    GEOSX_ERROR( "Target oil rate for well "<< getName() << " is negative" );
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

  // 5) check the flag for surface / reservoir conditions
  GEOSX_ERROR_IF( m_useSurfaceConditions != 0 && m_useSurfaceConditions != 1,
                  "The flag to select surface/reservoir conditions must be equal to 0 or 1" );

  // 6) check that at least one rate constraint has been defined
  GEOSX_ERROR_IF( m_targetPhaseRate <= 0.0 && m_targetTotalRate <= 0.0,
                  "You need to specify a phase rate constraint for producers, and a total rate constraint for injectors" );

}


void WellControls::InitializePostInitialConditions_PreSubGroups( Group * const GEOSX_UNUSED_PARAM( rootGroup ) )
{
  // for a producer, the solvers compute negative rates, so we adjust the input here
  if( GetType() == Type::PRODUCER
      && (m_targetPhaseRate > 0.0 || m_targetTotalRate > 0.0) )
  {
    m_targetTotalRate *= -1;
    m_targetPhaseRate *= -1;
  }
}

} //namespace geosx
