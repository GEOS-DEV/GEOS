/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 *  @file MohrCoulomb.cpp
 */

#include "MohrCoulomb.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace cxx_utilities;
namespace constitutive
{

MohrCoulomb::MohrCoulomb( std::string const & name, Group * const parent ):
  ContactRelationBase( name, parent ),
  m_postProcessed( false ),
  m_cohesion(),
  m_frictionInput(),
  m_frictionInputUnitOfMeasurement()
{
  registerWrapper( viewKeyStruct::cohesionString, &m_cohesion, false )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Cohesion" );

  registerWrapper( viewKeyStruct::frictionInputString, &m_frictionInput, false )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setRestartFlags( RestartFlags::NO_WRITE )->
    setDescription( "Friction Input" );

  registerWrapper( viewKeyStruct::frictionInputUnitOfMeasurementString, &m_frictionInputUnitOfMeasurement, false )->
    setApplyDefaultValue( "rad" )->
    setInputFlag( InputFlags::OPTIONAL )->
    setRestartFlags( RestartFlags::NO_WRITE )->
    setDescription( "Friction Input Unit of Measurement. Valid inputs are (upper/lower case does not matter):\n"
                    " - rad/radians (default);\n"
                    " - deg/degrees." );

  registerWrapper( viewKeyStruct::frictionCoefficientString, &m_frictionCoefficient, false )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Friction Coefficient" );
}


MohrCoulomb::~MohrCoulomb()
{}

void
MohrCoulomb::DeliverClone( string const & name,
                           Group * const parent,
                           std::unique_ptr< ConstitutiveBase > & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique< MohrCoulomb >( name, parent );
  }
  ConstitutiveBase::DeliverClone( name, parent, clone );
  MohrCoulomb * const newConstitutiveRelation = dynamic_cast< MohrCoulomb * >(clone.get());

  newConstitutiveRelation->m_postProcessed = false;
  newConstitutiveRelation->m_cohesion = m_cohesion;
  newConstitutiveRelation->m_frictionInput = m_frictionInput;
  newConstitutiveRelation->m_frictionCoefficient = m_frictionCoefficient;
  newConstitutiveRelation->m_frictionInputUnit = m_frictionInputUnit;
  newConstitutiveRelation->m_frictionInputUnitOfMeasurement = m_frictionInputUnitOfMeasurement;
}

real64 MohrCoulomb::limitTangentialTractionNorm( real64 const normalTraction ) const
{
  return ( m_cohesion - normalTraction * m_frictionCoefficient );
}

real64 MohrCoulomb::dLimitTangentialTractionNorm_dNormalTraction( real64 const GEOSX_UNUSED_PARAM( normalTraction ) ) const
{
  return ( m_frictionCoefficient );
}

void MohrCoulomb::PostProcessInput()
{
  if( !m_postProcessed )
  {
    GEOSX_ERROR_IF( m_frictionInput >= 0.0 && m_frictionCoefficient >= 0.0,
                    "Both friction input and friction coefficient provided. This is not allowed." );
    GEOSX_ERROR_IF( m_frictionInput < 0.0 && m_frictionCoefficient < 0.0,
                    "Neither friction input nor friction coefficient provided. This is not allowed." );

    if( m_frictionCoefficient < 0.0 )
    {
      string & frictionInputUnitOfMeasurement = getReference< string >( viewKeyStruct::frictionInputUnitOfMeasurementString );

      // convert string to lower case
      std::for_each( frictionInputUnitOfMeasurement.begin(), frictionInputUnitOfMeasurement.end(), []( char & c ) {
        c = ::tolower( c );
      } );

      if( frictionInputUnitOfMeasurement == "radians" || frictionInputUnitOfMeasurement == "rad" )
      {
        m_frictionInputUnit = AngleUnit::RADIANS;
      }
      else if( frictionInputUnitOfMeasurement == "degrees" || frictionInputUnitOfMeasurement == "deg" )
      {
        m_frictionInputUnit = AngleUnit::DEGREES;
      }
      else
      {
        GEOSX_ERROR( "Provided wrong unit of measurement for friction angle. It can be rad, deg or grad. Provided: " << frictionInputUnitOfMeasurement );
      }

      switch( m_frictionInputUnit )
      {
        case AngleUnit::RADIANS:
        {
          // Everything is ok!
          break;
        }
        case AngleUnit::DEGREES:
        {
          m_frictionInput *= M_PI / 180;
          break;
        }
      }

      GEOSX_ERROR_IF( m_frictionInput < 0.0,
                      "The provided friction angle is less than zero. Value: " << m_frictionInput << " [rad]" );

      GEOSX_ERROR_IF( m_frictionInput > M_PI_2,
                      "The provided friction angle is larger than pi/2. Value: " << m_frictionInput << " [rad]" );

      // Compute the tangent of the friction angle just once
      m_frictionCoefficient = std::tan( m_frictionInput );
    }
  }

  m_postProcessed = true;
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, MohrCoulomb, std::string const &, Group * const )
}
} /* namespace geosx */
