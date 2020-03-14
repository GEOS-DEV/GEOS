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
  m_frictionAngle(),
  m_frictionAngleUnitOfMeasurement()
{
  registerWrapper( viewKeyStruct::cohesionString, &m_cohesion, false )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Cohesion" );

  registerWrapper( viewKeyStruct::frictionAngleString, &m_frictionAngle, false )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Friction Angle" );

  registerWrapper( viewKeyStruct::frictionAngleUnitOfMeasurementString, &m_frictionAngleUnitOfMeasurement, false )->
    setApplyDefaultValue( "rad" )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Friction Angles Unit of Measurement. Valid inputs are (upper/lower case does not matter):\n"
                    " - rad/radiants (default);\n"
                    " - deg/degrees\n;"
                    " - grad/gradians." );
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
  newConstitutiveRelation->m_frictionAngle = m_frictionAngle;
  newConstitutiveRelation->m_frictionAngleTangent = m_frictionAngleTangent;
  newConstitutiveRelation->m_frictionAngleUnit = m_frictionAngleUnit;
  newConstitutiveRelation->m_frictionAngleUnitOfMeasurement = m_frictionAngleUnitOfMeasurement;
}

real64 MohrCoulomb::limitTangentialTractionNorm( real64 const normalTraction ) const
{
  return ( m_cohesion - normalTraction * m_frictionAngleTangent );
}

real64 MohrCoulomb::dLimitTangentialTractionNorm_dNormalTraction( real64 const GEOSX_UNUSED_PARAM( normalTraction ) ) const
{
  return ( m_frictionAngleTangent );
}

void MohrCoulomb::PostProcessInput()
{
  if( !m_postProcessed )
  {
    string & frictionAngleUnitOfMeasurement = getReference< string >( viewKeyStruct::frictionAngleUnitOfMeasurementString );

    // convert string to lower case
    std::for_each( frictionAngleUnitOfMeasurement.begin(), frictionAngleUnitOfMeasurement.end(), []( char & c ) {
      c = ::tolower( c );
    } );

    if( frictionAngleUnitOfMeasurement == "radians" || frictionAngleUnitOfMeasurement == "rad" )
    {
      m_frictionAngleUnit = AngleUnit::RADIANS;
    }
    else if( frictionAngleUnitOfMeasurement == "degrees" || frictionAngleUnitOfMeasurement == "deg" )
    {
      m_frictionAngleUnit = AngleUnit::DEGREES;
    }
    else if( frictionAngleUnitOfMeasurement == "gradians" || frictionAngleUnitOfMeasurement == "grad" )
    {
      m_frictionAngleUnit = AngleUnit::GRADIANS;
    }
    else
    {
      GEOSX_ERROR( "Provided wrong unit of measurement for friction angle. It can be rad, deg or grad. Provided: " << frictionAngleUnitOfMeasurement );
    }

    real64 & frictionAngle = getReference< real64 >( viewKeyStruct::frictionAngleString );

    switch( m_frictionAngleUnit )
    {
      case AngleUnit::RADIANS:
      {
        // Everything is ok!
        break;
      }
      case AngleUnit::DEGREES:
      {
        frictionAngle *= M_PI / 180;
        break;
      }
      case AngleUnit::GRADIANS:
      {
        frictionAngle *= M_PI / 200;
        break;
      }
    }

    GEOSX_ERROR_IF( frictionAngle < 0.0,
                    "The provided friction angle is less than zero. Value: " << frictionAngle << " [rad]" );

    GEOSX_ERROR_IF( frictionAngle > M_PI_2,
                    "The provided friction angle is larger than pi/2. Value: " << frictionAngle << " [rad]" );

    // Compute the tangent of the friction angle just once
    m_frictionAngleTangent = std::tan( frictionAngle );
  }
  m_postProcessed = true;
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, MohrCoulomb, std::string const &, Group * const )
}
} /* namespace geosx */
