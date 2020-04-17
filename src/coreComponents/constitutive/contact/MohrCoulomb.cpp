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
  m_frictionInputUnitOfMeasurement( "coefficient" )
{
  registerWrapper( viewKeyStruct::cohesionString, &m_cohesion, false )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Cohesion" );

  registerWrapper( viewKeyStruct::frictionInputString, &m_frictionInput, false )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::REQUIRED )->
    setRestartFlags( RestartFlags::NO_WRITE )->
    setDescription( "Friction Input" );

  registerWrapper( viewKeyStruct::frictionInputUnitOfMeasurementString, &m_frictionInputUnitOfMeasurement, false )->
    setApplyDefaultValue( "rad" )->
    setInputFlag( InputFlags::REQUIRED )->
    setRestartFlags( RestartFlags::NO_WRITE )->
    setDescription( "Friction Input Unit of Measurement. Valid inputs are (upper/lower case does not matter):\n"
                    " - coefficient: the friction coefficient (default);\n"
                    " - rad/radians: the angle in radians;\n"
                    " - deg/degrees: the angle in degrees." );

  registerWrapper( viewKeyStruct::frictionCoefficientString, &m_frictionCoefficient, false )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::FALSE )->
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
    string & frictionInputUnitOfMeasurement = getReference< string >( viewKeyStruct::frictionInputUnitOfMeasurementString );

    // convert string to lower case
    std::for_each( frictionInputUnitOfMeasurement.begin(), frictionInputUnitOfMeasurement.end(), []( char & c ) {
      c = ::tolower( c );
    } );

    if( frictionInputUnitOfMeasurement == "coefficient" )
    {
      m_frictionInputUnit = InputUnit::COEFFICIENT;
    }
    else if( frictionInputUnitOfMeasurement == "radians" || frictionInputUnitOfMeasurement == "rad" )
    {
      m_frictionInputUnit = InputUnit::RADIANS;
    }
    else if( frictionInputUnitOfMeasurement == "degrees" || frictionInputUnitOfMeasurement == "deg" )
    {
      m_frictionInputUnit = InputUnit::DEGREES;
    }
    else
    {
      GEOSX_ERROR( "Provided wrong unit of measurement for friction coefficient. It can be coefficient, rad or deg. Provided: "
                   << frictionInputUnitOfMeasurement );
    }

    switch( m_frictionInputUnit )
    {
      case InputUnit::COEFFICIENT:
      {
        // Everything is ok!
        m_frictionCoefficient = m_frictionInput;
        break;
      }
      case InputUnit::RADIANS:
      {
        GEOSX_ERROR_IF( m_frictionInput > M_PI_2,
                        "The provided friction angle is larger than pi/2. Value: " << m_frictionInput << " [rad]" );
        // Compute the tangent of the friction angle just once
        m_frictionCoefficient = std::tan( m_frictionInput );
        break;
      }
      case InputUnit::DEGREES:
      {
        GEOSX_ERROR_IF( m_frictionInput > 90,
                        "The provided friction angle is larger than 90 degrees. Value: " << m_frictionInput << " [deg]" );
        m_frictionInput *= M_PI / 180;
        // Compute the tangent of the friction angle just once
        m_frictionCoefficient = std::tan( m_frictionInput );
        break;
      }
    }

    GEOSX_ERROR_IF( m_frictionCoefficient < 0.0,
                    "The provided friction coefficient is less than zero. Value: " << m_frictionCoefficient );
  }

  m_postProcessed = true;
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, MohrCoulomb, std::string const &, Group * const )
}
} /* namespace geosx */
