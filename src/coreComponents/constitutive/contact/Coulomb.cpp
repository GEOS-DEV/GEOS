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
 *  @file Coulomb.cpp
 */

#include "Coulomb.hpp"

namespace geosx
{
using namespace dataRepository;
namespace constitutive
{

Coulomb::Coulomb( std::string const & name, Group * const parent ):
  ContactRelationBase( name, parent ),
  m_postProcessed( false ),
  m_cohesion(),
  m_frictionAngle(),
  m_frictionCoefficient()
{
  registerWrapper( viewKeyStruct::cohesionString, &m_cohesion )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Cohesion" );

  registerWrapper( viewKeyStruct::frictionAngleString, &m_frictionAngle )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setRestartFlags( RestartFlags::NO_WRITE )->
    setDescription( "Friction Angle (in radians)" );

  registerWrapper( viewKeyStruct::frictionCoefficientString, &m_frictionCoefficient )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Friction Coefficient" );
}


Coulomb::~Coulomb()
{}

void
Coulomb::DeliverClone( string const & name,
                       Group * const parent,
                       std::unique_ptr< ConstitutiveBase > & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique< Coulomb >( name, parent );
  }
  ConstitutiveBase::DeliverClone( name, parent, clone );
  Coulomb * const newConstitutiveRelation = dynamic_cast< Coulomb * >(clone.get());

  newConstitutiveRelation->m_postProcessed = false;
  newConstitutiveRelation->m_cohesion = m_cohesion;
  newConstitutiveRelation->m_frictionAngle = m_frictionAngle;
  newConstitutiveRelation->m_frictionCoefficient = m_frictionCoefficient;
}

real64 Coulomb::limitTangentialTractionNorm( real64 const normalTraction ) const
{
  return ( m_cohesion - normalTraction * m_frictionCoefficient );
}

real64 Coulomb::dLimitTangentialTractionNorm_dNormalTraction( real64 const GEOSX_UNUSED_PARAM( normalTraction ) ) const
{
  return ( m_frictionCoefficient );
}

static real64 const machinePrecision = std::numeric_limits< real64 >::epsilon();

void Coulomb::PostProcessInput()
{
  if( !m_postProcessed )
  {
    GEOSX_ERROR_IF( m_frictionCoefficient < 0.0 && m_frictionAngle < 0,
                    "Both friction angle and friction coefficient are less than zero. Values: "
                    << m_frictionAngle
                    << ", "
                    << m_frictionCoefficient
                    << ". Invalid input." );
    real64 frictionCoefficient = -1.0;
    if( m_frictionAngle >= 0.0 )
    {
      // Compute the tangent of the friction angle just once
      frictionCoefficient = std::tan( m_frictionAngle );
    }

    if( m_frictionCoefficient >= 0.0 )
    {
      if( frictionCoefficient >= 0.0 )
      {
        GEOSX_ERROR_IF( std::fabs( m_frictionCoefficient - frictionCoefficient ) > 1.e+1*machinePrecision,
                        "Provided friction angle and friction coefficient do not match: "
                        << m_frictionCoefficient
                        << ", "
                        << frictionCoefficient
                        << ". Invalid input." );
      }
    }
    else
    {
      m_frictionCoefficient = frictionCoefficient;
    }

    GEOSX_ERROR_IF( m_frictionCoefficient < 0.0,
                    "The provided friction coefficient is less than zero. Value: " << m_frictionCoefficient );
  }

  m_postProcessed = true;
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, Coulomb, std::string const &, Group * const )
}
} /* namespace geosx */
