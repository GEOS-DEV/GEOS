/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 *  @file CoulombContact.cpp
 */

#include "CoulombContact.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

CoulombContact::CoulombContact( string const & name, Group * const parent ):
  ContactBase( name, parent ),
  m_cohesion(),
  m_frictionAngle(),
  m_frictionCoefficient(),
  m_elasticSlip()
{
  registerWrapper( viewKeyStruct::cohesionString(), &m_cohesion ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Cohesion" );

  registerWrapper( viewKeyStruct::frictionAngleString(), &m_frictionAngle ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Friction angle (in radians)" );

  registerWrapper( viewKeyStruct::frictionCoefficientString(), &m_frictionCoefficient ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Friction coefficient" );

  registerWrapper( viewKeyStruct::elasticSlipString(), &m_elasticSlip ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Elastic Slip" );

}

CoulombContact::~CoulombContact()
{}

void CoulombContact::postProcessInput()
{
  GEOSX_THROW_IF( m_frictionCoefficient < 0.0 && m_frictionAngle < 0,
                  getCatalogName() << " " << getName() << ": Both friction angle and friction coefficient are less than zero. Values: "
                                   << m_frictionAngle << ", " << m_frictionCoefficient << ". Invalid input.",
                  InputError );

  // Step 1: compute a tentative tangent of the friction angle
  real64 const frictionCoefficient = ( m_frictionAngle >= 0.0 ) ? std::tan( m_frictionAngle ) : -1.0;

  // Step 2: check that inputs are consistent, and if so, set the definitive friction coefficient
  if( m_frictionCoefficient >= 0.0 )
  {
    if( frictionCoefficient >= 0.0 )
    {
      GEOSX_THROW_IF( LvArray::math::abs( m_frictionCoefficient - frictionCoefficient ) > 1e1*std::numeric_limits< real64 >::epsilon(),
                      getCatalogName() << " " << getName() << ": Provided friction angle and friction coefficient do not match: "
                                       << m_frictionCoefficient << ", " << frictionCoefficient << ". Invalid input.",
                      InputError );
    }
  }
  else
  {
    m_frictionCoefficient = frictionCoefficient;
  }

  GEOSX_THROW_IF( m_frictionCoefficient < 0.0,
                  getCatalogName() << " " << getName() << ": The provided friction coefficient is less than zero. Value: " << m_frictionCoefficient,
                  InputError );

}

void CoulombContact::allocateConstitutiveData( Group & parent,
                                               localIndex const numConstitutivePointsPerParentIndex )
{
  m_elasticSlip.resize( 0, 2 );

  ContactBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}


CoulombContactUpdates CoulombContact::createKernelWrapper() const
{
  return CoulombContactUpdates( m_penaltyStiffness,
                                m_shearStiffness,
                                *m_apertureTable,
                                m_cohesion,
                                m_frictionCoefficient,
                                m_elasticSlip );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CoulombContact, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geosx */
