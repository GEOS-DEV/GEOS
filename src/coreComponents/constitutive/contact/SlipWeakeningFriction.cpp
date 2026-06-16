/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SlipWeakeningFriction.cpp
 */

#include "SlipWeakeningFriction.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

SlipWeakeningFriction::SlipWeakeningFriction( string const & name, Group * const parent )
  : FrictionBase( name, parent )
{
  registerWrapper( viewKeyStruct::shearStiffnessString(), &m_shearStiffness ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value of the shear elastic stiffness. Units of Pressure/length" );

  registerWrapper( viewKeyStruct::elasticSlipString(), &m_elasticSlip ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Elastic slip" );

  registerWrapper( viewKeyStruct::defaultCohesionString(), &m_defaultCohesion ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default cohesion value" );

  registerWrapper( viewKeyStruct::defaultInitialFrictionCoefficientString(), &m_defaultInitialFrictionCoefficient ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default peak friction coefficient" );

  registerWrapper( viewKeyStruct::defaultResidualFrictionCoefficientString(), &m_defaultResidualFrictionCoefficient ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default residual friction coefficient (must satisfy muResidual <= muPeak)" );

  registerWrapper( viewKeyStruct::defaultDcString(), &m_defaultDc ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default characteristic slip distance over which friction weakens from muPeak to muResidual" );

  registerField< fields::contact::cohesion >( &m_cohesion );
  registerField< fields::contact::initialFrictionCoefficient >( &m_initialFrictionCoefficient );
  registerField< fields::contact::residualFrictionCoefficient >( &m_residualFrictionCoefficient );
  registerField< fields::contact::characteristicSlipDistance >( &m_Dc );
  registerField< fields::contact::cumulativeSlip >( &m_cumulativeSlip );

  registerWrapper( "cumulativeSlipSaved", &m_cumulativeSlipSaved ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Snapshot of cumulative slip at the start of the current time step" );
}

void SlipWeakeningFriction::postInputInitialization()
{
  GEOS_THROW_IF( m_defaultCohesion < 0.0,
                 GEOS_FMT( ": The provided default cohesion is less than zero. Value: {}",
                           m_defaultCohesion ),
                 InputError, getDataContext() );

  GEOS_THROW_IF( m_defaultInitialFrictionCoefficient < 0.0,
                 GEOS_FMT( ": The provided default peak friction coefficient is less than zero. Value: {}",
                           m_defaultInitialFrictionCoefficient ),
                 InputError, getDataContext() );

  GEOS_THROW_IF( m_defaultResidualFrictionCoefficient < 0.0 || m_defaultResidualFrictionCoefficient > m_defaultInitialFrictionCoefficient,
                 GEOS_FMT( ": The residual friction coefficient must satisfy 0 <= muResidual <= muPeak. "
                           "muResidual = {}, muPeak = {}",
                           m_defaultResidualFrictionCoefficient, m_defaultInitialFrictionCoefficient ),
                 InputError, getDataContext() );

  GEOS_THROW_IF( m_defaultDc <= 0.0,
                 GEOS_FMT( ": The characteristic slip distance must be strictly positive. Value: {}",
                           m_defaultDc ),
                 InputError, getDataContext() );

  getField< fields::contact::cohesion >().setApplyDefaultValue( m_defaultCohesion );
  getField< fields::contact::initialFrictionCoefficient >().setApplyDefaultValue( m_defaultInitialFrictionCoefficient );
  getField< fields::contact::residualFrictionCoefficient >().setApplyDefaultValue( m_defaultResidualFrictionCoefficient );
  getField< fields::contact::characteristicSlipDistance >().setApplyDefaultValue( m_defaultDc );
}

void SlipWeakeningFriction::initializePostInitialConditionsPreSubGroups()
{
  FrictionBase::initializePostInitialConditionsPreSubGroups();

  for( localIndex k = 0; k < m_initialFrictionCoefficient.size(); ++k )
  {
    GEOS_THROW_IF( m_initialFrictionCoefficient[k] < 0.0 ||
                   m_residualFrictionCoefficient[k] < 0.0 ||
                   m_residualFrictionCoefficient[k] > m_initialFrictionCoefficient[k] ||
                   m_Dc[k] <= 0.0,
                   GEOS_FMT( "Invalid slip-weakening parameters at element {}: "
                             "muPeak = {}, muResidual = {}, Dc = {}",
                             k, m_initialFrictionCoefficient[k], m_residualFrictionCoefficient[k], m_Dc[k] ),
                   InputError, getDataContext() );
  }
}

void SlipWeakeningFriction::allocateConstitutiveData( Group & parent, localIndex const numPts )
{
  m_elasticSlip.resize( 0, 2 );
  FrictionBase::allocateConstitutiveData( parent, numPts );
}

SlipWeakeningFrictionUpdates SlipWeakeningFriction::createKernelUpdates() const
{
  return SlipWeakeningFrictionUpdates( m_displacementJumpThreshold,
                                      m_shearStiffness,
                                      m_cohesion.toViewConst(),
                                      m_initialFrictionCoefficient.toViewConst(),
                                      m_residualFrictionCoefficient.toViewConst(),
                                      m_Dc.toViewConst(),
                                      m_cumulativeSlip,
                                      m_cumulativeSlipSaved,
                                      m_elasticSlip );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, SlipWeakeningFriction, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geos */
