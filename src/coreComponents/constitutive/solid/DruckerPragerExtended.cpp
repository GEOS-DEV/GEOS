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
 *  @file DruckerPragerExtended.cpp
 */

#include "DruckerPragerExtended.hpp"
#include "SolidFields.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

DruckerPragerExtended::DruckerPragerExtended( string const & name, Group * const parent ):
  ElasticIsotropic( name, parent )
{
  // register default values

  registerWrapper( viewKeyStruct::defaultInitialFrictionAngleString(), &m_defaultInitialFrictionAngle ).
    setApplyDefaultValue( 30.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Initial friction angle (degrees)" );

  registerWrapper( viewKeyStruct::defaultResidualFrictionAngleString(), &m_defaultResidualFrictionAngle ).
    setApplyDefaultValue( 30.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Residual friction angle (degrees)" );

  registerWrapper( viewKeyStruct::defaultDilationRatioString(), &m_defaultDilationRatio ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Dilation ratio [0,1] (ratio = tan dilationAngle / tan frictionAngle)" );

  registerWrapper( viewKeyStruct::defaultHardeningString(), &m_defaultHardening ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Hardening parameter (hardening rate is faster for smaller values)" );

  registerWrapper( viewKeyStruct::defaultCohesionString(), &m_defaultCohesion ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Initial cohesion" );

  // register fields

  registerField< fields::solid::initialFriction >( &m_initialFriction );

  registerField< fields::solid::residualFriction >( &m_residualFriction );

  registerField< fields::solid::dilation >( &m_dilationRatio );

  registerField< fields::solid::pressureIntercept >( &m_pressureIntercept );

  registerField< fields::solid::hardening >( &m_hardening );

  registerField< fields::solid::state >( &m_newState );

  registerField< fields::solid::oldState >( &m_oldState );
}


void DruckerPragerExtended::allocateConstitutiveData( Group & parent, localIndex const numPts )
{
  m_newState.resize( 0, numPts );
  m_oldState.resize( 0, numPts );

  ElasticIsotropic::allocateConstitutiveData( parent, numPts );
}


void DruckerPragerExtended::postInputInitialization()
{
  ElasticIsotropic::postInputInitialization();

  GEOS_THROW_IF( m_defaultCohesion < 0,
                 getFullName() << ": Negative cohesion value detected",
                 InputError, getDataContext() );
  GEOS_THROW_IF( m_defaultInitialFrictionAngle < 0,
                 getFullName() << ": Negative initial friction angle detected",
                 InputError, getDataContext() );
  GEOS_THROW_IF( m_defaultResidualFrictionAngle < 0,
                 getFullName() << ": Negative residual friction angle detected",
                 InputError, getDataContext() );
  GEOS_THROW_IF( m_defaultDilationRatio < 0,
                 getFullName() << ": Dilation ratio out of [0,1] range detected",
                 InputError, getDataContext() );
  GEOS_THROW_IF( m_defaultDilationRatio > 1,
                 getFullName() << ": Dilation ratio out of [0,1] range detected",
                 InputError, getDataContext() );
  GEOS_THROW_IF( m_defaultHardening < 0,
                 getFullName() << ": Negative hardening parameter detected",
                 InputError, getDataContext() );

  // convert from Mohr-Coulomb constants to Drucker-Prager constants, assuming DP
  // passes through the triaxial tension corners of the MC surface.
  // see Borja (2013) p. 75

  real64 phi_i = m_defaultInitialFrictionAngle * M_PI / 180;
  real64 phi_r = m_defaultResidualFrictionAngle * M_PI / 180;
  real64 F_i = 6 * sin( phi_i ) / ( 3 - sin( phi_i ) );
  real64 F_r = 6 * sin( phi_r ) / ( 3 - sin( phi_r ) );
  real64 C = 6 * m_defaultCohesion * cos( phi_i ) / ( 3 - sin( phi_i ) );
  real64 P = C / F_i;

  // set results as array default values

  getWrapper< array1d< real64 > >( fields::solid::initialFriction::key()).
    setApplyDefaultValue( F_i );

  getWrapper< array1d< real64 > >( fields::solid::residualFriction::key() ).
    setApplyDefaultValue( F_r );

  getWrapper< array1d< real64 > >( fields::solid::pressureIntercept::key() ).
    setApplyDefaultValue( P );

  getWrapper< array1d< real64 > >( fields::solid::dilation::key() ).
    setApplyDefaultValue( m_defaultDilationRatio );

  getWrapper< array1d< real64 > >( fields::solid::hardening::key() ).
    setApplyDefaultValue( m_defaultHardening );
}


void DruckerPragerExtended::saveConvergedState() const
{
  SolidBase::saveConvergedState(); // TODO: not ideal, as we have separate loops for base and derived data

  localIndex const numE = numElem();
  localIndex const numQ = numQuad();

  arrayView2d< real64 const > newState = m_newState;
  arrayView2d< real64 > oldState = m_oldState;

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < numQ; ++q )
    {
      oldState( k, q ) = newState( k, q );
    }
  } );
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, DruckerPragerExtended, string const &, Group * const )
}
} /* namespace geos */
