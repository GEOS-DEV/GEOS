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
 *  @file DruckerPragerExtended.cpp
 */

#include "DruckerPragerExtended.hpp"

namespace geosx
{
using namespace dataRepository;
namespace constitutive
{

DruckerPragerExtended::DruckerPragerExtended( string const & name, 
                                              Group * const parent ):
  ElasticIsotropic( name, parent ),
  m_defaultInitialFrictionAngle(),
  m_defaultResidualFrictionAngle(),
  m_defaultDilationRatio(),
  m_defaultCohesion(),
  m_defaultHardening(),
  m_initialFriction(),
  m_residualFriction(),
  m_dilationRatio(),
  m_pressureIntercept(),
  m_hardening(),
  m_state()
{
  // register default values

  registerWrapper( viewKeyStruct::defaultInitialFrictionAngleString, &m_defaultInitialFrictionAngle )->
    setApplyDefaultValue( 30.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Initial friction angle (degrees)" );

  registerWrapper( viewKeyStruct::defaultResidualFrictionAngleString, &m_defaultResidualFrictionAngle )->
    setApplyDefaultValue( 30.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Residual friction angle (degrees)" );

  registerWrapper( viewKeyStruct::defaultDilationRatioString, &m_defaultDilationRatio )->
    setApplyDefaultValue( 1.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Dilation ratio [0,1] (ratio = tan dilationAngle / tan frictionAngle)" );

  registerWrapper( viewKeyStruct::defaultHardeningString, &m_defaultHardening )->
    setApplyDefaultValue( 1.0e6 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Hardening parameter (hardening rate is faster for smaller values)" );

  registerWrapper( viewKeyStruct::defaultCohesionString, &m_defaultCohesion )->
    setApplyDefaultValue( 0.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Initial cohesion" );

  // register fields

  registerWrapper( viewKeyStruct::initialFrictionString, &m_initialFriction )->
    setApplyDefaultValue( -1 )->
    setDescription( "Initial yield surface slope" );

  registerWrapper( viewKeyStruct::residualFrictionString, &m_residualFriction )->
    setApplyDefaultValue( -1 )->
    setDescription( "Residual yield surface slope" );

  registerWrapper( viewKeyStruct::dilationRatioString, &m_dilationRatio )->
    setApplyDefaultValue( -1 )->
    setDescription( "Plastic potential slope ratio" );

  registerWrapper( viewKeyStruct::pressureInterceptString, &m_pressureIntercept )->
    setApplyDefaultValue( -1 )->
    setDescription( "Pressure point at cone vertex" );

  registerWrapper( viewKeyStruct::hardeningString, &m_hardening )->
    setApplyDefaultValue( -1 )->
    setDescription( "Hardening parameter" );

  registerWrapper( viewKeyStruct::stateString, &m_state )->
    setApplyDefaultValue( 0.0 )->
    setPlotLevel( dataRepository::PlotLevel::LEVEL_0 )->
    setDescription( "Equivalent plastic shear strain" );
}

DruckerPragerExtended::~DruckerPragerExtended()
{}

void DruckerPragerExtended::allocateConstitutiveData( dataRepository::Group * const parent,
                                                      localIndex const numConstitutivePointsPerParentIndex )
{
  m_state.resize( 0, numConstitutivePointsPerParentIndex );

  ElasticIsotropic::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void DruckerPragerExtended::postProcessInput()
{
  ElasticIsotropic::postProcessInput();

  GEOSX_ASSERT_MSG( m_defaultCohesion >= 0, "Negative cohesion value detected" );
  GEOSX_ASSERT_MSG( m_defaultInitialFrictionAngle >= 0, "Negative initial friction angle detected" );
  GEOSX_ASSERT_MSG( m_defaultResidualFrictionAngle >= 0, "Negative residual friction angle detected" );
  GEOSX_ASSERT_MSG( m_defaultDilationRatio >= 0, "Dilation ratio out of [0,1] range detected" );
  GEOSX_ASSERT_MSG( m_defaultDilationRatio <= 1, "Dilation ratio out of [0,1] range detected" );
  GEOSX_ASSERT_MSG( m_defaultHardening >= 0, "Negative hardening parameter detected" );

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

  this->getWrapper< array1d< real64 > >( viewKeyStruct::initialFrictionString )->
    setApplyDefaultValue( F_i );
  this->getWrapper< array1d< real64 > >( viewKeyStruct::residualFrictionString )->
    setApplyDefaultValue( F_r );
  this->getWrapper< array1d< real64 > >( viewKeyStruct::pressureInterceptString )->
    setApplyDefaultValue( P );
  this->getWrapper< array1d< real64 > >( viewKeyStruct::dilationRatioString )->
    setApplyDefaultValue( m_defaultDilationRatio );
  this->getWrapper< array1d< real64 > >( viewKeyStruct::hardeningString )->
    setApplyDefaultValue( m_defaultHardening );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, DruckerPragerExtended, string const &, Group * const )
}
} /* namespace geosx */
