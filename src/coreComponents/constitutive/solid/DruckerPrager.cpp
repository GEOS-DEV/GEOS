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
 *  @file DruckerPrager.cpp
 */

#include "DruckerPrager.hpp"

namespace geosx
{
using namespace dataRepository;
namespace constitutive
{

DruckerPrager::DruckerPrager( std::string const & name, Group * const parent ):
  ElasticIsotropic( name, parent ),
  m_defaultFrictionAngle(),
  m_defaultDilationAngle(),
  m_defaultCohesion(),
  m_defaultHardening(),
  m_friction(),
  m_dilation(),
  m_hardening(),
  m_cohesion()
{
  // register default values

  registerWrapper( viewKeyStruct::defaultFrictionAngleString, &m_defaultFrictionAngle )->
    setApplyDefaultValue( 30.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Friction angle (degrees)" );

  registerWrapper( viewKeyStruct::defaultDilationAngleString, &m_defaultDilationAngle )->
    setApplyDefaultValue( 30.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Dilation angle (degrees)" );

  registerWrapper( viewKeyStruct::defaultHardeningString, &m_defaultHardening )->
    setApplyDefaultValue( 0.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Cohesion hardening/softening rate" );

  registerWrapper( viewKeyStruct::defaultCohesionString, &m_defaultCohesion )->
    setApplyDefaultValue( 0.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Initial cohesion" );

  // register fields

  registerWrapper( viewKeyStruct::frictionString, &m_friction )->
    setApplyDefaultValue( -1 )->
    setDescription( "Yield surface slope" );

  registerWrapper( viewKeyStruct::dilationString, &m_dilation )->
    setApplyDefaultValue( -1 )->
    setDescription( "Plastic potential slope" );

  registerWrapper( viewKeyStruct::hardeningString, &m_hardening )->
    setApplyDefaultValue( -1 )->
    setDescription( "Hardening rate" );

  registerWrapper( viewKeyStruct::cohesionString, &m_cohesion )->
    setApplyDefaultValue( -1 )->
    setPlotLevel( dataRepository::PlotLevel::LEVEL_3 )->
    setDescription( "Current cohesion" );
}

DruckerPrager::~DruckerPrager()
{}

void DruckerPrager::allocateConstitutiveData( dataRepository::Group * const parent,
                                              localIndex const numConstitutivePointsPerParentIndex )
{
  m_cohesion.resize( 0, numConstitutivePointsPerParentIndex );

  ElasticIsotropic::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}


void DruckerPrager::postProcessInput()
{
  ElasticIsotropic::postProcessInput();

  GEOSX_ASSERT_MSG( m_defaultCohesion >= 0, "Negative cohesion value detected" );
  GEOSX_ASSERT_MSG( m_defaultFrictionAngle >= 0, "Negative friction angle detected" );
  GEOSX_ASSERT_MSG( m_defaultDilationAngle >= 0, "Negative dilation angle detected" );
  GEOSX_ASSERT_MSG( m_defaultFrictionAngle >= m_defaultDilationAngle, "Dilation angle should not exceed friction angle" );

  // convert from Mohr-Coulomb constants to Drucker-Prager constants, assuming DP
  // passes through the triaxial compression corners of the MC surface.
  // see Borja (2013) p. 75

  real64 phi = m_defaultFrictionAngle * M_PI / 180;
  real64 psi = m_defaultDilationAngle * M_PI / 180;

  real64 C = 6 * m_defaultCohesion * cos( phi ) / ( 3 - sin( phi ) );
  real64 F = 6 * sin( phi ) / ( 3 - sin( phi ) );
  real64 D = 6 * sin( psi ) / ( 3 - sin( psi ) );

  // set results as array default values

  this->getWrapper< array2d< real64 > >( viewKeyStruct::cohesionString )->
    setApplyDefaultValue( C );
  this->getWrapper< array1d< real64 > >( viewKeyStruct::dilationString )->
    setApplyDefaultValue( D );
  this->getWrapper< array1d< real64 > >( viewKeyStruct::frictionString )->
    setApplyDefaultValue( F );
  this->getWrapper< array1d< real64 > >( viewKeyStruct::hardeningString )->
    setApplyDefaultValue( m_defaultHardening );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, DruckerPrager, std::string const &, Group * const )
}
} /* namespace geosx */

