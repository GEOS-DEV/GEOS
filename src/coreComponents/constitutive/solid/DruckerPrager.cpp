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
  m_defaultTanFrictionAngle(),
  m_defaultTanDilationAngle(),
  m_defaultCohesion(),
  m_defaultHardeningRate(),
  m_tanFrictionAngle(),
  m_tanDilationAngle(),
  m_hardeningRate(),
  m_newCohesion(),
  m_oldCohesion()
{
  // register default values

  registerWrapper( viewKeyStruct::defaultTanFrictionAngleString, &m_defaultTanFrictionAngle )->
    setApplyDefaultValue( 1.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Yield surface slope parameter tan(phi)" );
  
  registerWrapper( viewKeyStruct::defaultTanDilationAngleString, &m_defaultTanDilationAngle )->
    setApplyDefaultValue( 0.5 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Plastic potential slope parameter tan(psi)" );
  
  registerWrapper( viewKeyStruct::defaultHardeningRateString, &m_defaultHardeningRate )->
    setApplyDefaultValue( 1e8 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Cohesion hardening/softening rate parameter" );
  
  registerWrapper( viewKeyStruct::defaultCohesionString, &m_defaultCohesion )->
    setApplyDefaultValue( 5e6 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Initial cohesion parameter" );

  // register fields
  
  registerWrapper( viewKeyStruct::tanFrictionAngleString, &m_tanFrictionAngle )->
    setApplyDefaultValue( -1 )->
    setDescription( "Yield surface slope tan(phi) field" );
  
  registerWrapper( viewKeyStruct::tanDilationAngleString, &m_tanDilationAngle )->
    setApplyDefaultValue( -1 )->
    setDescription( "Plastic potential slope tan(psi) field" );
  
  registerWrapper( viewKeyStruct::hardeningRateString, &m_hardeningRate )->
    setApplyDefaultValue( -1 )->
    setDescription( "Hardening rate field" );
  
  registerWrapper( viewKeyStruct::newCohesionString, &m_newCohesion )->
    setApplyDefaultValue( -1 )->
    setDescription( "New cohesion field" );
  
  registerWrapper( viewKeyStruct::oldCohesionString, &m_newCohesion )->
    setApplyDefaultValue( -1 )->
    setDescription( "Old cohesion field" );
}


DruckerPrager::~DruckerPrager()
{}


void DruckerPrager::allocateConstitutiveData( dataRepository::Group * const parent,
                                              localIndex const numConstitutivePointsPerParentIndex )
{
  m_newCohesion.resize( 0, numConstitutivePointsPerParentIndex );
  m_oldCohesion.resize( 0, numConstitutivePointsPerParentIndex );
  
  ElasticIsotropic::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void DruckerPrager::PostProcessInput()
{
  GEOSX_ASSERT_MSG(m_defaultCohesion >= 0, "Negative cohesion value detected");
  GEOSX_ASSERT_MSG(m_defaultTanFrictionAngle >= 0, "Negative friction angle detected");
  GEOSX_ASSERT_MSG(m_defaultTanDilationAngle >= 0, "Negative dilation angle detected");
  GEOSX_ASSERT_MSG(m_defaultTanFrictionAngle >= m_defaultTanDilationAngle, "Friction angle should exceed dilation angle");
  
  // set results as array default values
  this->getWrapper< array2d< real64 > >( viewKeyStruct::oldCohesionString )->
    setApplyDefaultValue( m_defaultCohesion );
  this->getWrapper< array2d< real64 > >( viewKeyStruct::newCohesionString )->
    setApplyDefaultValue( m_defaultCohesion );
  this->getWrapper< array1d< real64 > >( viewKeyStruct::tanDilationAngleString )->
    setApplyDefaultValue( m_defaultTanDilationAngle );
  this->getWrapper< array1d< real64 > >( viewKeyStruct::tanFrictionAngleString )->
    setApplyDefaultValue( m_defaultTanFrictionAngle );
  this->getWrapper< array1d< real64 > >( viewKeyStruct::hardeningRateString )->
    setApplyDefaultValue( m_defaultHardeningRate );
  
  m_postProcessed = true; // TODO: add parameter conversion helper class for more flexible input
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, DruckerPrager, std::string const &, Group * const )
}
} /* namespace geosx */
