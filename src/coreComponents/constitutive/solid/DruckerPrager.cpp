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


void
DruckerPrager::DeliverClone( string const & name,
                                      Group * const parent,
                                      std::unique_ptr< ConstitutiveBase > & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique< DruckerPrager >( name, parent );
  }
  ElasticIsotropic::DeliverClone( name, parent, clone );
  
  DruckerPrager * const newConstitutiveRelation = dynamic_cast< DruckerPrager * >(clone.get());

  newConstitutiveRelation->m_defaultTanFrictionAngle = m_defaultTanFrictionAngle;
  newConstitutiveRelation->m_defaultTanDilationAngle = m_defaultTanDilationAngle;
  newConstitutiveRelation->m_defaultCohesion         = m_defaultCohesion;
  newConstitutiveRelation->m_defaultHardeningRate    = m_defaultHardeningRate;
  
  newConstitutiveRelation->m_tanFrictionAngle = m_tanFrictionAngle;
  newConstitutiveRelation->m_tanDilationAngle = m_tanDilationAngle;
  newConstitutiveRelation->m_hardeningRate    = m_hardeningRate;
  newConstitutiveRelation->m_newCohesion      = m_newCohesion;
  newConstitutiveRelation->m_oldCohesion      = m_oldCohesion;
}

void DruckerPrager::AllocateConstitutiveData( dataRepository::Group * const parent,
                                              localIndex const numConstitutivePointsPerParentIndex )
{
  ElasticIsotropic::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
  
  localIndex const numElems = parent->size();
  this->resize( numElems );
  
  // 1d arrays
  m_tanFrictionAngle.resize( numElems );
  m_tanDilationAngle.resize( numElems );
  m_hardeningRate.resize( numElems );
  
  // 2d arrays
  m_newCohesion.resize( numElems, numConstitutivePointsPerParentIndex );
  m_oldCohesion.resize( numElems, numConstitutivePointsPerParentIndex );
}

void DruckerPrager::PostProcessInput()
{
  GEOSX_ASSERT_MSG(m_defaultCohesion >= 0, "Negative cohesion value detected");
  GEOSX_ASSERT_MSG(m_defaultTanFrictionAngle >= 0, "Negative friction angle detected");
  GEOSX_ASSERT_MSG(m_defaultTanDilationAngle >= 0, "Negative dilation angle detected");
  GEOSX_ASSERT_MSG(m_defaultTanFrictionAngle >= m_defaultTanDilationAngle, "Friction angle should exceed dilation angle");
  
  // set results as array default values
  this->getWrapper< array1d< real64 > >( viewKeyStruct::oldCohesionString )->
    setApplyDefaultValue( m_defaultCohesion );
  this->getWrapper< array1d< real64 > >( viewKeyStruct::newCohesionString )->
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
