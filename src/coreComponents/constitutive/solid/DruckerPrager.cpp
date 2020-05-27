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
  SolidBase( name, parent ),
  m_defaultBulkModulus(),
  m_defaultPoissonRatio(),
  m_defaultTanFrictionAngle(),
  m_defaultCohesion(),
  m_defaultHardeningRate(),
  m_bulkModulus(),
  m_poissonRatio(),
  m_tanFrictionAngle(),
  m_hardeningRate(),
  m_newCohesion(),
  m_oldCohesion(),
  m_newStress(),
  m_oldStress()
{
  // register default values
  
  registerWrapper( viewKeyStruct::defaultBulkModulusString, &m_defaultBulkModulus )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Elastic bulk modulus parameter" );

  registerWrapper( viewKeyStruct::defaultPoissonRatioString, &m_defaultPoissonRatio )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Elastic Poisson ratio parameter" );
  
  registerWrapper( viewKeyStruct::defaultTanFrictionAngleString, &m_defaultTanFrictionAngle )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Yield surface slope parameter tan(phi)" );
  
  registerWrapper( viewKeyStruct::defaultHardeningRateString, &m_defaultHardeningRate )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Cohesion hardening/softening rate parameter" );
  
  registerWrapper( viewKeyStruct::defaultCohesionString, &m_defaultCohesion )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Initial cohesion parameter" );

  // register fields
  
  registerWrapper( viewKeyStruct::bulkModulusString, &m_bulkModulus )->
    setApplyDefaultValue( -1 )->
    setDescription( "Elastic bulk modulus field" );

  registerWrapper( viewKeyStruct::poissonRatioString, &m_poissonRatio )->
    setApplyDefaultValue( -1 )->
    setDescription( "Elastic Poisson ration field" );
  
  registerWrapper( viewKeyStruct::tanFrictionAngleString, &m_tanFrictionAngle )->
    setApplyDefaultValue( -1 )->
    setDescription( "Yield surface slope tan(phi) field" );
  
  registerWrapper( viewKeyStruct::hardeningRateString, &m_hardeningRate )->
    setApplyDefaultValue( -1 )->
    setDescription( "Hardening rate field" );
  
  registerWrapper( viewKeyStruct::newCohesionString, &m_newCohesion )->
    setApplyDefaultValue( -1 )->
    setDescription( "New cohesion field" );
  
  registerWrapper( viewKeyStruct::oldCohesionString, &m_newCohesion )->
    setApplyDefaultValue( -1 )->
    setDescription( "Old cohesion field" );
  
  registerWrapper( viewKeyStruct::newStressString, &m_newStress )->
    setApplyDefaultValue( -1 )->
    setDescription( "New stress field" );
  
  registerWrapper( viewKeyStruct::oldStressString, &m_oldStress )->
    setApplyDefaultValue( -1 )->
    setDescription( "Old stress field" );
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
  SolidBase::DeliverClone( name, parent, clone );
  DruckerPrager * const newConstitutiveRelation = dynamic_cast< DruckerPrager * >(clone.get());

  newConstitutiveRelation->m_defaultBulkModulus      = m_defaultBulkModulus;
  newConstitutiveRelation->m_defaultPoissonRatio     = m_defaultPoissonRatio;
  newConstitutiveRelation->m_defaultTanFrictionAngle = m_defaultTanFrictionAngle;
  newConstitutiveRelation->m_defaultCohesion         = m_defaultCohesion;
  newConstitutiveRelation->m_defaultHardeningRate    = m_defaultHardeningRate;
  
  newConstitutiveRelation->m_bulkModulus = m_bulkModulus;
  newConstitutiveRelation->m_poissonRatio = m_poissonRatio;
  newConstitutiveRelation->m_tanFrictionAngle = m_tanFrictionAngle;
  newConstitutiveRelation->m_hardeningRate = m_hardeningRate;
  newConstitutiveRelation->m_newCohesion = m_newCohesion;
  newConstitutiveRelation->m_oldCohesion = m_oldCohesion;
  newConstitutiveRelation->m_newStress = m_newStress;
  newConstitutiveRelation->m_oldStress = m_oldStress;
}

void DruckerPrager::AllocateConstitutiveData( dataRepository::Group * const parent,
                                                       localIndex const numConstitutivePointsPerParentIndex )
{
  SolidBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
  this->resize( parent->size() );
  
  m_bulkModulus.resize( parent->size() );
  m_poissonRatio.resize( parent->size() );
  m_tanFrictionAngle.resize( parent->size() );
  m_hardeningRate.resize( parent->size() );
  
  m_newCohesion.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_oldCohesion.resize( parent->size(), numConstitutivePointsPerParentIndex );
  
  m_newStress.resize( parent->size(), numConstitutivePointsPerParentIndex, 6 );
  m_oldStress.resize( parent->size(), numConstitutivePointsPerParentIndex, 6 ); // TODO: figure out how to set initial stress
  
  // set arrays to default values
  m_bulkModulus = m_defaultBulkModulus;
  m_poissonRatio = m_defaultPoissonRatio;
  m_tanFrictionAngle = m_defaultTanFrictionAngle;
  m_hardeningRate = m_defaultHardeningRate;
  m_newCohesion = m_defaultCohesion;
  m_oldCohesion = m_defaultCohesion;
}

void DruckerPrager::PostProcessInput()
{
  m_postProcessed = true; // TODO: add parameter conversion helper class for more flexible input
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, DruckerPrager, std::string const &, Group * const )
}
} /* namespace geosx */
