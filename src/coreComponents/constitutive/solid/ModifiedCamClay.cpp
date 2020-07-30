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
 *  @file ModifiedCamClay.cpp
 */

#include "ModifiedCamClay.hpp"

namespace geosx
{
using namespace dataRepository;
namespace constitutive
{

  ModifiedCamClay::ModifiedCamClay( std::string const & name, Group * const parent ):
  SolidBase( name, parent ),
  m_defaultRefPInvariant(),
  m_defaultRefElasticStrainVolumetric(),
  m_defaultRefShearModulus(),
  m_defaultShearModulusEvolution(),
  m_defaultVirginCompressionIndex(),
  m_defaultRecompressionIndex(),
  m_defaultCriticalStateSlope(),
  m_defaultAssociativity(),
  m_defaultPreconsolidationPressure(),
  m_referencePInvariant(),
  m_refElasticStrainVolumetric(),
  m_referenceShearModulus(),
  m_shearModulusEvolution(),
  m_virginCompressionIndex(),
  m_recompressionIndex(),
  m_criticalStateSlope(),
  m_associativity(),
  m_newPreconsolidationPressure(),
  m_oldPreconsolidationPressure(),
  m_newElasticStrain(),
  m_oldElasticStrain()
{
  // register default values
  
  registerWrapper( viewKeyStruct::defaultRefPInvariantString, &m_defaultRefPInvariant )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Reference stress invariant P" ); //Q: set default value

  registerWrapper( viewKeyStruct::defaultRefElasticStrainVolumetricString, &m_defaultRefElasticStrainVolumetric )->
    setApplyDefaultValue( 0.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Reference volumetric invariant of the elastic strain" );
  
  registerWrapper( viewKeyStruct::defaultRefShearModulusString, &m_defaultRefShearModulus )->
    setApplyDefaultValue( 1.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Reference shear modulus parameter" );
  
  registerWrapper( viewKeyStruct::defaultShearModulusEvolutionString, &m_defaultShearModulusEvolution )->
    setApplyDefaultValue( 0.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Parameter responsable for shear modulus evolution" );
  
  registerWrapper( viewKeyStruct::defaultVirginCompressionIndexString, &m_defaultVirginCompressionIndex )->
    setApplyDefaultValue( 0.1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Virgin compression index parameter (Cc)" );
  
  registerWrapper( viewKeyStruct::defaultRecompressionIndexString, &m_defaultRecompressionIndex )->
    setApplyDefaultValue( 0.01 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Recompression index parameter (Cr)" );

  registerWrapper( viewKeyStruct::defaultCriticalStateSlopeString, &m_defaultCriticalStateSlope )->
    setApplyDefaultValue( 1.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Slope of the critical state line (M)" );

  registerWrapper( viewKeyStruct::defaultAssociativityString, &m_defaultAssociativity )->
    setApplyDefaultValue( 1.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Associativity parameter" );

  registerWrapper( viewKeyStruct::defaultPreconsolidationPressureString, &m_defaultPreconsolidationPressure )->
    setApplyDefaultValue( -1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Initial preconsolidation pressure" ); //Q: set default value

  // register fields
  
  registerWrapper( viewKeyStruct::refPInvariantString, &m_referencePInvariant )->
    setApplyDefaultValue( -1 )->
    setDescription( "Reference stress invariant P field" );

  registerWrapper( viewKeyStruct::refElasticStrainVolumetricString, &m_refElasticStrainVolumetric )->
    setApplyDefaultValue( -1 )->
    setDescription( "Reference volumetric invariant of the elastic strain field" );
  
  registerWrapper( viewKeyStruct::refShearModulusString, &m_referenceShearModulus )->
    setApplyDefaultValue( -1 )->
    setDescription( "Reference shear modulus field" );
  
  registerWrapper( viewKeyStruct::shearModulusEvolutionString, &m_shearModulusEvolution )->
    setApplyDefaultValue( -1 )->
    setDescription( "Shear modulus evolution field" );
  
  registerWrapper( viewKeyStruct::virginCompressionIndexString, &m_virginCompressionIndex )->
    setApplyDefaultValue( -1 )->
    setDescription( "Virgin compression index field" );
  
  registerWrapper( viewKeyStruct::recompressionIndexString, &m_recompressionIndex )->
    setApplyDefaultValue( -1 )->
    setDescription( "Recompression index field" );
  
  registerWrapper( viewKeyStruct::criticalStateSlopeString, &m_criticalStateSlope )->
    setApplyDefaultValue( -1 )->
    setDescription( "Slope of the critical state line field" );

  registerWrapper( viewKeyStruct::associativityString, &m_associativity )->
    setApplyDefaultValue( -1 )->
    setDescription( "Associativity field" );
  
  registerWrapper( viewKeyStruct::newPreconsolidationPressureString, &m_newPreconsolidationPressure )->
    setApplyDefaultValue( -1 )->
    setDescription( "New preconsolidation pressure field" );
  
  registerWrapper( viewKeyStruct::oldPreconsolidationPressureString, &m_oldPreconsolidationPressure )->
    setApplyDefaultValue( -1 )->
    setDescription( "Old preconsolidation pressure field" );

  registerWrapper( viewKeyStruct::newElasticStrainString, &m_newElasticStrain )->
    setApplyDefaultValue( -1 )->
    setDescription( "New elastic strain field" );

  registerWrapper( viewKeyStruct::oldElasticStrainString, &m_oldElasticStrain )->
    setApplyDefaultValue( -1 )->
    setDescription( "Old elastic strain field" );
}


ModifiedCamClay::~ModifiedCamClay()
{}


void
ModifiedCamClay::DeliverClone( string const & name,
                                      Group * const parent,
                                      std::unique_ptr< ConstitutiveBase > & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique< ModifiedCamClay >( name, parent );
  }
  SolidBase::DeliverClone( name, parent, clone );
  ModifiedCamClay * const newConstitutiveRelation = dynamic_cast< ModifiedCamClay * >(clone.get());

  newConstitutiveRelation->m_defaultRefPInvariant = m_defaultRefPInvariant;
  newConstitutiveRelation->m_defaultRefElasticStrainVolumetric = m_defaultRefElasticStrainVolumetric;
  newConstitutiveRelation->m_defaultRefShearModulus = m_defaultRefShearModulus;
  newConstitutiveRelation->m_defaultShearModulusEvolution = m_defaultShearModulusEvolution;
  newConstitutiveRelation->m_defaultVirginCompressionIndex = m_defaultVirginCompressionIndex;
  newConstitutiveRelation->m_defaultRecompressionIndex = m_defaultRecompressionIndex;
  newConstitutiveRelation->m_defaultCriticalStateSlope = m_defaultCriticalStateSlope;
  newConstitutiveRelation->m_defaultAssociativity = m_defaultAssociativity;
  newConstitutiveRelation->m_defaultPreconsolidationPressure = m_defaultPreconsolidationPressure;

  newConstitutiveRelation->m_referencePInvariant = m_referencePInvariant;
  newConstitutiveRelation->m_refElasticStrainVolumetric = m_refElasticStrainVolumetric;
  newConstitutiveRelation->m_referenceShearModulus = m_referenceShearModulus;
  newConstitutiveRelation->m_shearModulusEvolution = m_shearModulusEvolution;
  newConstitutiveRelation->m_virginCompressionIndex = m_virginCompressionIndex;
  newConstitutiveRelation->m_recompressionIndex = m_recompressionIndex;
  newConstitutiveRelation->m_criticalStateSlope = m_criticalStateSlope;
  newConstitutiveRelation->m_associativity = m_associativity;
  newConstitutiveRelation->m_newPreconsolidationPressure = m_newPreconsolidationPressure;
  newConstitutiveRelation->m_oldPreconsolidationPressure = m_oldPreconsolidationPressure;
  newConstitutiveRelation->m_newElasticStrain = m_newElasticStrain;
  newConstitutiveRelation->m_oldElasticStrain = m_oldElasticStrain;
}

void ModifiedCamClay::AllocateConstitutiveData( dataRepository::Group * const parent,
                                                       localIndex const numConstitutivePointsPerParentIndex )
{
  SolidBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
  this->resize( parent->size() );
  
  m_referencePInvariant.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_refElasticStrainVolumetric.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_referenceShearModulus.resize( parent->size(), numConstitutivePointsPerParentIndex );

  m_shearModulusEvolution.resize( parent->size() );
  m_virginCompressionIndex.resize( parent->size() );
  m_recompressionIndex.resize( parent->size() );
  m_criticalStateSlope.resize( parent->size() );
  m_associativity.resize( parent->size() );

  m_newPreconsolidationPressure.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_oldPreconsolidationPressure.resize( parent->size(), numConstitutivePointsPerParentIndex );
  
  m_newElasticStrain.resize( parent->size(), numConstitutivePointsPerParentIndex, 6 );
  m_oldElasticStrain.resize( parent->size(), numConstitutivePointsPerParentIndex, 6 );// TODO: figure out how to set initial strain
  
  // set arrays to default values
  m_referencePInvariant = m_defaultRefPInvariant;
  m_refElasticStrainVolumetric = m_defaultRefElasticStrainVolumetric;
  m_referenceShearModulus = m_defaultRefShearModulus;
  m_shearModulusEvolution = m_defaultShearModulusEvolution;
  m_virginCompressionIndex = m_defaultVirginCompressionIndex;
  m_recompressionIndex = m_defaultRecompressionIndex;
  m_criticalStateSlope = m_defaultCriticalStateSlope;
  m_associativity = m_defaultAssociativity;
  m_newPreconsolidationPressure = m_defaultPreconsolidationPressure;
  m_oldPreconsolidationPressure = m_defaultPreconsolidationPressure;
}

void ModifiedCamClay::PostProcessInput()
{
  GEOSX_ASSERT_MSG(m_defaultShearModulusEvolution >= 0, "Negative shear modulus evolution parameter detected");
  GEOSX_ASSERT_MSG(m_defaultVirginCompressionIndex >= 0, "Negative virgin compression index detected");
  GEOSX_ASSERT_MSG(m_defaultRecompressionIndex >= 0, "Negative recompression index detected");
  GEOSX_ASSERT_MSG(m_defaultCriticalStateSlope >= 0, "Negative critical state line slope detected");

  // Check this assessment of negative Poisson ratio
  // Q: should we assess with the default values?

  // real64 initialP = m_defaultRefPInvariant;
  // real64 Cr = m_defaultRecompressionIndex;
  // real64 shear = m_defaultRefShearModulus;
  // real64 bulk = -initialP / Cr;
  // real64 poissonRatio = ( 3 * bulk - 2 * shear) / ( 2 * ( 3 * bulk + shear ) );

  // GEOSX_ERROR_IF( poissonRatio < 0 , "Negative poisson ratio produced" );
  
  m_postProcessed = true; // TODO: add parameter conversion helper class for more flexible input
} 

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ModifiedCamClay, std::string const &, Group * const )
}
} /* namespace geosx */
