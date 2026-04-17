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
 * @file Damage.cpp
 */

#include "Damage.hpp"
#include "SolidFields.hpp"

namespace geos
{

using namespace dataRepository;
namespace constitutive
{

template< typename BASE >
Damage< BASE >::Damage( string const & name, Group * const parent ):
  BASE( name, parent )
{
  this->registerWrapper( viewKeyStruct::lengthScaleString(), &m_lengthScale ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Length scale l in the phase-field equation" );

  this->registerWrapper( viewKeyStruct::defaultCriticalFractureEnergyString(), &m_defaultCriticalFractureEnergy ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default critical fracture energy" );

  this->registerWrapper( viewKeyStruct::criticalStrainEnergyString(), &m_criticalStrainEnergy ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Critical stress in a 1d tension test" );

  this->registerWrapper( viewKeyStruct::degradationLowerLimitString(), &m_degradationLowerLimit ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The lower limit of the degradation function" );

  this->registerWrapper( viewKeyStruct::extDrivingForceFlagString(), &m_extDrivingForceFlag ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Whether to have external driving force. Can be 0 or 1" );

  this->registerWrapper( viewKeyStruct::defaultTensileStrengthString(), &m_defaultTensileStrength ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default tensile strength from the uniaxial tension test" );

  this->registerWrapper( viewKeyStruct::defaultCompressiveStrengthString(), &m_defaultCompressiveStrength ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default compressive strength from the uniaxial compression test" );

  this->registerWrapper( viewKeyStruct::defaultDeltaCoefficientString(), &m_defaultDeltaCoefficient ).
    setApplyDefaultValue( -1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default coefficient in the calculation of the external driving force" );

  // fields

  this->template registerField< fields::solid::damage >( &m_newDamage );

  this->template registerField< fields::solid::oldDamage >( &m_oldDamage );

  this->template registerField< fields::solid::damageGrad >( &m_damageGrad );

  this->template registerField< fields::solid::strainEnergyDensity >( &m_strainEnergyDensity );

  this->template registerField< fields::solid::volStrain >( &m_volStrain );

  this->template registerField< fields::solid::extDrivingForce >( &m_extDrivingForce );

  this->template registerField< fields::solid::criticalFractureEnergy >( &m_criticalFractureEnergy );

  this->template registerField< fields::solid::tensileStrength >( &m_tensileStrength );

  this->template registerField< fields::solid::compressiveStrength >( &m_compressiveStrength );

  this->template registerField< fields::solid::deltaCoefficient >( &m_deltaCoefficient );

  this->template registerField< fields::solid::biotCoefficient >( &m_biotCoefficient );
}


template< typename BASE >
void Damage< BASE >::postInputInitialization()
{
  BASE::postInputInitialization();

  GEOS_ERROR_IF( m_extDrivingForceFlag != 0 && m_extDrivingForceFlag!= 1,
                 "invalid external driving force flag option - must"
                 " be 0 or 1",
                 this->getDataContext() );
  GEOS_ERROR_IF( m_extDrivingForceFlag == 1 && m_defaultTensileStrength <= 0.0,
                 "tensile strength must be input and positive when the"
                 " external driving force flag is turned on",
                 this->getDataContext()  );
  GEOS_ERROR_IF( m_extDrivingForceFlag == 1 && m_defaultCompressiveStrength  <= 0.0,
                 "compressive strength must be input and positive when the"
                 " external driving force flag is turned on",
                 this->getDataContext()  );
  GEOS_ERROR_IF( m_extDrivingForceFlag == 1 && m_defaultDeltaCoefficient < 0.0,
                 "delta coefficient must be input and non-negative when the"
                 " external driving force flag is turned on",
                 this->getDataContext()  );

  // set results as array default values
  this->template getField< fields::solid::criticalFractureEnergy >().
    setApplyDefaultValue( m_defaultCriticalFractureEnergy );

  this->template getField< fields::solid::tensileStrength >().
    setApplyDefaultValue( m_defaultTensileStrength );

  this->template getField< fields::solid::compressiveStrength >().
    setApplyDefaultValue( m_defaultCompressiveStrength );

  this->template getField< fields::solid::deltaCoefficient >().
    setApplyDefaultValue( m_defaultDeltaCoefficient );
}

template< typename BASE >
void Damage< BASE >::allocateConstitutiveData( Group & parent, localIndex const numPts )
{
  m_newDamage.resize( 0, numPts );
  m_oldDamage.resize( 0, numPts );
  m_damageGrad.resize( 0, numPts, 3 );
  m_strainEnergyDensity.resize( 0, numPts );
  m_volStrain.resize( 0, numPts );
  m_extDrivingForce.resize( 0, numPts );
  m_biotCoefficient.resize( parent.size() );
  m_criticalFractureEnergy.resize( parent.size() );
  m_tensileStrength.resize( parent.size() );
  m_compressiveStrength.resize( parent.size() );
  m_deltaCoefficient.resize( parent.size() );

  BASE::allocateConstitutiveData( parent, numPts );
}

template< typename BASE >
void Damage< BASE >::saveConvergedState() const
{
  SolidBase::saveConvergedState(); // TODO: not ideal, as we have separate loops for base and derived data

  localIndex const numE = SolidBase::numElem();
  localIndex const numQ = SolidBase::numQuad();

  arrayView2d< real64 const > newDamage = m_newDamage;
  arrayView2d< real64 > oldDamage = m_oldDamage;

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < numQ; ++q )
    {
      oldDamage( k, q ) = newDamage( k, q );
    }
  } );
}

typedef Damage< ElasticIsotropic > DamageElasticIsotropic;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, DamageElasticIsotropic, string const &, Group * const )

}
} /* namespace geos */
