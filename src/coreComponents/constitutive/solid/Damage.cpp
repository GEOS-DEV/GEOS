
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
 * @file Damage.cpp
 */

#include "Damage.hpp"
#include "ElasticIsotropic.hpp"

namespace geosx
{

using namespace dataRepository;
namespace constitutive
{

template< typename BASE >
Damage< BASE >::Damage( string const & name, Group * const parent ):
  BASE( name, parent ),
  m_damage(),
  m_strainEnergyDensity(),
  m_extDrivingForce(),
  m_lengthScale(),
  m_criticalFractureEnergy(),
  m_criticalStrainEnergy(),
  m_degradationLowerLimit( 0.0 ),
  m_extDrivingForceFlag( 0 ),
  m_tensileStrength(),
  m_compressStrength(),
  m_deltaCoefficient()
{
  this->registerWrapper( viewKeyStruct::damageString(), &m_damage ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Material Damage Variable" );

  this->registerWrapper( viewKeyStruct::strainEnergyDensityString(), &m_strainEnergyDensity ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Strain Energy Density" );

  this->registerWrapper( viewKeyStruct::extDrivingForceString(), &m_extDrivingForce ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "External Driving Force" );

  this->registerWrapper( viewKeyStruct::lengthScaleString(), &m_lengthScale ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Length scale l in the phase-field equation" );

  this->registerWrapper( viewKeyStruct::criticalFractureEnergyString(), &m_criticalFractureEnergy ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Critical fracture energy" );

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

  this->registerWrapper( viewKeyStruct::tensileStrengthString(), &m_tensileStrength ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Tensile strength from the uniaxial tension test" );

  this->registerWrapper( viewKeyStruct::compressStrengthString(), &m_compressStrength ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Compressive strength from the uniaxial compression test" );

  this->registerWrapper( viewKeyStruct::deltaCoefficientString(), &m_deltaCoefficient ).
    setApplyDefaultValue( -1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Coefficient in the calculation of the external driving force" );
}


template< typename BASE >
void Damage< BASE >::postProcessInput()
{
  BASE::postProcessInput();

  GEOSX_ERROR_IF( m_extDrivingForceFlag != 0 && m_extDrivingForceFlag!= 1, "invalid external driving force flag option - must be 0 or 1" );
  GEOSX_ERROR_IF( m_extDrivingForceFlag == 1 && m_tensileStrength <= 0.0, "tensile strength must be input and positive when the external driving force flag is turned on" );
  GEOSX_ERROR_IF( m_extDrivingForceFlag == 1 && m_compressStrength <= 0.0, "compressive strength must be input and positive when the external driving force flag is turned on" );
  GEOSX_ERROR_IF( m_extDrivingForceFlag == 1 && m_deltaCoefficient < 0.0, "delta coefficient must be input and non-negative when the external driving force flag is turned on" );
}

template< typename BASE >
void Damage< BASE >::allocateConstitutiveData( dataRepository::Group & parent,
                                               localIndex const numConstitutivePointsPerParentIndex )
{
  m_damage.resize( 0, numConstitutivePointsPerParentIndex );
  m_strainEnergyDensity.resize( 0, numConstitutivePointsPerParentIndex );
  m_extDrivingForce.resize( 0, numConstitutivePointsPerParentIndex );
  BASE::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

typedef Damage< ElasticIsotropic > DamageElasticIsotropic;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, DamageElasticIsotropic, string const &, Group * const )

}
} /* namespace geosx */
