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
 *  @file JohnsonCook.cpp
 */

#include "JohnsonCook.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

JohnsonCook::JohnsonCook( string const & name, Group * const parent ) :
  ElasticIsotropic( name, parent ),
  
  m_hardeningModulus(),
  m_strainRateCoeff(),
  m_thermalSofteningExponent(),
  m_meltingTemperature(),
  m_referenceTemperature(),
  m_deformationGradient(),
  m_velocityGradient(),
  m_plasticStrain()
  
{
  registerWrapper( viewKeyStruct::defaultYieldStrengthString(), &m_defaultYieldStrength )
    .setInputFlag( InputFlags::REQUIRED )
    .setDescription( "Default yield strength parameter" );

  registerWrapper( viewKeyStruct::hardeningModulusString(), &m_hardeningModulus )
    .setInputFlag( InputFlags::REQUIRED )
    .setDescription( "Hardening modulus parameter" );

  registerWrapper( viewKeyStruct::strainRateCoeffString(), &m_strainRateCoeff )
    .setInputFlag( InputFlags::REQUIRED )
    .setDescription( "Strain rate sensitivity coefficient" );

  registerWrapper( viewKeyStruct::thermalSofteningExponentString(), &m_thermalSofteningExponent)
    .setInputFlag( InputFlags::REQUIRED )
    .setDescription( "Thermal softening coefficient" );

  registerWrapper( viewKeyStruct::referenceStrainRateString(), &m_referenceStrainRate )
    .setInputFlag( InputFlags::REQUIRED )
    .setDescription( "Reference strain rate" );

  registerWrapper( viewKeyStruct::meltingTemperatureString(), &m_meltingTemperature )
    .setInputFlag( InputFlags::REQUIRED )
    .setDescription( "Melting temperature" );

  registerWrapper( viewKeyStruct::referenceTemperatureString(), &m_referenceTemperature )
    .setInputFlag( InputFlags::REQUIRED )
    .setDescription( "Reference temperature" );

  registerWrapper( viewKeyStruct::yieldStrengthString(), &m_yieldStrength )
    .setApplyDefaultValue( -1 )
    .setDescription( "Yield strength field" );

  registerWrapper( viewKeyStruct::deformationGradientString(), &m_deformationGradient )
    .setApplyDefaultValue( 1.0 )
    .setPlotLevel( PlotLevel::NOPLOT )
    .setDescription( "Deformation gradient field" );

  registerWrapper( viewKeyStruct::velocityGradientString(), &m_velocityGradient )
    .setApplyDefaultValue( 0.0 )
    .setPlotLevel( PlotLevel::NOPLOT )
    .setDescription( "Velocity gradient field" );

  registerWrapper( viewKeyStruct::plasticStrainString(), &m_plasticStrain )
    .setApplyDefaultValue( 0.0 )
    .setPlotLevel( PlotLevel::NOPLOT )
    .setDescription( "Plastic strain field" );

}

JohnsonCook::~JohnsonCook()
{}

void JohnsonCook::allocateConstitutiveData( dataRepository::Group & parent,
                                            localIndex const numConstitutivePointsPerParentIndex )
{
  ElasticIsotropic::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_deformationGradient.resize( 0, 3, 3 );
  m_velocityGradient.resize( 0, 3, 3 );
  m_plasticStrain.resize( 0, numConstitutivePointsPerParentIndex, 6 );
  m_temperature.resize( 0, numConstitutivePointsPerParentIndex );
  m_yieldStrength.resize( 0 );
}

void JohnsonCook::postInputInitialization()
{
  ElasticIsotropic::postInputInitialization();

  this->getWrapper< array1d< real64 > >( viewKeyStruct::yieldStrengthString() )
    .setApplyDefaultValue( m_defaultInitialYieldStrength );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::temperatureString() )
    .setApplyDefaultValue( m_referenceTemperature );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, JohnsonCook, string const &, Group * const )

} // namespace constitutive
} // namespace geos


