/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 *  @file StrainHardeningPolymer.cpp
 */

#include "StrainHardeningPolymer.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

StrainHardeningPolymer::StrainHardeningPolymer( string const & name, Group * const parent ):
  ElasticIsotropic( name, parent ),
  m_deformationGradient(),
  m_plasticStrain(),
  m_damage(),
  m_jacobian(),
  m_bulkModulusA( 0.0 ),
  m_bulkModulusB(0.0),
  m_bulkModulusC(0.0),
  m_bulkModulusD(0.0),
  m_shearModulusA(0.0),
  m_shearModulusB(0.0),
  m_shearModulusC(0.0),
  m_shearModulusD(0.0),
  m_yieldStrength(0.0),
  m_yieldStrengthB(0.0),
  m_yieldStrengthC(0.0),
  m_yieldStrengthD(0.0),
  m_strainHardeningSlope(0.0),
  m_strainHardeningSlopeB(0.0),
  m_strainHardeningSlopeC(0.0),
  m_strainHardeningSlopeD(0.0),
  m_shearSofteningMagnitude(0.0),
  m_shearSofteningMagnitudeB(0.0),
  m_shearSofteningMagnitudeC(0.0),
  m_shearSofteningMagnitudeD(0.0),
  m_shearSofteningShapeParameter1(0.0),
  m_shearSofteningShapeParameter2(0.0),
  m_initialTemperature( 0.0 ),
  m_temperature(0.0),
  m_maximumStretch(0.0),
  m_maximumStretchB(0.0),
  m_maximumStretchC(0.0),
  m_maximumStretchD(0.0)
{
  // register default values

  registerWrapper( viewKeyStruct::bulkModulusAString(), &m_bulkModulusA ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "bulk modulus  parameterA" );

  registerWrapper( viewKeyStruct::bulkModulusBString(), &m_bulkModulusB ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "bulk modulus  parameterB" );

  registerWrapper( viewKeyStruct::bulkModulusCString(), &m_bulkModulusC ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "bulk modulus  parameterC" );

  registerWrapper( viewKeyStruct::bulkModulusDString(), &m_bulkModulusD ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "bulk modulus  parameterD" );

  registerWrapper( viewKeyStruct::shearModulusAString(), &m_shearModulusA ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "shear modulus  parameterA" );

  registerWrapper( viewKeyStruct::shearModulusBString(), &m_shearModulusB ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "shear modulus  parameterB" );

  registerWrapper( viewKeyStruct::shearModulusCString(), &m_shearModulusC ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "shear modulus  parameterC" );

  registerWrapper( viewKeyStruct::shearModulusDString(), &m_shearModulusD ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "shear modulus  parameterD" );


  registerWrapper( viewKeyStruct::yieldStrengthString(), &m_yieldStrength ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "yield strength  parameterA" );

  registerWrapper( viewKeyStruct::yieldStrengthBString(), &m_yieldStrengthB ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "yield strength parameterB" );

  registerWrapper( viewKeyStruct::yieldStrengthCString(), &m_yieldStrengthC ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "yield strength parameterC" );

  registerWrapper( viewKeyStruct::yieldStrengthDString(), &m_yieldStrengthD ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "yield strength parameterD" );




  registerWrapper( viewKeyStruct::strainHardeningSlopeString(), &m_strainHardeningSlope ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Strain hardening slope" );

  registerWrapper( viewKeyStruct::strainHardeningSlopeBString(), &m_strainHardeningSlopeB ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Strain hardening slope" );

  registerWrapper( viewKeyStruct::strainHardeningSlopeCString(), &m_strainHardeningSlopeC ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Strain hardening slope" );

  registerWrapper( viewKeyStruct::strainHardeningSlopeDString(), &m_strainHardeningSlopeD ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Strain hardening slope" );

  registerWrapper( viewKeyStruct::shearSofteningMagnitudeString(), &m_shearSofteningMagnitude ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Shear softening magnitude" );

  registerWrapper( viewKeyStruct::shearSofteningMagnitudeBString(), &m_shearSofteningMagnitudeB ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Shear softening magnitude" );

  registerWrapper( viewKeyStruct::shearSofteningMagnitudeCString(), &m_shearSofteningMagnitudeC ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Shear softening magnitude" );

  registerWrapper( viewKeyStruct::shearSofteningMagnitudeDString(), &m_shearSofteningMagnitudeD ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Shear softening magnitude" );

  registerWrapper( viewKeyStruct::shearSofteningShapeParameter1String(), &m_shearSofteningShapeParameter1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Shear softening shape parameter 1" );

  registerWrapper( viewKeyStruct::shearSofteningShapeParameter2String(), &m_shearSofteningShapeParameter2 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Shear softening shape parameter 2" );

  registerWrapper( viewKeyStruct::initialTemperatureString(), &m_initialTemperature ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "initial Temperature" );

  registerWrapper( viewKeyStruct::temperatureString(), &m_temperature ).
    setApplyDefaultValue( 300.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Array of quadrature point temperature values" );

  registerWrapper( viewKeyStruct::defaultYieldStrengthString(), &m_defaultYieldStrength ).
    setApplyDefaultValue( DBL_MAX ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default yield strength" );

  registerWrapper( viewKeyStruct::maximumStretchString(), &m_maximumStretch ).
    setApplyDefaultValue( DBL_MAX ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum stretch" );

  registerWrapper( viewKeyStruct::maximumStretchBString(), &m_maximumStretchB ).
    setApplyDefaultValue( DBL_MAX ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum stretch" );

  registerWrapper( viewKeyStruct::maximumStretchCString(), &m_maximumStretchC ).
    setApplyDefaultValue( DBL_MAX ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum stretch" );

  registerWrapper( viewKeyStruct::maximumStretchDString(), &m_maximumStretchD ).
    setApplyDefaultValue( DBL_MAX ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum stretch" );

  // CC: TODO add defaults for plasticStrain, how to apply default for array (voigt notation)
  // Check if defaults for matrices and tensors are correctly set

  // register fields
  registerWrapper( viewKeyStruct::deformationGradientString(), &m_deformationGradient ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Array of element/particle deformation gradient values" );

  registerWrapper( viewKeyStruct::plasticStrainString(), &m_plasticStrain ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Array of element/particle plastic strain values" );

  registerWrapper( viewKeyStruct::damageString(), &m_damage ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Array of quadrature point damage values" );

  registerWrapper( viewKeyStruct::jacobianString(), &m_jacobian ).
    setApplyDefaultValue( 1.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Array of quadrature point jacobian values" );

  registerWrapper( viewKeyStruct::yieldStrengthString(), &m_yieldStrength ).
    setApplyDefaultValue( -1.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Array of element/particle yield strength values" );

  registerWrapper( viewKeyStruct::yieldStrengthBString(), &m_yieldStrengthB ).
    setApplyDefaultValue( -1.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Array of element/particle yield strength values" );

  registerWrapper( viewKeyStruct::yieldStrengthCString(), &m_yieldStrengthC ).
    setApplyDefaultValue( -1.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Array of element/particle yield strength values" );

  registerWrapper( viewKeyStruct::yieldStrengthDString(), &m_yieldStrengthD ).
    setApplyDefaultValue( -1.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Array of element/particle yield strength values" );
}


StrainHardeningPolymer::~StrainHardeningPolymer()
{}


void StrainHardeningPolymer::allocateConstitutiveData( dataRepository::Group & parent,
                                                       localIndex const numConstitutivePointsPerParentIndex )
{
  ElasticIsotropic::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_deformationGradient.resize( 0, 3, 3);
  m_plasticStrain.resize( 0, numConstitutivePointsPerParentIndex, 6);
  m_damage.resize( 0, numConstitutivePointsPerParentIndex );
  m_jacobian.resize( 0, numConstitutivePointsPerParentIndex );
  m_yieldStrength.resize( 0 );
}


void StrainHardeningPolymer::postInputInitialization()
{
  ElasticIsotropic::postInputInitialization();

  // CC: need checks for strain hardening and softening inputs
  GEOS_THROW_IF( m_strainHardeningSlope < 0.0, "Strain hardening slope must be a positive number.", InputError ); // CC: Check that these are the rules for inputs
  GEOS_THROW_IF( m_shearSofteningMagnitude < 0.0, "Shear softening magnitude must be a positive number.", InputError );
  GEOS_THROW_IF( m_shearSofteningShapeParameter1 < 0.0, "Shear softening shape paraemter 1 must be a positive number.", InputError );
  GEOS_THROW_IF( m_shearSofteningShapeParameter2 < 0.0, "Shear softening shape paraemter 2 must be a positive number.", InputError );
  GEOS_THROW_IF( m_defaultYieldStrength < 0.0, "Yield strength must be a positive number.", InputError );
  GEOS_THROW_IF( m_maximumStretch <= 1.0, "Max stretch must be greater than 1", InputError );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::yieldStrengthString() ).setApplyDefaultValue( m_defaultYieldStrength );
}


void StrainHardeningPolymer::saveConvergedState() const
{
  SolidBase::saveConvergedState();
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, StrainHardeningPolymer, std::string const &, Group * const )
}
} /* namespace geos */
