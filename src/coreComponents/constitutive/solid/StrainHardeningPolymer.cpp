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
  m_temperature(),
  m_jacobian(),
  m_yieldStrength(),
  //m_defaultBulkModulus(),
  m_bulkModulusA(),
  m_bulkModulusB(),
  m_bulkModulusT0(),
  //m_defaultShearModulus(),
  m_shearModulusA(),
  m_shearModulusB(),
  m_shearModulusT0(),
  m_defaultYieldStrength(),
  m_yieldStrengthA(),
  m_yieldStrengthB(),
  m_yieldStrengthT0(),
  m_strainHardeningSlope(),
  m_strainHardeningSlopeA(),
  m_strainHardeningSlopeB(),
  m_strainHardeningSlopeT0(),
  m_shearSofteningMagnitude(),
  m_shearSofteningMagnitudeA(),
  m_shearSofteningMagnitudeB(),
  m_shearSofteningMagnitudeT0(),
  m_shearSofteningShapeParameter1(),
  m_shearSofteningShapeParameter2(),
  m_maximumStretch(),
  m_maximumStretchA(),
  m_maximumStretchB(),
  m_maximumStretchT0()
{
  // register default values

  registerWrapper( viewKeyStruct::strainHardeningSlopeString(), &m_strainHardeningSlope ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Strain hardening slope" );

  registerWrapper( viewKeyStruct::strainHardeningSlopeAString(), &m_strainHardeningSlopeA ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Strain hardening slope A" );

  registerWrapper( viewKeyStruct::strainHardeningSlopeBString(), &m_strainHardeningSlopeB ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Strain hardening slope B" );

  registerWrapper( viewKeyStruct::strainHardeningSlopeT0String(), &m_strainHardeningSlopeT0 ).
    setApplyDefaultValue( 300.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Strain hardening slope T0" );

  registerWrapper( viewKeyStruct::shearSofteningMagnitudeString(), &m_shearSofteningMagnitude ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Shear softening magnitude" );

  registerWrapper( viewKeyStruct::shearSofteningMagnitudeAString(), &m_shearSofteningMagnitudeA ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Shear softening magnitude A" );

  registerWrapper( viewKeyStruct::shearSofteningMagnitudeBString(), &m_shearSofteningMagnitudeB ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Shear softening magnitude B" );

  registerWrapper( viewKeyStruct::shearSofteningMagnitudeT0String(), &m_shearSofteningMagnitudeT0 ).
    setApplyDefaultValue( 300 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Shear softening magnitude T0" );

  registerWrapper( viewKeyStruct::shearSofteningShapeParameter1String(), &m_shearSofteningShapeParameter1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Shear softening shape parameter 1" );

  registerWrapper( viewKeyStruct::shearSofteningShapeParameter2String(), &m_shearSofteningShapeParameter2 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Shear softening shape parameter 2" );

  registerWrapper( viewKeyStruct::defaultYieldStrengthString(), &m_defaultYieldStrength ).
    setApplyDefaultValue( DBL_MAX ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default yield strength" );

  registerWrapper( viewKeyStruct::yieldStrengthAString(), &m_yieldStrengthA ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Yield strength A" );

  registerWrapper( viewKeyStruct::yieldStrengthBString(), &m_yieldStrengthB ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Yield strength B" );

  registerWrapper( viewKeyStruct::yieldStrengthT0String(), &m_yieldStrengthT0 ).
    setApplyDefaultValue( 300 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Yield strength T0" );

  registerWrapper( viewKeyStruct::bulkModulusAString(), &m_bulkModulusA ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Temperature dependent bulk modulus A parameter" );

  registerWrapper( viewKeyStruct::bulkModulusBString(), &m_bulkModulusB ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Temperature dependent bulk modulus B parameter" );

  registerWrapper( viewKeyStruct::bulkModulusT0String(), &m_bulkModulusT0 ).
    setApplyDefaultValue( 300.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Temperature dependent bulk modulus T0 parameter" );

  registerWrapper( viewKeyStruct::shearModulusAString(), &m_shearModulusA ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Temperature dependent shear modulus A parameter" );

  registerWrapper( viewKeyStruct::shearModulusBString(), &m_shearModulusB ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Temperature dependent shear modulus B parameter" );

  registerWrapper( viewKeyStruct::shearModulusT0String(), &m_shearModulusT0 ).
    setApplyDefaultValue( 300.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Temperature dependent shear modulus T0 parameter" );

  registerWrapper( viewKeyStruct::maximumStretchString(), &m_maximumStretch ).
    setApplyDefaultValue( DBL_MAX ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum stretch" );

  registerWrapper( viewKeyStruct::maximumStretchAString(), &m_maximumStretchA ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "max Stretch A" );

  registerWrapper( viewKeyStruct::maximumStretchBString(), &m_maximumStretchB ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum stretch B" );

  registerWrapper( viewKeyStruct::maximumStretchT0String(), &m_maximumStretchT0 ).
    setApplyDefaultValue( 300.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum stretch T0" );

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

  registerWrapper( viewKeyStruct::temperatureString(), &m_temperature ).
    setApplyDefaultValue( 300.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Array of quadrature point temperature values" );

  registerWrapper( viewKeyStruct::jacobianString(), &m_jacobian ).
    setApplyDefaultValue( 1.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Array of quadrature point jacobian values" );

  registerWrapper( viewKeyStruct::yieldStrengthString(), &m_yieldStrength ).
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

  m_deformationGradient.resize( 0, 3, 3 );
  m_plasticStrain.resize( 0, numConstitutivePointsPerParentIndex, 6 );
  m_damage.resize( 0, numConstitutivePointsPerParentIndex );
  m_jacobian.resize( 0, numConstitutivePointsPerParentIndex );
  m_yieldStrength.resize( 0 );
  m_temperature.resize( 0 );
}


void StrainHardeningPolymer::postInputInitialization()
{
  ElasticIsotropic::postInputInitialization();

  // CC: need checks for strain hardening and softening inputs
  GEOS_THROW_IF( m_strainHardeningSlope < 0.0, "Strain hardening slope must be a positive number.", InputError ); // CC: Check that these
                                                                                                                  // are the rules for
                                                                                                                  // inputs
  GEOS_THROW_IF( m_shearSofteningMagnitude < 0.0, "Shear softening magnitude must be a positive number.", InputError );
  GEOS_THROW_IF( m_shearSofteningShapeParameter1 <= 0.0, "Shear softening shape paraemter 1 must be a positive number.", InputError );
  GEOS_THROW_IF( m_shearSofteningShapeParameter2 <= 0.0, "Shear softening shape paraemter 2 must be a positive number.", InputError );
  GEOS_THROW_IF( m_defaultYieldStrength < 0.0, "Yield strength must be a positive number.", InputError );
  GEOS_THROW_IF( m_maximumStretch <= 1.0, "Max stretch must be greater than 1", InputError );

  // We are using default yield strength to be the temperature-independent value, that will be modified by thermal softening.
  this->getWrapper< array1d< real64 > >( viewKeyStruct::yieldStrengthString() ).setApplyDefaultValue( m_defaultYieldStrength );
}


void StrainHardeningPolymer::saveConvergedState() const
{
  SolidBase::saveConvergedState();
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, StrainHardeningPolymer, std::string const &, Group * const )
}
} /* namespace geos */
