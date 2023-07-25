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
  m_damage(),
  m_jacobian(),
  m_yieldStrength(),
  m_maximumStretch(),
{
  // register default values
  registerWrapper( viewKeyStruct::yieldStrengthString(), &m_yieldStrength ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Yield strength" );

  registerWrapper( viewKeyStruct::maximumStretchString(), &m_maximumStrength ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Maximum stretch" );

  registerWrapper( viewKeyStruct::strainHardeningSlopeString(), &m_strainHardeningSlope ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Strain hardening slope" );

  registerWrapper( viewKeyStruct::shearSofteningMagnitudeString(), &m_shearSofteningMagnitude ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Shear softening magnitude" );

  registerWrapper( viewKeyStruct::shearSofteningShapeParameter1String(), &m_shearSofteningShapeParameter1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Shear softening shape parameter 1" );

  registerWrapper( viewKeyStruct::shearSofteningShapeParameter2String(), &m_shearSofteningShapeParameter2 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Shear softening shape parameter 2" );

  // CC: TODO add defaults for plasticStrain, how to apply default for array (voigt notation)

  // register fields
  registerWrapper( viewKeyStruct::damageString(), &m_damage ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Array of quadrature point damage values" );

  registerWrapper( viewKeyStruct::jacobianString(), &m_jacobian ).
    setApplyDefaultValue( 1.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Array of quadrature point jacobian values" );
}


StrainHardeningPolymer::~StrainHardeningPolymer()
{}


void StrainHardeningPolymer::allocateConstitutiveData( dataRepository::Group & parent,
                                                       localIndex const numConstitutivePointsPerParentIndex )
{
  ElasticIsotropic::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_strainHardeningSlope.resize( 0, numConstitutivePointsPerParentIndex );
  m_shearSofteningMagnitude.resize( 0, numConstitutivePointsPerParentIndex );
  m_shearSofteningShapeParameter1.resize( 0, numConstitutivePointsPerParentIndex );
  m_shearSofteningShapeParameter2.resize( 0, numConstitutivePointsPerParentIndex );
  m_plasticStrain.( 0, numConstitutivePointsPerParentIndex, 6);
  m_damage.resize( 0, numConstitutivePointsPerParentIndex );
  m_jacobian.resize( 0, numConstitutivePointsPerParentIndex );
  m_yieldStrength.resize( 0, numConstitutivePointsPerParentIndex );
  m_maximumStretch.resize( 0, numConstitutivePointsPerParentIndex );
}


void StrainHardeningPolymer::postProcessInput()
{
  ElasticIsotropic::postProcessInput();

  // CC: need checks for strain hardening and softening inputs
  GEOS_THROW_IF( m_yieldStrength < 0.0, "Yield strength must be a positive number.", InputError );
  GEOS_THROW_IF( m_maxStrength <= 1.0, "Max stretch must be greater than 1", InputError );
}


void StrainHardeningPolymer::saveConvergedState() const
{
  SolidBase::saveConvergedState();
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, StrainHardeningPolymer, std::string const &, Group * const )
}
} /* namespace geos */
