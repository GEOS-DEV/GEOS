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
 *  @file CeramicDamage.cpp
 */

#include "CeramicDamage.hpp"
#include "SolidFields.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

CeramicDamage::CeramicDamage( string const & name, Group * const parent ):
  ElasticIsotropic( name, parent ),
  m_damage(),
  m_jacobian(),
  m_lengthScale(),
  m_strengthScale(),
  m_porosity(),
  m_referencePorosity(),
  m_tensileStrength(),
  m_compressiveStrength(),
  m_maximumStrength(),
  m_crackSpeed(),
  m_damagedMaterialFrictionSlope(),
  m_thirdInvariantDependence(),
  m_velocityGradient(),
  m_plasticStrain(),
  m_crackTipStressConcentration(),
  m_accumulatedModeIWork(),
  m_accumulatedModeIIWork(),
  m_distanceToCrackTip(),
  m_surfaceFlag()
{
  // register default values
  registerWrapper( viewKeyStruct::tensileStrengthString(), &m_tensileStrength ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Tensile strength" );

  registerWrapper( viewKeyStruct::compressiveStrengthString(), &m_compressiveStrength ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Compressive strength" );

  registerWrapper( viewKeyStruct::maximumStrengthString(), &m_maximumStrength ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Maximum theoretical strength" );

  registerWrapper( viewKeyStruct::crackSpeedString(), &m_crackSpeed ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.e16 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Crack speed" );

  // register fields
  registerWrapper( viewKeyStruct::strengthScaleString(), &m_strengthScale ).
    setApplyDefaultValue( 1.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Strength scale" );

  registerWrapper( viewKeyStruct::porosityString(), &m_porosity ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Porosity" );

  registerWrapper( viewKeyStruct::referencePorosityString(), &m_referencePorosity ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Reference porosity" );

  registerWrapper( viewKeyStruct::damageString(), &m_damage ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Array of quadrature point damage values" );

  registerWrapper( viewKeyStruct::jacobianString(), &m_jacobian ).
    setApplyDefaultValue( 1.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Jacobian state; in MPM this is synchronized to end-of-step det(F) before update" );

  registerWrapper( viewKeyStruct::lengthScaleString(), &m_lengthScale ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Array of quadrature point damage values" );

  registerWrapper( viewKeyStruct::damagedMaterialFrictionSlopeString(), &m_damagedMaterialFrictionSlope ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.5773502691896258 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Value of the damaged material friction slope" );

  registerWrapper( viewKeyStruct::thirdInvariantDependenceString(), &m_thirdInvariantDependence ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Flag to enable third invariant dependence" );

  registerWrapper( viewKeyStruct::velocityGradientString(), &m_velocityGradient ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Velocity gradient" );

  registerWrapper( viewKeyStruct::plasticStrainString(), &m_plasticStrain ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Plastic strain" );

  registerWrapper( viewKeyStruct::enableCrackTipStressConcentrationString(), &m_enableCrackTipStressConcentration ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Use crack-tip stress concentration" );

  registerWrapper( viewKeyStruct::fractureToughnessString(), &m_fractureToughness ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Fracture toughness to compute fracture radius from crack-tip stress concentration" );

  registerWrapper( viewKeyStruct::enableEnergyFailureCriterionString(), &m_enableEnergyFailureCriterion ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Enable energy failure criterion" );

  registerWrapper( viewKeyStruct::fractureEnergyReleaseRateString(), &m_fractureEnergyReleaseRate ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( DBL_MAX ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Fracture energy release rate" );

  registerWrapper( viewKeyStruct::crackTipStressConcentrationString(), &m_crackTipStressConcentration ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Crack tip stress concentration" );

  registerWrapper( viewKeyStruct::accumulatedModeIWorkString(), &m_accumulatedModeIWork ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Accumulated mode I work" );

  registerWrapper( viewKeyStruct::accumulatedModeIIWorkString(), &m_accumulatedModeIIWork ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Accumulated mode II work" );

  registerWrapper( viewKeyStruct::distanceToCrackTipString(), &m_distanceToCrackTip ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Distance to crack tip" );

  registerWrapper( viewKeyStruct::surfaceFlagString(), &m_surfaceFlag ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( 0 ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Particle surface flag" );
}

void CeramicDamage::allocateConstitutiveData( Group & parent, localIndex const numPts )
{
  ElasticIsotropic::allocateConstitutiveData( parent, numPts );

  // Either need to resize state variable arrays or restrict model to 1 pt per particle (q=1)
  GEOS_THROW_IF( numPts > 1,
               "State variables for this model are not 2D arrays that accomodate numPts > 1 per particle", InputError );

  m_strengthScale.resize( 0 );
  m_porosity.resize( 0 );
  m_referencePorosity.resize( 0 );
  m_lengthScale.resize( 0 );

  m_damage.resize( 0, numPts );
  m_jacobian.resize( 0, numPts );
  m_velocityGradient.resize( 0, 3, 3 );
  m_plasticStrain.resize( 0, numPts, 6 );

  m_crackTipStressConcentration.resize( 0 );
  m_accumulatedModeIWork.resize( 0 );
  m_accumulatedModeIIWork.resize( 0 );
  m_distanceToCrackTip.resize( 0 );
  m_surfaceFlag.resize( 0 ); 
}


void CeramicDamage::postInputInitialization()
{
  ElasticIsotropic::postInputInitialization();

  GEOS_THROW_IF( m_tensileStrength <= 0.0,
               "Tensile strength must be strictly positive.", InputError );

  GEOS_THROW_IF( m_compressiveStrength <= m_tensileStrength,
                "Compressive strength must be strictly greater than tensile strength.", InputError );

  GEOS_THROW_IF( m_maximumStrength <= m_compressiveStrength,
                "Maximum theoretical strength must be strictly greater than compressive strength.", InputError );

  GEOS_THROW_IF( m_crackSpeed <= 0.0 && m_enableEnergyFailureCriterion == 0,
                "Crack speed must be strictly positive when using time-to-failure damage.", InputError );

  GEOS_THROW_IF( m_damagedMaterialFrictionSlope <= 0.0,
                "Damaged material friction slope must be strictly positive.", InputError );

  GEOS_THROW_IF( m_thirdInvariantDependence != 0 && m_thirdInvariantDependence != 1,
                "thirdInvariantDependence must be 0 or 1.", InputError );

  GEOS_THROW_IF( m_enableCrackTipStressConcentration != 0 &&
                m_enableCrackTipStressConcentration != 1,
                "enableCrackTipStressConcentration must be 0 or 1.", InputError );

  GEOS_THROW_IF( m_enableEnergyFailureCriterion != 0 &&
                m_enableEnergyFailureCriterion != 1,
                "enableEnergyFailureCriterion must be 0 or 1.", InputError );

  GEOS_THROW_IF( m_enableCrackTipStressConcentration == 1 &&
                m_fractureToughness <= 0.0,
                "fractureToughness must be strictly positive when crack-tip stress concentration is enabled.",
                InputError );

  GEOS_THROW_IF( m_enableEnergyFailureCriterion == 1 &&
                m_fractureEnergyReleaseRate <= 0.0,
                "fractureEnergyReleaseRate must be strictly positive when the energy failure criterion is enabled.",
                InputError );

  real64 const p1 = m_compressiveStrength / 3.0;
  real64 const p2 = m_maximumStrength / m_damagedMaterialFrictionSlope;
  GEOS_THROW_IF( p2 <= p1,
               "Invalid strength surface: maximumStrength / damagedMaterialFrictionSlope "
               "must be greater than compressiveStrength / 3.",
               InputError );

}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CeramicDamage, std::string const &, Group * const )
}
} /* namespace geos */
