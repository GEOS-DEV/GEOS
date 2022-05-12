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
 * @file ThermalSinglePhaseFluid.cpp
 */

#include "ThermalSinglePhaseFluid.hpp"

#include "SingleFluidExtrinsicData.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

ThermalSinglePhaseFluid::ThermalSinglePhaseFluid( string const & name, Group * const parent ):
  SingleFluidBase( name, parent ), 
  m_densityPressureModelType( ExponentApproximationType::Full ), 
  m_densityTemperatureModelType( ExponentApproximationType::Full ), 
  m_viscosityModelType( ExponentApproximationType::Full ), 
  m_internalEnergyModelType( ExponentApproximationType::Linear )
{

  registerWrapper( viewKeyStruct::defaultDensityString(), &m_defaultDensity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default value for density." );

  registerWrapper( viewKeyStruct::defaultViscosityString(), &m_defaultViscosity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default value for viscosity." );

  registerWrapper( viewKeyStruct::compressibilityString(), &m_compressibility ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fluid compressibility" );

  registerWrapper( viewKeyStruct::thermalExpansionCoeffString(), &m_thermalExpansionCoeff ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fluid thermal expansion coefficient" );

  registerWrapper( viewKeyStruct::viscosibilityString(), &m_viscosibility ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fluid viscosity exponential coefficient" );

  registerWrapper( viewKeyStruct::volumetricHeatCapacityString(), &m_volumetricHeatCapacity ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Fluid volumetric heat capacity" );

  registerWrapper( viewKeyStruct::referencePressureString(), &m_referencePressure ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference pressure" );

  registerWrapper( viewKeyStruct::referenceTemperatureString(), &m_referenceTemperature ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference temperature" );

  registerWrapper( viewKeyStruct::referenceDensityString(), &m_referenceDensity ).
    setApplyDefaultValue( 1000.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference fluid density" );

  registerWrapper( viewKeyStruct::referenceViscosityString(), &m_referenceViscosity ).
    setApplyDefaultValue( 0.001 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference fluid viscosity" );

  registerWrapper( viewKeyStruct::referenceInternalEnergyString(), &m_referenceInternalEnergy ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference fluid internal energy" );

  registerWrapper( viewKeyStruct::densityPressureModelTypeString(), &m_densityPressureModelType ).
    setApplyDefaultValue( m_densityPressureModelType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Type of density model in terms of pressure . Valid options:\n* " + EnumStrings< ExponentApproximationType >::concat( "\n* " ) );

  registerWrapper( viewKeyStruct::densityTemperatureModelTypeString(), &m_densityTemperatureModelType ).
    setApplyDefaultValue( m_densityTemperatureModelType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Type of density model in terms of temperature . Valid options:\n* " + EnumStrings< ExponentApproximationType >::concat( "\n* " ) );

  registerWrapper( viewKeyStruct::viscosityModelTypeString(), &m_viscosityModelType ).
    setApplyDefaultValue( m_viscosityModelType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Type of viscosity model. Valid options:\n* " + EnumStrings< ExponentApproximationType >::concat( "\n* " ) );

  registerWrapper( viewKeyStruct::internalEnergyModelTypeString(), &m_internalEnergyModelType ).
    setApplyDefaultValue( m_internalEnergyModelType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Type of internal energy model. Valid options:\n* " + EnumStrings< ExponentApproximationType >::concat( "\n* " ) );
}

ThermalSinglePhaseFluid::~ThermalSinglePhaseFluid() = default;

void ThermalSinglePhaseFluid::allocateConstitutiveData( dataRepository::Group & parent,
                                                             localIndex const numConstitutivePointsPerParentIndex )
{
  SingleFluidBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  getExtrinsicData< extrinsicMeshData::singlefluid::density >().setApplyDefaultValue( m_defaultDensity );
  getExtrinsicData< extrinsicMeshData::singlefluid::viscosity >().setApplyDefaultValue( m_defaultViscosity );

  m_density.setValues< serialPolicy >( m_referenceDensity );
  m_viscosity.setValues< serialPolicy >( m_referenceViscosity );
  m_internalEnergy.setValues< serialPolicy >( m_referenceInternalEnergy ); 
}

void ThermalSinglePhaseFluid::postProcessInput()
{
  SingleFluidBase::postProcessInput();

  auto const checkNonnegative = [&]( real64 const value, auto const & attribute )
  {
    GEOSX_THROW_IF_LT_MSG( value, 0.0,
                           GEOSX_FMT( "{}: invalid value of attribute '{}'", getFullName(), attribute ),
                           InputError );
  };
  checkNonnegative( m_compressibility, viewKeyStruct::compressibilityString() );
  checkNonnegative( m_viscosibility, viewKeyStruct::viscosibilityString() );
  checkNonnegative( m_thermalExpansionCoeff, viewKeyStruct::thermalExpansionCoeffString() ); 
  checkNonnegative( m_volumetricHeatCapacity, viewKeyStruct::volumetricHeatCapacityString() ); 

  auto const checkPositive = [&]( real64 const value, auto const & attribute )
  {
    GEOSX_THROW_IF_LE_MSG( value, 0.0,
                           GEOSX_FMT( "{}: invalid value of attribute '{}'", getFullName(), attribute ),
                           InputError );
  };
  checkPositive( m_referenceDensity, viewKeyStruct::referenceDensityString() );
  checkPositive( m_referenceViscosity, viewKeyStruct::referenceViscosityString() );
  checkNonnegative( m_referenceInternalEnergy, viewKeyStruct::referenceInternalEnergyString() ); 

  // // Set default values for derivatives (cannot be done in base class)
  // // TODO: reconsider the necessity of this

  // real64 dRho_dP, dRho_dT;
  // real64 dVisc_dP, dVisc_dT;
  // real64 dIntE_dP, dIntE_dT; 

  // createKernelWrapper().compute( m_referencePressure, m_referenceTemperature, m_referenceDensity, dRho_dP, dRho_dT, m_referenceViscosity, dVisc_dP, dVisc_dT, m_referenceInternalEnergy, dIntE_dP, dIntE_dT );
  // getExtrinsicData< extrinsicMeshData::singlefluid::dDensity_dPressure >().setDefaultValue( dRho_dP );
  // getExtrinsicData< extrinsicMeshData::singlefluid::dDensity_dTemperature >().setDefaultValue( dRho_dT ); 
  // getExtrinsicData< extrinsicMeshData::singlefluid::dViscosity_dPressure >().setDefaultValue( dVisc_dP );
  // getExtrinsicData< extrinsicMeshData::singlefluid::dViscosity_dTemperature >().setDefaultValue( dVisc_dT ); 
  // getExtrinsicData< extrinsicMeshData::singlefluid::dInternalEnergy_dPressure >().setDefaultValue( dIntE_dP ); 
  // getExtrinsicData< extrinsicMeshData::singlefluid::dInternalEnergy_dTemperature >().setDefaultValue( dIntE_dT ); 
}

ThermalSinglePhaseFluid::KernelWrapper
ThermalSinglePhaseFluid::createKernelWrapper()
{
  return KernelWrapper( KernelWrapper::DensPresRelationType( m_referencePressure, m_referenceDensity, m_compressibility ),
                        KernelWrapper::DensTempRelationType( m_referenceTemperature, 1.0, m_thermalExpansionCoeff ),
                        KernelWrapper::ViscRelationType( m_referencePressure, m_referenceViscosity, m_viscosibility ), 
                        KernelWrapper::IntEnergyRelationType( m_referenceTemperature, m_referenceInternalEnergy, m_volumetricHeatCapacity/m_referenceInternalEnergy ), 
                        m_density,
                        m_dDensity_dPressure,
                        m_dDensity_dTemperature, 
                        m_viscosity,
                        m_dViscosity_dPressure,
                        m_dViscosity_dTemperature, 
                        m_internalEnergy, 
                        m_dInternalEnergy_dPressure,
                        m_dInternalEnergy_dTemperature );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermalSinglePhaseFluid, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geosx */
