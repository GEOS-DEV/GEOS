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
 * @file ThermalCompressibleSinglePhaseFluid.cpp
 */

#include "ThermalCompressibleSinglePhaseFluid.hpp"

#include "SingleFluidExtrinsicData.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

ThermalCompressibleSinglePhaseFluid::ThermalCompressibleSinglePhaseFluid( string const & name, Group * const parent ):
  CompressibleSinglePhaseFluid( name, parent ),
  m_densityPressureModelType( ExponentApproximationType::Full ),
  m_densityTemperatureModelType( ExponentApproximationType::Full ),
  m_internalEnergyModelType( ExponentApproximationType::Linear )
{

  registerWrapper( viewKeyStruct::thermalExpansionCoeffString(), &m_thermalExpansionCoeff ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fluid thermal expansion coefficient. Unit: 1/K" );

  registerWrapper( viewKeyStruct::volumetricHeatCapacityString(), &m_volumetricHeatCapacity ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fluid volumetric heat capacity. Unit: J/kg/K" );

  registerWrapper( viewKeyStruct::referenceTemperatureString(), &m_referenceTemperature ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference temperature" );

  registerWrapper( viewKeyStruct::referenceInternalEnergyString(), &m_referenceInternalEnergy ).
    setApplyDefaultValue( 0.001 ).
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

  registerWrapper( viewKeyStruct::internalEnergyModelTypeString(), &m_internalEnergyModelType ).
    setApplyDefaultValue( m_internalEnergyModelType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Type of internal energy model. Valid options:\n* " + EnumStrings< ExponentApproximationType >::concat( "\n* " ) );

}

ThermalCompressibleSinglePhaseFluid::~ThermalCompressibleSinglePhaseFluid() = default;

void ThermalCompressibleSinglePhaseFluid::allocateConstitutiveData( dataRepository::Group & parent,
                                                                    localIndex const numConstitutivePointsPerParentIndex )
{
  CompressibleSinglePhaseFluid::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_internalEnergy.setValues< serialPolicy >( m_referenceInternalEnergy );
}

void ThermalCompressibleSinglePhaseFluid::postProcessInput()
{
  CompressibleSinglePhaseFluid::postProcessInput();

  auto const checkNonnegative = [&]( real64 const value, auto const & attribute )
  {
    GEOSX_THROW_IF_LT_MSG( value, 0.0,
                           GEOSX_FMT( "{}: invalid value of attribute '{}'", getFullName(), attribute ),
                           InputError );
  };

  checkNonnegative( m_thermalExpansionCoeff, viewKeyStruct::thermalExpansionCoeffString() );
  checkNonnegative( m_volumetricHeatCapacity, viewKeyStruct::volumetricHeatCapacityString() );
  checkNonnegative( m_referenceInternalEnergy, viewKeyStruct::referenceInternalEnergyString() );

  // Due to the way update wrapper is currently implemented, we can only support one model type
  auto const checkModelType = [&]( ExponentApproximationType const value, auto const & attribute )
  {
    GEOSX_THROW_IF( value != ExponentApproximationType::Linear && value != ExponentApproximationType::Full,
                    GEOSX_FMT( "{}: invalid model type in attribute '{}' (only linear or fully exponential currently supported)", getFullName(), attribute ),
                    InputError );
  };
  checkModelType( m_densityPressureModelType, viewKeyStruct::densityPressureModelTypeString() );
  checkModelType( m_densityTemperatureModelType, viewKeyStruct::densityTemperatureModelTypeString() );
  checkModelType( m_internalEnergyModelType, viewKeyStruct::internalEnergyModelTypeString() );
}

ThermalCompressibleSinglePhaseFluid::KernelWrapper
ThermalCompressibleSinglePhaseFluid::createKernelWrapper()
{
  return KernelWrapper( KernelWrapper::DensPresRelationType( m_referencePressure, m_referenceDensity, m_compressibility ),
                        KernelWrapper::DensTempRelationType( m_referenceTemperature, 1.0, -m_thermalExpansionCoeff ),
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
                        m_dInternalEnergy_dTemperature,
                        m_enthalpy,
                        m_dEnthalpy_dPressure,
                        m_dEnthalpy_dTemperature,
                        m_referenceInternalEnergy );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermalCompressibleSinglePhaseFluid, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geosx */
