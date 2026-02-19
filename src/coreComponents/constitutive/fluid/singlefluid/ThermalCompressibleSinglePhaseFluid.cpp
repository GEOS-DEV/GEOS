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
 * @file ThermalCompressibleSinglePhaseFluid.cpp
 */

#include "ThermalCompressibleSinglePhaseFluid.hpp"

#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

ThermalCompressibleSinglePhaseFluid::ThermalCompressibleSinglePhaseFluid( string const & name, Group * const parent ):
  CompressibleSinglePhaseFluid( name, parent )
{
  m_numDOF=2;
  registerWrapper( viewKeyStruct::thermalExpansionCoeffString(), &m_thermalExpansionCoeff ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fluid thermal expansion coefficient. Unit: 1/K" );

  registerWrapper( viewKeyStruct::viscosityExpansivityString(), &m_viscosityExpansivity ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fluid viscosity thermal expansion coefficient at the reference temperature. Unit: 1/K" );

  registerWrapper( viewKeyStruct::specificHeatCapacityString(), &m_specificHeatCapacity ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fluid heat capacity. Unit: J/kg/K" );

  registerWrapper( viewKeyStruct::referenceTemperatureString(), &m_referenceTemperature ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference temperature" );

  registerWrapper( viewKeyStruct::referenceInternalEnergyString(), &m_referenceInternalEnergy ).
    setApplyDefaultValue( 0.001 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference fluid internal energy" );

  registerWrapper( viewKeyStruct::internalEnergyModelTypeString(), &m_internalEnergyModelType ).
    setApplyDefaultValue( ExponentApproximationType::Linear ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Type of internal energy model. Valid options:\n* " + EnumStrings< ExponentApproximationType >::concat( "\n* " ) );

}

void ThermalCompressibleSinglePhaseFluid::allocateConstitutiveData( Group & parent,
                                                                    localIndex const numPts )
{
  CompressibleSinglePhaseFluid::allocateConstitutiveData( parent, numPts );

  m_internalEnergy.value.setValues< serialPolicy >( m_referenceInternalEnergy );
}

void ThermalCompressibleSinglePhaseFluid::postInputInitialization()
{
  CompressibleSinglePhaseFluid::postInputInitialization();

  auto const checkNonnegative = [&]( real64 const value, auto const & attribute )
  {
    GEOS_THROW_IF_LT_MSG( value, 0.0,
                          GEOS_FMT( "invalid value of attribute '{}'", attribute ),
                          InputError, getDataContext() );
  };

  checkNonnegative( m_thermalExpansionCoeff, viewKeyStruct::thermalExpansionCoeffString() );
  checkNonnegative( m_viscosityExpansivity, viewKeyStruct::viscosityExpansivityString() );
  checkNonnegative( m_specificHeatCapacity, viewKeyStruct::specificHeatCapacityString() );
  checkNonnegative( m_referenceInternalEnergy, viewKeyStruct::referenceInternalEnergyString() );

  // Due to the way update wrapper is currently implemented, we can only support one model type
  auto const checkModelType = [&]( ExponentApproximationType const value, ExponentApproximationType const expectedValue, auto const & attribute )
  {
    GEOS_THROW_IF( value != expectedValue,
                   GEOS_FMT( "invalid model type in attribute '{}' (only {} currently supported)",
                             attribute, EnumStrings< ExponentApproximationType >::toString( expectedValue ) ),
                   InputError, getDataContext() );
  };
  checkModelType( m_internalEnergyModelType, ExponentApproximationType::Linear, viewKeyStruct::internalEnergyModelTypeString() );
}

ThermalCompressibleSinglePhaseFluid::KernelWrapper
ThermalCompressibleSinglePhaseFluid::createKernelWrapper()
{
  return KernelWrapper( KernelWrapper::DensRelationType( m_referencePressure, m_referenceTemperature, m_referenceDensity, m_compressibility, -m_thermalExpansionCoeff ),
                        KernelWrapper::ViscRelationType( m_referencePressure, m_referenceTemperature, m_referenceViscosity, m_viscosibility, -m_viscosityExpansivity ),
                        KernelWrapper::IntEnergyRelationType( m_referenceTemperature, m_referenceInternalEnergy, m_specificHeatCapacity/m_referenceInternalEnergy ),
                        m_density.value,
                        m_density.derivs,
                        m_viscosity.value,
                        m_viscosity.derivs,
                        m_internalEnergy.value,
                        m_internalEnergy.derivs,
                        m_enthalpy.value,
                        m_enthalpy.derivs,
                        m_referenceInternalEnergy );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermalCompressibleSinglePhaseFluid, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geos */
