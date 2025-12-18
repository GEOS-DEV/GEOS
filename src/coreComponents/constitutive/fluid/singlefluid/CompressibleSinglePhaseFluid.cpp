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
 * @file CompressibleSinglePhaseFluid.cpp
 */

#include "CompressibleSinglePhaseFluid.hpp"

#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

CompressibleSinglePhaseFluid::CompressibleSinglePhaseFluid( string const & name, Group * const parent ):
  SingleFluidBase( name, parent )
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
    setDescription( "Fluid compressibility [Pa^-1]" );

  registerWrapper( viewKeyStruct::viscosibilityString(), &m_viscosibility ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fluid viscosity exponential coefficient" );

  registerWrapper( viewKeyStruct::referencePressureString(), &m_referencePressure ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference pressure [Pa]" );

  registerWrapper( viewKeyStruct::referenceDensityString(), &m_referenceDensity ).
    setApplyDefaultValue( 1000.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference fluid density" );

  registerWrapper( viewKeyStruct::referenceViscosityString(), &m_referenceViscosity ).
    setApplyDefaultValue( 0.001 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference fluid viscosity" );

  registerWrapper( viewKeyStruct::densityModelTypeString(), &m_densityModelType ).
    setApplyDefaultValue( ExponentApproximationType::Full ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Type of density model. Valid options:\n* " + EnumStrings< ExponentApproximationType >::concat( "\n* " ) );

  registerWrapper( viewKeyStruct::viscosityModelTypeString(), &m_viscosityModelType ).
    setApplyDefaultValue( ExponentApproximationType::Linear ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Type of viscosity model. Valid options:\n* " + EnumStrings< ExponentApproximationType >::concat( "\n* " ) );

}

void CompressibleSinglePhaseFluid::allocateConstitutiveData( Group & parent, localIndex const numPts )
{
  SingleFluidBase::allocateConstitutiveData( parent, numPts );

  getField< fields::singlefluid::density >().setApplyDefaultValue( m_referenceDensity );
  getField< fields::singlefluid::viscosity >().setApplyDefaultValue( m_referenceViscosity );
}

void CompressibleSinglePhaseFluid::postInputInitialization()
{
  SingleFluidBase::postInputInitialization();

  auto const checkNonnegative = [&]( real64 const value, auto const & attribute )
  {
    GEOS_THROW_IF_LT_MSG( value, 0.0,
                          GEOS_FMT( "invalid value of attribute '{}'", attribute ),
                          InputError, getDataContext() );
  };
  checkNonnegative( m_compressibility, viewKeyStruct::compressibilityString() );
  checkNonnegative( m_viscosibility, viewKeyStruct::viscosibilityString() );

  auto const checkPositive = [&]( real64 const value, auto const & attribute )
  {
    GEOS_THROW_IF_LE_MSG( value, 0.0,
                          GEOS_FMT( "invalid value of attribute '{}'", attribute ),
                          InputError, getDataContext() );
  };
  checkPositive( m_referenceDensity, viewKeyStruct::referenceDensityString() );
  checkPositive( m_referenceViscosity, viewKeyStruct::referenceViscosityString() );

  // Due to the way update wrapper is currently implemented, we can only support one model type
  auto const checkModelType = [&]( ExponentApproximationType const value, ExponentApproximationType const expectedValue, auto const & attribute )
  {
    GEOS_THROW_IF( value != expectedValue,
                   GEOS_FMT( "invalid model type in attribute '{}' (only {} currently supported)",
                             attribute, EnumStrings< ExponentApproximationType >::toString( expectedValue ) ),
                   InputError, getDataContext() );
  };
  checkModelType( m_densityModelType, ExponentApproximationType::Full, viewKeyStruct::densityModelTypeString() );
  checkModelType( m_viscosityModelType, ExponentApproximationType::Linear, viewKeyStruct::viscosityModelTypeString() );

  // Set default values for derivatives (cannot be done in base class)
  // TODO: reconsider the necessity of this

  real64 dRho_dP;
  real64 dVisc_dP;
  createKernelWrapper().compute( m_referencePressure, m_referenceDensity, dRho_dP, m_referenceViscosity, dVisc_dP );

  for( integer i=0; i<m_density.value.size(); i++ )
  {
    m_density.derivs[0][i][DerivOffset::dP] = dRho_dP;
  }
}

CompressibleSinglePhaseFluid::KernelWrapper
CompressibleSinglePhaseFluid::createKernelWrapper()
{
  return KernelWrapper( KernelWrapper::DensRelationType( m_referencePressure, m_referenceDensity, m_compressibility ),
                        KernelWrapper::ViscRelationType( m_referencePressure, m_referenceViscosity, m_viscosibility ),
                        m_density.value,
                        m_density.derivs,
                        m_viscosity.value,
                        m_viscosity.derivs );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompressibleSinglePhaseFluid, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geos */
