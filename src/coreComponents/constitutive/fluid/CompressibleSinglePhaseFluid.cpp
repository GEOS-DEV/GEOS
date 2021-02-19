/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CompressibleSinglePhaseFluid.cpp
 */

#include "CompressibleSinglePhaseFluid.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

CompressibleSinglePhaseFluid::CompressibleSinglePhaseFluid( string const & name, Group * const parent ):
  SingleFluidBase( name, parent ),
  m_densityModelType( ExponentApproximationType::Linear ),
  m_viscosityModelType( ExponentApproximationType::Linear )
{
  registerWrapper( viewKeyStruct::compressibilityString(), &m_compressibility ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fluid compressibility" );

  registerWrapper( viewKeyStruct::viscosibilityString(), &m_viscosibility ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fluid viscosity exponential coefficient" );

  registerWrapper( viewKeyStruct::referencePressureString(), &m_referencePressure ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference pressure" );

  registerWrapper( viewKeyStruct::referenceDensityString(), &m_referenceDensity ).
    setApplyDefaultValue( 1000.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference fluid density" );

  registerWrapper( viewKeyStruct::referenceViscosityString(), &m_referenceViscosity ).
    setApplyDefaultValue( 0.001 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference fluid viscosity" );

  registerWrapper( viewKeyStruct::densityModelTypeString(), &m_densityModelType ).
    setApplyDefaultValue( m_densityModelType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Type of density model. Valid options:\n* " + EnumStrings< ExponentApproximationType >::concat( "\n* " ) );

  registerWrapper( viewKeyStruct::viscosityModelTypeString(), &m_viscosityModelType ).
    setApplyDefaultValue( m_viscosityModelType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Type of viscosity model. Valid options:\n* " + EnumStrings< ExponentApproximationType >::concat( "\n* " ) );
}

CompressibleSinglePhaseFluid::~CompressibleSinglePhaseFluid() = default;

void CompressibleSinglePhaseFluid::allocateConstitutiveData( dataRepository::Group & parent,
                                                             localIndex const numConstitutivePointsPerParentIndex )
{
  SingleFluidBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_density.setValues< serialPolicy >( m_referenceDensity );
  m_viscosity.setValues< serialPolicy >( m_referenceViscosity );
}

void CompressibleSinglePhaseFluid::postProcessInput()
{
  SingleFluidBase::postProcessInput();

  GEOSX_ERROR_IF_LT_MSG( m_compressibility, 0.0,
                         getName() << ": invalid value of " << viewKeyStruct::compressibilityString() );

  GEOSX_ERROR_IF_LT_MSG( m_viscosibility, 0.0,
                         getName() << ": invalid value of " << viewKeyStruct::viscosibilityString() );

  GEOSX_ERROR_IF_LE_MSG( m_referenceDensity, 0.0,
                         getName() << ": invalid value of " << viewKeyStruct::referenceDensityString() );

  GEOSX_ERROR_IF_LE_MSG( m_referenceViscosity, 0.0,
                         getName() << ": invalid value of " << viewKeyStruct::referenceViscosityString() );

  // Due to the way update wrapper is currently implemented, we can only support one model type

  GEOSX_ERROR_IF( m_densityModelType != ExponentApproximationType::Linear,
                  getName() << ": model type currently not supported: " << m_densityModelType );

  GEOSX_ERROR_IF( m_viscosityModelType != ExponentApproximationType::Linear,
                  getName() << ": model type currently not supported: " << m_viscosityModelType );

  // Set default values for derivatives (cannot be done in base class)
  // TODO: reconsider the necessity of this

  real64 dRho_dP;
  real64 dVisc_dP;
  createKernelWrapper().compute( m_referencePressure, m_referenceDensity, dRho_dP, m_referenceViscosity, dVisc_dP );
  this->getWrapper< array2d< real64 > >( viewKeyStruct::dDens_dPresString() ).setDefaultValue( dRho_dP );
  this->getWrapper< array2d< real64 > >( viewKeyStruct::dVisc_dPresString() ).setDefaultValue( dVisc_dP );
}

CompressibleSinglePhaseFluid::KernelWrapper
CompressibleSinglePhaseFluid::createKernelWrapper()
{
  return KernelWrapper( KernelWrapper::DensRelationType( m_referencePressure, m_referenceDensity, m_compressibility ),
                        KernelWrapper::ViscRelationType( m_referencePressure, m_referenceViscosity, m_viscosibility ),
                        m_density,
                        m_dDensity_dPressure,
                        m_viscosity,
                        m_dViscosity_dPressure );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompressibleSinglePhaseFluid, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geosx */
