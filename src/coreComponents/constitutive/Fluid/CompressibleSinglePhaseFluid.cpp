/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
  * @file CompressibleSinglePhaseFluid.cpp
  */

#include "CompressibleSinglePhaseFluid.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

namespace constitutive
{

static ExponentApproximationType stringToExponentType( string const & model )
{
  if (model == "linear")
  {
    return ExponentApproximationType::Linear;
  }
  else if (model == "quadratic")
  {
    return ExponentApproximationType::Quadratic;
  }
  else if (model == "exponential")
  {
    return ExponentApproximationType::Full;
  }
  GEOS_ERROR("Model type not supported: " << model);

  // otherwise compilers complain about reaching the end of non-void function
  return ExponentApproximationType::Full;
}

CompressibleSinglePhaseFluid::CompressibleSinglePhaseFluid( std::string const & name, ManagedGroup * const parent ):
  SingleFluidBase( name, parent )
{
  RegisterViewWrapper( viewKeyStruct::compressibilityString, &m_compressibility, false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Fluid compressibility");

  RegisterViewWrapper( viewKeyStruct::viscosibilityString, &m_viscosibility, false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Fluid viscosity exponential coefficient");

  RegisterViewWrapper( viewKeyStruct::referencePressureString, &m_referencePressure, false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Reference pressure");

  RegisterViewWrapper( viewKeyStruct::referenceDensityString, &m_referenceDensity, false )->
    setApplyDefaultValue(1000.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Reference fluid density");

  RegisterViewWrapper( viewKeyStruct::referenceViscosityString, &m_referenceViscosity, false )->
    setApplyDefaultValue(0.001)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Reference fluid viscosity");

  RegisterViewWrapper( viewKeyStruct::referenceViscosityString, &m_referenceViscosity, false )->
    setApplyDefaultValue(0.001)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Reference fluid viscosity");

  RegisterViewWrapper( viewKeyStruct::densityModelString, &m_densityModelString, false )->
    setApplyDefaultValue("linear")->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Type of density model (linear, quadratic, exponential)");

  RegisterViewWrapper( viewKeyStruct::viscosityModelString, &m_viscosityModelString, false )->
    setApplyDefaultValue("linear")->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Type of viscosity model (linear, quadratic, exponential)");
}

CompressibleSinglePhaseFluid::~CompressibleSinglePhaseFluid() = default;

void CompressibleSinglePhaseFluid::AllocateConstitutiveData(dataRepository::ManagedGroup * const parent,
                                                            localIndex const numConstitutivePointsPerParentIndex)
{
  SingleFluidBase::AllocateConstitutiveData(parent, numConstitutivePointsPerParentIndex);

  m_density = m_referenceDensity;
  m_viscosity = m_referenceViscosity;
}

std::unique_ptr<ConstitutiveBase>
CompressibleSinglePhaseFluid::DeliverClone( string const & name,
                                            ManagedGroup * const parent ) const
{
  std::unique_ptr<CompressibleSinglePhaseFluid> newFluid =
    std::make_unique<CompressibleSinglePhaseFluid>( name, parent );

  newFluid->m_compressibility      = this->m_compressibility;
  newFluid->m_viscosibility        = this->m_viscosibility;
  newFluid->m_referencePressure    = this->m_referencePressure;
  newFluid->m_referenceDensity     = this->m_referenceDensity;
  newFluid->m_referenceViscosity   = this->m_referenceViscosity;
  newFluid->m_densityModelString   = this->m_densityModelString;
  newFluid->m_viscosityModelString = this->m_viscosityModelString;
  newFluid->m_densityModelType     = this->m_densityModelType;
  newFluid->m_viscosityModelType   = this->m_viscosityModelType;

  return std::move(newFluid);
}

void CompressibleSinglePhaseFluid::PostProcessInput()
{
  SingleFluidBase::PostProcessInput();

  GEOS_ERROR_IF( m_compressibility < 0.0, "An invalid value of fluid compressibility ("
                                          << m_compressibility << ") is specified" );

  GEOS_ERROR_IF( m_viscosibility < 0.0, "An invalid value of fluid viscosibility ("
                                        << m_compressibility << ") is specified" );

  GEOS_ERROR_IF( m_referenceDensity <= 0.0, "An invalid value of fluid reference density ("
                                            << m_compressibility << ") is specified" );

  GEOS_ERROR_IF( m_referenceViscosity <= 0.0, "An invalid value of fluid reference viscosity ("
                                              << m_compressibility << ") is specified" );

  m_densityModelType   = stringToExponentType( m_densityModelString );
  m_viscosityModelType = stringToExponentType( m_viscosityModelString );
}

void CompressibleSinglePhaseFluid::PointUpdate( real64 const & pressure, localIndex const k, localIndex const q )
{
  Compute( pressure, m_density[k][q], m_dDensity_dPressure[k][q], m_viscosity[k][q], m_dViscosity_dPressure[k][q] );
}

void CompressibleSinglePhaseFluid::BatchUpdate( arrayView1d<double const> const & pressure )
{
  makeExponentialRelation( m_densityModelType, m_referencePressure, m_referenceDensity, m_compressibility, [&] ( auto & relation )
  {
    SingleFluidBase::BatchDensityUpdateKernel<CompressibleSinglePhaseFluid>( pressure, relation );
  } );
  makeExponentialRelation( m_viscosityModelType, m_referencePressure, m_referenceViscosity, m_viscosibility, [&] ( auto & relation )
  {
    SingleFluidBase::BatchDensityUpdateKernel<CompressibleSinglePhaseFluid>( pressure, relation );
  } );
}

void CompressibleSinglePhaseFluid::Compute( real64 const & pressure,
                                            real64 & density, real64 & dDensity_dPressure,
                                            real64 & viscosity, real64 & dViscosity_dPressure ) const
{
  makeExponentialRelation( m_densityModelType, m_referencePressure, m_referenceDensity, m_compressibility, [&] ( auto relation )
  {
    Compute( pressure, density, dDensity_dPressure, relation );
  } );
  makeExponentialRelation( m_viscosityModelType, m_referencePressure, m_referenceViscosity, m_viscosibility, [&] ( auto relation )
  {
    Compute( pressure, viscosity, dViscosity_dPressure, relation );
  } );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompressibleSinglePhaseFluid, std::string const &, ManagedGroup * const )

} /* namespace constitutive */

} /* namespace geosx */
