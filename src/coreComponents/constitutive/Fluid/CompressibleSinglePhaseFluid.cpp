/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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


CompressibleSinglePhaseFluid::CompressibleSinglePhaseFluid( std::string const & name, ManagedGroup * const parent ):
  SingleFluidBase( name, parent ),
  m_densityRelation( ExponentApproximationType::Linear ),
  m_viscosityRelation( ExponentApproximationType::Linear )
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
  std::unique_ptr<CompressibleSinglePhaseFluid> newConstitutiveRelation =
    std::make_unique<CompressibleSinglePhaseFluid>( name, parent );

  newConstitutiveRelation->m_compressibility    = this->m_compressibility;
  newConstitutiveRelation->m_viscosibility      = this->m_viscosibility;
  newConstitutiveRelation->m_referencePressure  = this->m_referencePressure;
  newConstitutiveRelation->m_referenceDensity   = this->m_referenceDensity;
  newConstitutiveRelation->m_referenceViscosity = this->m_referenceViscosity;

  newConstitutiveRelation->m_densityRelation    = this->m_densityRelation;
  newConstitutiveRelation->m_viscosityRelation  = this->m_viscosityRelation;

  return std::move(newConstitutiveRelation);
}

void CompressibleSinglePhaseFluid::ProcessInputFile_PostProcess()
{
  GEOS_ERROR_IF( m_compressibility < 0.0, "An invalid value of fluid compressibility ("
                                          << m_compressibility << ") is specified" );

  GEOS_ERROR_IF( m_viscosibility < 0.0, "An invalid value of fluid viscosibility ("
                                        << m_compressibility << ") is specified" );

  GEOS_ERROR_IF( m_referenceDensity <= 0.0, "An invalid value of fluid reference density ("
                                            << m_compressibility << ") is specified" );

  GEOS_ERROR_IF( m_referenceViscosity <= 0.0, "An invalid value of fluid reference viscosity ("
                                              << m_compressibility << ") is specified" );

  m_densityRelation.SetCoefficients( m_referencePressure, m_referenceDensity, m_compressibility );
  m_viscosityRelation.SetCoefficients( m_referencePressure, m_referenceViscosity, m_viscosibility );
}

void CompressibleSinglePhaseFluid::PointUpdate( real64 const & pressure, localIndex const k, localIndex const q )
{
  Compute( pressure, m_density[k][q], m_dDensity_dPressure[k][q], m_viscosity[k][q], m_dViscosity_dPressure[k][q] );
}

void CompressibleSinglePhaseFluid::BatchUpdate( arrayView1d<double const> const & pressure )
{
#if 1
  m_densityRelation.Compute<elemPolicy>( pressure, m_density, m_dDensity_dPressure );
  m_viscosityRelation.Compute<elemPolicy>( pressure, m_viscosity, m_dViscosity_dPressure );
#else
  SingleFluidBase::BatchUpdateKernel<CompressibleSinglePhaseFluid>( pressure, m_densityRelation, m_viscosityRelation );
#endif
}

void CompressibleSinglePhaseFluid::Compute( real64 const & pressure,
                                            real64 & density, real64 & dDensity_dPressure,
                                            real64 & viscosity, real64 & dViscosity_dPressure )
{
  Compute( pressure, density, dDensity_dPressure, viscosity, dViscosity_dPressure,
           m_densityRelation, m_viscosityRelation );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompressibleSinglePhaseFluid, std::string const &, ManagedGroup * const )

} /* namespace constitutive */

} /* namespace geosx */
