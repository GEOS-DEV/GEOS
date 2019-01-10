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
  ConstitutiveBase( name, parent ),
  m_compressibility(),
  m_viscosibility(),
  m_referencePressure(),
  m_referenceDensity(),
  m_density(),
  m_dDensity_dPressure(),
  m_viscosity(),
  m_dViscosity_dPressure(),
  m_densityRelation( ExponentApproximationType::Linear ),
  m_viscosityRelation( ExponentApproximationType::Linear )
{
  RegisterViewWrapper( viewKeys.compressibility.Key(), &m_compressibility, 0 )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Fluid compressibility");

  RegisterViewWrapper( viewKeys.viscosibility.Key(), &m_viscosibility, 0 )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Fluid viscosity exponential coefficient");

  RegisterViewWrapper( viewKeys.referencePressure.Key(), &m_referencePressure, 0 )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Reference pressure");

  RegisterViewWrapper( viewKeys.referenceDensity.Key(), &m_referenceDensity, 0 )->
    setApplyDefaultValue(1000.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Reference fluid density");

  RegisterViewWrapper( viewKeys.referenceViscosity.Key(), &m_referenceViscosity, 0 )->
    setApplyDefaultValue(0.001)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Reference fluid viscosity");

  RegisterViewWrapper( viewKeyStruct::densityString, &m_density, 0 )->setPlotLevel(PlotLevel::LEVEL_0);
  RegisterViewWrapper( viewKeyStruct::dDens_dPresString, &m_dDensity_dPressure, 0 );

  RegisterViewWrapper( viewKeyStruct::viscosityString, &m_viscosity, 0 );
  RegisterViewWrapper( viewKeyStruct::dVisc_dPresString, &m_dViscosity_dPressure, 0 );
}

CompressibleSinglePhaseFluid::~CompressibleSinglePhaseFluid() = default;

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

void CompressibleSinglePhaseFluid::AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                                             localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );

  m_density.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dDensity_dPressure.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_density = this->m_referenceDensity;

  m_viscosity.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dViscosity_dPressure.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_viscosity = this->m_referenceViscosity;
}



void CompressibleSinglePhaseFluid::PostProcessInput()
{
  if( m_compressibility < 0.0 )
  {
    string const message = "An invalid value of fluid compressibility ("+std::to_string( m_compressibility )+") is specified";
    GEOS_ERROR( message );
  }

  if( m_viscosibility < 0.0 )
  {
    string const message = "An invalid value of fluid viscosibility ("+std::to_string( m_viscosibility )+") is specified";
    GEOS_ERROR( message );
  }

  if( m_referenceDensity <= 0.0 )
  {
    string const message = "An invalid value of reference density ("+std::to_string( m_referenceDensity )+") is specified";
    GEOS_ERROR( message );
  }

  if( m_referenceViscosity <= 0.0 )
  {
    string const message = "An invalid value of reference viscosity ("+std::to_string( m_referenceViscosity )+") is specified";
    GEOS_ERROR( message );
  }

  m_densityRelation.SetCoefficients( m_referencePressure, m_referenceDensity, m_compressibility );
  m_viscosityRelation.SetCoefficients( m_referencePressure, m_referenceViscosity, m_viscosibility );

}

void CompressibleSinglePhaseFluid::FluidDensityCompute( real64 const & pres,
                                                        localIndex const i,
                                                        real64 & dens,
                                                        real64 & dDens_dPres )
{
  m_densityRelation.Compute( pres, dens, dDens_dPres );
}


void CompressibleSinglePhaseFluid::FluidViscosityCompute( real64 const & pres,
                                                          localIndex const i,
                                                          real64 & visc,
                                                          real64 & dVisc_dPres )
{
  m_viscosityRelation.Compute( pres, visc, dVisc_dPres );
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompressibleSinglePhaseFluid, std::string const &, ManagedGroup * const )

} /* namespace constitutive */

} /* namespace geosx */
