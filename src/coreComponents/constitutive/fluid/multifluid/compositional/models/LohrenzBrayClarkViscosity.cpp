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
 * @file LohrenzBrayClarkViscosity.cpp
 */

#include "LohrenzBrayClarkViscosity.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/compositional/CompositionalMultiphaseFluid.hpp"
#include "constitutive/fluid/multifluid/compositional/models/NullFlashModel.hpp"
#include "constitutive/fluid/multifluid/compositional/models/PhaseModel.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

LohrenzBrayClarkViscosityUpdate::LohrenzBrayClarkViscosityUpdate( MixingType const mixing_type,
  arrayView1d<real64 const> const & componentCriticalVolume )
  : m_mixing_type( mixing_type ),
  m_componentCriticalVolume(componentCriticalVolume)
{}

LohrenzBrayClarkViscosity::LohrenzBrayClarkViscosity( string const & name,
                                                      ComponentProperties const & componentProperties,
                                                      integer const phaseIndex,
                             ModelParameters const & modelParameters ):
  FunctionBase( name, componentProperties ),
  m_parameters(dynamic_cast<Parameters const*>(&modelParameters))
{
    GEOS_UNUSED_VAR(phaseIndex);
}

LohrenzBrayClarkViscosity::KernelWrapper
LohrenzBrayClarkViscosity::createKernelWrapper() const
{
  return KernelWrapper( m_mixing_type, m_parameters->m_componentCriticalVolume );
}

  void LohrenzBrayClarkViscosity::Parameters::registerParameters( MultiFluidBase * fluid )
  {
  fluid->registerWrapper( viewKeyStruct::componentCriticalVolumeString(), &m_componentCriticalVolume ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Component critical volumnes" );
  }

  void LohrenzBrayClarkViscosity::Parameters::postProcessInput( MultiFluidBase const * fluid )
  {
  integer const numComponents = fluid->numFluidComponents();
    
  if( m_componentCriticalVolume.empty() )
  {
    using Fluid = CompositionalMultiphaseFluid<NullFlashModel,NullPhaseModel,NullPhaseModel>;
    m_componentCriticalVolume.resize( numComponents );

auto const & componentCriticalPressure = fluid->getWrapper< array1d<real64> >( Fluid::viewKeyStruct::componentCriticalPressureString() ).reference();
auto const & componentCriticalTemperature = fluid->getWrapper< array1d<real64> >( Fluid::viewKeyStruct::componentCriticalTemperatureString() ).reference();

    calculateCriticalVolume( numComponents,
    componentCriticalPressure,
                             componentCriticalTemperature,
                             m_componentCriticalVolume );
  }

GEOS_THROW_IF_NE_MSG( m_componentCriticalVolume.size(), numComponents,
                          GEOS_FMT( "{}: invalid number of values in attribute '{}'", fluid->getFullName(), 
                          viewKeyStruct::componentCriticalVolumeString() ),
                          InputError );
  }

void LohrenzBrayClarkViscosity::Parameters::calculateCriticalVolume(
    integer const numComponents,
  arrayView1d< const real64 > const criticalPressure,
  arrayView1d< const real64 > const criticalTemperature,
  arrayView1d< real64 > const criticalVolume )
{
  for( integer ic=0; ic<numComponents; ++ic )
  {
    criticalVolume[ic] = 2.215e-6 * criticalTemperature[ic] / (0.025 + 1e-6*criticalPressure[ic] );   // m^3/mol
  }
}

} // end namespace compositional

} // namespace constitutive

} // end namespace geos
