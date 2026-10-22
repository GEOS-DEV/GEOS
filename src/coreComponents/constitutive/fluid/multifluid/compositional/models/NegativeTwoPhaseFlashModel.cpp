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
 * @file NegativeTwoPhaseFlashModel.cpp
 */

#include "NegativeTwoPhaseFlashModel.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/CriticalVolume.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/BrineSalinity.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/FlashParameters.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/PhaseType.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

// Naming conventions
string NegativeTwoPhaseFlashModel::catalogName()
{
  return "TwoPhase";
}

NegativeTwoPhaseFlashModel::NegativeTwoPhaseFlashModel( string const & name,
                                                        ComponentProperties const & componentProperties,
                                                        ModelParameters const & modelParameters,
                                                        arrayView1d< integer const > const phaseTypes ):
  FunctionBase( name, componentProperties ),
  m_parameters( modelParameters ),
  m_phaseTypes( phaseTypes )
{}

NegativeTwoPhaseFlashModel::KernelWrapper
NegativeTwoPhaseFlashModel::createKernelWrapper() const
{
  array1d< integer > phaseOrder( KernelWrapper::getNumberOfPhases());
  calculatePhaseOrdering( m_phaseTypes, phaseOrder );
  integer const liquidIndex = phaseOrder[0];
  integer const vapourIndex = phaseOrder[1];

  EquationOfState const * equationOfState = m_parameters.get< EquationOfState >();
  EquationOfStateType const liquidEos =  EnumStrings< EquationOfStateType >::fromString( equationOfState->m_equationsOfStateNames[liquidIndex] );
  EquationOfStateType const vapourEos =  EnumStrings< EquationOfStateType >::fromString( equationOfState->m_equationsOfStateNames[vapourIndex] );

  BrineSalinity const * brineSalinity = m_parameters.get< BrineSalinity >();
  real64 const salinity = (brineSalinity == nullptr) ? 0.0 : brineSalinity->m_salinity;

  CriticalVolume const * criticalVolume = m_parameters.get< CriticalVolume >();

  FlashParameters const * flashParameters = m_parameters.get< FlashParameters >();

  return KernelWrapper( m_componentProperties.getNumberOfComponents(),
                        liquidIndex,
                        vapourIndex,
                        liquidEos,
                        vapourEos,
                        salinity,
                        criticalVolume->m_componentCriticalVolume,
                        flashParameters->m_continuousParameters,
                        flashParameters->m_discreteParameters );
}

NegativeTwoPhaseFlashModelUpdate::NegativeTwoPhaseFlashModelUpdate(
  integer const numComponents,
  integer const liquidIndex,
  integer const vapourIndex,
  EquationOfStateType const liquidEos,
  EquationOfStateType const vapourEos,
  real64 const salinity,
  arrayView1d< real64 const > const componentCriticalVolume,
  arrayView1d< real64 const > const continuousFlashParameters,
  arrayView1d< integer const > const discreteFlashParameters ):
  m_numComponents( numComponents ),
  m_liquidIndex( liquidIndex ),
  m_vapourIndex( vapourIndex ),
  m_componentCriticalVolume( componentCriticalVolume ),
  m_continuousFlashParameters( continuousFlashParameters ),
  m_discreteFlashParameters( discreteFlashParameters )
{
  m_flashData.liquidEos = liquidEos;
  m_flashData.vapourEos = vapourEos;
  m_flashData.salinity = salinity;
}

std::unique_ptr< ModelParameters >
NegativeTwoPhaseFlashModel::createParameters( std::unique_ptr< ModelParameters > parameters )
{
  std::unique_ptr< ModelParameters > params = FlashParameters::create( std::move( parameters ) );
  params = EquationOfState::create( std::move( params ) );
  params = CriticalVolume::create( std::move( params ) );
  return params;
}

void NegativeTwoPhaseFlashModel::calculatePhaseOrdering( arrayView1d< integer const > const & phaseTypes,
                                                         arrayView1d< integer > const & phaseOrder )
{
  integer li = -1;
  integer vi = -1;
  integer ai = -1;
  for( integer ip = 0; ip < KernelWrapper::getNumberOfPhases(); ++ip )
  {
    if( isPhaseType( phaseTypes[ip], PhaseType::LIQUID ))
    {
      li = ip;
    }
    if( isPhaseType( phaseTypes[ip], PhaseType::VAPOUR ))
    {
      vi = ip;
    }
    if( isPhaseType( phaseTypes[ip], PhaseType::AQUEOUS ))
    {
      ai = ip;
    }
  }

  integer liquidIndex = -1;
  integer vapourIndex = -1;
  if( vi < 0 )
  {
    // If there's no gas then there is oil and water: oil is the wetting phase
    vapourIndex = li;
    liquidIndex = ai;
  }
  else if( 0 <= li )
  {
    // Oil and gas both present: oil is the wetting phase
    vapourIndex = vi;
    liquidIndex = li;
  }
  else
  {
    // Gas and water: water is the wetting phase
    vapourIndex = vi;
    liquidIndex = ai;
  }

  phaseOrder[0] = liquidIndex;
  phaseOrder[1] = vapourIndex;
}

} // end namespace compositional

} // namespace constitutive

} // end namespace geos
