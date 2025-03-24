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
 * @file ImmiscibleWaterFlashModel.cpp
 */

#include "ImmiscibleWaterFlashModel.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/ImmiscibleWaterParameters.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/EquationOfState.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/CriticalVolume.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/BrineSalinity.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

// Naming conventions
string ImmiscibleWaterFlashModel::catalogName()
{
  return "ThreePhase";
}

ImmiscibleWaterFlashModel::ImmiscibleWaterFlashModel( string const & name,
                                                      ComponentProperties const & componentProperties,
                                                      ModelParameters const & modelParameters ):
  FunctionBase( name, componentProperties ),
  m_parameters( modelParameters )
{
  m_waterComponentIndex = ImmiscibleWaterParameters::getWaterComponentIndex( componentProperties );
}

ImmiscibleWaterFlashModel::KernelWrapper
ImmiscibleWaterFlashModel::createKernelWrapper() const
{
  constexpr integer liquidIndex = 0;
  constexpr integer vapourIndex = 1;
  constexpr integer aqueousIndex = 2;
  EquationOfState const * equationOfState = m_parameters.get< EquationOfState >();
  EquationOfStateType const liquidEos =  EnumStrings< EquationOfStateType >::fromString( equationOfState->m_equationsOfStateNames[liquidIndex] );
  EquationOfStateType const vapourEos =  EnumStrings< EquationOfStateType >::fromString( equationOfState->m_equationsOfStateNames[vapourIndex] );

  BrineSalinity const * brineSalinity = m_parameters.get< BrineSalinity >();
  real64 const salinity = (brineSalinity == nullptr) ? 0.0 : brineSalinity->m_salinity;

  CriticalVolume const * criticalVolume = m_parameters.get< CriticalVolume >();

  return KernelWrapper( m_componentProperties.getNumberOfComponents(),
                        liquidIndex,
                        vapourIndex,
                        aqueousIndex,
                        m_waterComponentIndex,
                        liquidEos,
                        vapourEos,
                        salinity,
                        criticalVolume->m_componentCriticalVolume );
}

ImmiscibleWaterFlashModelUpdate::ImmiscibleWaterFlashModelUpdate(
  integer const numComponents,
  integer const liquidIndex,
  integer const vapourIndex,
  integer const aqueousIndex,
  integer const waterComponentIndex,
  EquationOfStateType const liquidEos,
  EquationOfStateType const vapourEos,
  real64 const salinity,
  arrayView1d< real64 const > const componentCriticalVolume ):
  m_twoPhaseModel( numComponents,
                   liquidIndex,
                   vapourIndex,
                   liquidEos,
                   vapourEos,
                   salinity,
                   componentCriticalVolume ),
  m_numComponents( numComponents ),
  m_liquidIndex( liquidIndex ),
  m_vapourIndex( vapourIndex ),
  m_aquoesIndex( aqueousIndex ),
  m_waterComponentIndex( waterComponentIndex )
{}

std::unique_ptr< ModelParameters >
ImmiscibleWaterFlashModel::createParameters( std::unique_ptr< ModelParameters > parameters )
{
  auto params = NegativeTwoPhaseFlashModel::createParameters( std::move( parameters ) );
  params = ImmiscibleWaterParameters::create( std::move( params ) );
  return params;
}

} // end namespace compositional

} // namespace constitutive

} // end namespace geos
