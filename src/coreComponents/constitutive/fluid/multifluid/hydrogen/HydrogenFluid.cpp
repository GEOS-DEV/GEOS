/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HydrogenFluid.cpp
 */

#include "HydrogenFluid.hpp"

#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "common/format/StringUtilities.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

string HydrogenFluid::catalogName()
{
  return "HydrogenFluid";
}

HydrogenFluid::HydrogenFluid( string const & name, Group * const parent ):
  MultiFluidBase( name, parent )
{
  enableLogLevelInput();

  // if this is a thermal model, we need to make sure that the arrays will be properly displayed and saved to restart
  if( isThermal() )
  {
    getField< fields::multifluid::phaseEnthalpy >().
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setRestartFlags( RestartFlags::WRITE_AND_READ );

    getField< fields::multifluid::phaseInternalEnergy >().
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setRestartFlags( RestartFlags::WRITE_AND_READ );
  }
}

std::unique_ptr< ConstitutiveBase >
HydrogenFluid::deliverClone( string const & name, Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = MultiFluidBase::deliverClone( name, parent );

  HydrogenFluid & hydrogenFluid = dynamicCast< HydrogenFluid & >( *clone );

  hydrogenFluid.m_gasPhaseIndex = m_gasPhaseIndex;
  hydrogenFluid.m_watPhaseIndex = m_watPhaseIndex;
  hydrogenFluid.m_h2ComponentIndex = m_h2ComponentIndex;
  hydrogenFluid.m_h2oComponentIndex = m_h2oComponentIndex;

  hydrogenFluid.createModels();

  return clone;
}

integer HydrogenFluid::getWaterPhaseIndex() const
{
  return m_watPhaseIndex;
}

void HydrogenFluid::postInputInitialization()
{
  MultiFluidBase::postInputInitialization();

  integer const numPhases = numFluidPhases();
  integer const numComps = numFluidComponents();

  GEOS_THROW_IF_NE_MSG( numPhases, 2,
                        GEOS_FMT( "{}: invalid number of phases", getFullName() ),
                        InputError );
  GEOS_THROW_IF_LT_MSG( numComps, min_n_components,
                        GEOS_FMT( "{}: invalid number of components. Should be between {}  and {}.",
                                  min_n_components, max_n_components, getFullName() ),
                        InputError );
  GEOS_THROW_IF_GT_MSG( numComps, max_n_components,
                        GEOS_FMT( "{}: invalid number of components. Should be between {}  and {}.",
                                  min_n_components, max_n_components, getFullName() ),
                        InputError );

  m_gasPhaseIndex = -1;
  m_watPhaseIndex = -1;
  arrayView1d< string const > const names = phaseNames();
  for( integer ip = 0; ip < numPhases; ++ip )
  {
    string const phaseName = stringutilities::toLower( names[ip] );
    if( phaseName == "gas" )
    {
      m_gasPhaseIndex = ip;
    }
    if( phaseName == "water" )
    {
      m_watPhaseIndex = ip;
    }
  }

  GEOS_THROW_IF_LT_MSG( m_gasPhaseIndex, 0,
                        GEOS_FMT( "{}: Gas phas component not found. There should be a phase named gas.", getFullName() ),
                        InputError );

  GEOS_THROW_IF_LT_MSG( m_watPhaseIndex, 0,
                        GEOS_FMT( "{}: Water phase not found. There should be a component named water.", getFullName() ),
                        InputError );

  m_h2ComponentIndex = -1;
  m_h2oComponentIndex = -1;

  arrayView1d< string const > const compNames = componentNames();
  for( integer ic = 0; ic < numComps; ++ic )
  {
    string const compName = stringutilities::toLower( compNames[ic] );
    if( compName == "h2" )
    {
      m_h2ComponentIndex = ic;
    }
    if( compName == "h2o" )
    {
      m_h2oComponentIndex = ic;
    }
  }

  GEOS_THROW_IF_LT_MSG( m_h2ComponentIndex, 0,
                        GEOS_FMT( "{}: Hydrogen component not found. There should be a component named H2.", getFullName() ),
                        InputError );
  GEOS_THROW_IF_LT_MSG( m_h2oComponentIndex, 0,
                        GEOS_FMT( "{}: Water component not found. There should be a component named H2O.", getFullName() ),
                        InputError );

  createModels();
}

void HydrogenFluid::checkTablesParameters( real64 pressure, real64 temperature ) const
{
  GEOS_UNUSED_VAR( pressure );
  GEOS_UNUSED_VAR( temperature );
}

void HydrogenFluid::createModels()
{
  m_flash = std::make_unique< HydrogenFlash >( GEOS_FMT( "{}_HydrogenFlash", getName()),
                                               m_componentMolarWeight,
                                               m_h2ComponentIndex,
                                               m_h2oComponentIndex,
                                               m_gasPhaseIndex,
                                               m_watPhaseIndex,
                                               getLogLevel() > 0 && logger::internal::rank==0 );
}

typename HydrogenFluid::KernelWrapper HydrogenFluid::createKernelWrapper()
{
  return KernelWrapper( *m_flash,
                        m_componentMolarWeight.toViewConst(),
                        m_useMass,
                        isThermal(),
                        m_phaseFraction.toView(),
                        m_phaseDensity.toView(),
                        m_phaseMassDensity.toView(),
                        m_phaseViscosity.toView(),
                        m_phaseEnthalpy.toView(),
                        m_phaseInternalEnergy.toView(),
                        m_phaseCompFraction.toView(),
                        m_totalDensity.toView() );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, HydrogenFluid, string const &, Group * const )

} //namespace constitutive

} //namespace geos
