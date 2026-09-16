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
 * @file DamageCloggingPermeability.cpp
 */

#include "DamageCloggingPermeability.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{


DamageCloggingPermeability::DamageCloggingPermeability( string const & name, Group * const parent ):
  DamagePermeability( name, parent )
{
  registerWrapper( viewKeyStruct::cloggingExponentString(), &m_cloggingExponent ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 3.0 ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Exponent n of the clogging multiplier (1 - cloggedPoreFraction)^n" );

  registerWrapper( viewKeyStruct::minCloggingMultiplierString(), &m_minCloggingMultiplier ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0e-6 ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Lower bound of the clogging multiplier, so a fully clogged cell keeps a residual permeability" );

  registerWrapper( viewKeyStruct::cloggedPoreFractionString(), &m_cloggedPoreFraction ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Fraction of the pre-precipitation pore space filled by precipitated minerals" );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, DamageCloggingPermeability, string const &, Group * const )

}
} /* namespace geos */
