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
 * @file CompositionalMultiphaseFluidParameters.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_COMPOSITIONALMULTIPHASEFLUIDPARAMETERS_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_COMPOSITIONALMULTIPHASEFLUIDPARAMETERS_HPP_

#include "codingUtilities/EnumStrings.hpp"

namespace geos
{

namespace constitutive
{

/// Phase types
enum class PhaseType : integer
{
  liquid,             ///< Liquid (oil)
  vapour,             ///< Vapor (gas)
  aqueous             ///< Aqueous (water)
};

/// Equation of state types
enum class EquationOfStateType : integer
{
  pr,                 ///< Peng-Robinson
  srk                 ///< Soave-Redlich-Kwong
};

/// Declare strings associated with enumeration values.
ENUM_STRINGS( PhaseType,
              "liquid",
              "vapor",
              "aqueous" );

ENUM_STRINGS( EquationOfStateType,
              "pr",
              "srk" );

} /* namespace constitutive */

} /* namespace geos */

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_COMPOSITIONALMULTIPHASEFLUIDPARAMETERS_HPP_
