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
 * @file PhaseType.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_PARAMETERS_PHASETYPE_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_PARAMETERS_PHASETYPE_HPP_

#include "common/DataTypes.hpp"
#include "common/format/EnumStrings.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

/**
 * @brief A list of fluid phase types
 */
enum class PhaseType : integer
{
  LIQUID = 0,
  VAPOUR = 1,
  AQUEOUS = 2
};

ENUM_STRINGS( PhaseType,
              "Liquid",
              "Vapour",
              "Aqueous" );

/**
 * @brief Get the phase type from the name of the phase
 * @param[in] name The name of the phase
 * @return returns The phase type. Throws if the phase is not known
 */
PhaseType getPhaseTypeFromName( string const & name );

/**
 * @brief Check if a phase is a specific type
 * @param[in] type The type of the phase
 * @param[in] targetType The target type to check
 * @return returns true if the type matches the target type
 */
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
static bool isPhaseType( integer const type, PhaseType const targetType )
{
  return (type == static_cast< int >(targetType));
}

} // namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_PARAMETERS_PHASETYPE_HPP_
