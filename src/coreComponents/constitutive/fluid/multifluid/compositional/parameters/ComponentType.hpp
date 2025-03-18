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
 * @file ComponentType.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_PARAMETERS_COMPONENTTYPE_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_PARAMETERS_COMPONENTTYPE_HPP_

#include "common/DataTypes.hpp"
#include "common/format/EnumStrings.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

/**
 * @brief Class holding standard component properties for a compositional fluid model.
 */
enum class ComponentType : integer
{
  HydroCarbon,
  Water,
  Hydrogen,
  HydrogenSulphide,
  CarbonDioxide,
  Nitrogen
};

ENUM_STRINGS( ComponentType,
              "hc",
              "h2o",
              "h2",
              "h2s",
              "co2",
              "n2" );

/**
 * @brief Get the component type from the name of the component
 * @param[in] name The name of the component
 * @return returns The component type. Default value is HydroCarbon
 */
GEOS_FORCE_INLINE
static ComponentType getComponentTypeFromName( string const & name )
{
  string const s = stringutilities::toLower( name );
  auto const & knownNames = EnumStrings< ComponentType >::get();
  auto const it = std::find( std::begin( knownNames ), std::end( knownNames ), s );
  if( it == std::end( knownNames ))
  {
    return ComponentType::HydroCarbon;
  }
  return static_cast< ComponentType >( LvArray::integerConversion< integer >( std::distance( std::begin( knownNames ), it ) ) );
}

/**
 * @brief Check if a component has a specific type
 * @param[in] type The type of the component
 * @param[in] targetType The target type to check
 * @return returns true if the type matches the target type
 */
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
static bool isComponentType( integer const type, ComponentType const targetType )
{
  return (type == static_cast< int >(targetType));
}

} // namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_PARAMETERS_COMPONENTTYPE_HPP_
