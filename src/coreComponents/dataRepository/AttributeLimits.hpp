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


#ifndef GEOS_DATAREPOSITORY_ATTRIBUTELIMITS_HPP_
#define GEOS_DATAREPOSITORY_ATTRIBUTELIMITS_HPP_

#include "common/DataTypes.hpp"
#include "common/format/EnumStrings.hpp"
#include <optional>
#include <type_traits>

namespace geos
{

namespace dataRepository
{

/**
 * @enum LimitsMode
 * @brief Enforcement mode associated with the limits of an attribute
 *
 * - Indicative: the limits are documentation only, no runtime check is performed.
 * - Warning:    a value outside the limits emits a runtime warning.
 * - Error:      a value outside the limits throws.
 */
enum class LimitsMode : integer
{
  Indicative,
  Warning,
  Error
};

ENUM_STRINGS( LimitsMode,
              "Indicative",
              "Warning",
              "Error" );

/**
 * @struct is_limitable
 * @tparam T type to check
 * @brief Trait determining whether attribute limits can be applied to type @p T
 *
 * Limits apply to scalar numeric types (integer, real32, real64, etc.)
 */
template< typename T >
struct is_limitable
{
  static constexpr bool value = std::is_arithmetic< T >::value;
};

/**
 * @brief Convenience variable template alias for is_limitable< T >::value
 */
template< typename T >
inline constexpr bool is_limitable_v = is_limitable< T >::value;

/**
 * @struct Bound
 * @brief Structure containing informations about an attribute limit.
 */
template< typename T >
struct Bound
{
  T value;
  bool isInclusive = true;

  /**
   * @brief Bound constructor to write a limit without the "Bound{ ... }" syntax
   * @param value The limit value to set
   * @param isInclusive Wether the limit should be inclusive or not
   *
   * @code
   *   .setLimits( 0.0, 1.0 )  // where setLimits takes `Bound` parameters, those parameters can
   *                           // be written only with the value. The isInclusive property will
   *                           // default to true.
   * @endcode
   */
  Bound( T v, bool inclusive = true )
    : value( v ), isInclusive( inclusive )
  {}
};

/**
 * @brief Creates an inclusive limit of @p value
 * @param value The inclusive limit value to set
 */
template< typename T >
Bound< T > inclusive( T value )
{
  return Bound< T >{ value, /*isInclusive*/ true };
}

/**
 * @brief Creates an exclusive limit of @p value
 * @param value The inclusive limit value to set
 */
template< typename T >
Bound< T > exclusive( T value )
{
  return Bound< T >{ value, /*isInclusive*/ false };
}

/**
 * @struct Limits
 * @brief Storage for the optional min/max bounds of a wrapped value.
 *
 * Specialized so that the members (std::optional< T >) are only instanciated
 * for limitable types. Preventing instantiation of non-limitable types, especially
 * abstract types that can't be instantiated with std::optional< absT >.
 */
template< typename T, bool = is_limitable_v< T > >
struct Limits
{};

template< typename T >
struct Limits< T, true >
{
  std::optional< Bound< T > > min;
  std::optional< Bound< T > > max;
};


// Helper methods

/**
 * @brief Compare the given value with the min limit, taking account for the inclusive xor exclusive
 *        property of the limit.
 * @param value The value to compare to the limit
 * @param minLimit The min limit containing the inclusive
 * @return True if the value is below the min limit, false otherwise
 */
template< typename T >
static bool isValueBelowMin( T const & value, Bound< T > const & minLimit )
{
  return minLimit.isInclusive ? ( value <  minLimit.value )
                              : ( value <= minLimit.value );
}

/**
 * @brief Compare the given value with the max limit, taking account for the inclusive xor exclusive
 *        property of the limit.
 * @param value The value to compare to the limit
 * @param maxLimit The max limit containing the inclusive
 * @return True if the value is above the max limit, false otherwise
 */
template< typename T >
static bool isValueAboveMax( T const & value, Bound< T > const & maxLimit )
{
  return maxLimit.isInclusive ? ( value >  maxLimit.value )
                              : ( value >= maxLimit.value );
}

} /* namespace dataRepository */

} /* namespace geos */

#endif /* GEOS_DATAREPOSITORY_ATTRIBUTELIMITS_HPP_ */
