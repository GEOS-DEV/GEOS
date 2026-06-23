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


#ifndef GEOS_DATAREPOSITORY_WRAPPERLIMITS_HPP_
#define GEOS_DATAREPOSITORY_WRAPPERLIMITS_HPP_

#include "codingUtilities/traits.hpp"
#include "common/DataTypes.hpp"
#include "common/format/EnumStrings.hpp"

namespace geos
{

namespace dataRepository
{

/**
 * @enum WrapperLimitsMode
 * @brief Enforcement mode associated with the limits of an attribute
 *
 * - Indicative: the limits are documentation only, no runtime check is performed.
 * - Warning:    a value outside the limits emits a runtime warning.
 * - Error:      a value outside the limits throws.
 */
enum class WrapperLimitsMode : integer
{
  Indicative,
  Warning,
  Error
};

ENUM_STRINGS( WrapperLimitsMode,
              "Indicative",
              "Warning",
              "Error" );

/**
 * @struct LimitValueType
 * @brief Structure giving the underlying value type that limits apply to.
 *
 * For a scalar T is T itself. For an array T is the type of the values in
 * the array (Array::value_type).
 */
template< typename T, bool = traits::is_array_type< T > >
struct LimitValueType
{};

template< typename T >
struct LimitValueType< T, true >
{
  using type = typename T::value_type;
};

template< typename T >
struct LimitValueType< T, false >
{
  using type = T;
};

/**
 * @brief Alias resolving to the type that limits apply to
 *
 * For a scalar value it is the scalar type itself.
 * For an array it is the type of the values in the array.
 */
template< typename T >
using limit_value_type_t = typename LimitValueType< T >::type;

/**
 * @struct is_limitable
 * @tparam T type to check
 * @brief Trait determining whether attribute limits can be applied to type @p T
 *
 * WrapperLimits apply to numeric types (integer, real32, real64, etc.) including arrays
 * of numeric types (array1d< integer >, array2d< real64 >, etc.)
 */
template< typename T >
struct is_limitable
{
  static constexpr bool value = std::is_arithmetic< limit_value_type_t< T > >::value;
};

/**
 * @brief Convenience variable template alias for is_limitable< T >::value
 */
template< typename T >
inline constexpr bool is_limitable_v = is_limitable< T >::value;

/**
 * @struct WrapperBound
 * @brief Structure containing informations about an attribute limit.
 */
template< typename T >
struct WrapperBound
{
  T value;
  bool isInclusive = true;

  /**
   * @brief WrapperBound constructor to write a limit without the "WrapperBound{ ... }" syntax
   * @param value The limit value to set
   * @param isInclusive Whether the limit should be inclusive or not
   *
   * @code
   *   .setLimits( 0.0, 1.0 )  // where setLimits takes `WrapperBound` parameters, those parameters
   *                           // can be written only with the value. The isInclusive property will
   *                           // default to true.
   * @endcode
   */
  WrapperBound( T v, bool inclusive = true )
    : value( v ), isInclusive( inclusive )
  {}

  /// @cond DO_NOT_DOCUMENT
  friend bool operator<( WrapperBound const & left, WrapperBound const & right )
  { return left.value < right.value; }

  friend bool operator>( WrapperBound const & left, WrapperBound const & right )
  { return operator<( right, left ); }

  friend bool operator<=( WrapperBound const & left, WrapperBound const & right )
  { return !( operator>( left, right ) ); }

  friend bool operator>=( WrapperBound const & left, WrapperBound const & right )
  { return !( operator<( left, right ) ); }

  friend inline std::ostream & operator<<( std::ostream & os, WrapperBound const & bound )
  {
    os << bound.value;
    return os;
  }
  ///@endcond
};

/**
 * @brief Creates an inclusive limit of @p value
 * @param value The inclusive limit value to set
 */
template< typename T >
WrapperBound< T > inclusive( T value )
{
  return WrapperBound< T >{ value, /*isInclusive*/ true };
}

/**
 * @brief Creates an exclusive limit of @p value
 * @param value The inclusive limit value to set
 */
template< typename T >
WrapperBound< T > exclusive( T value )
{
  return WrapperBound< T >{ value, /*isInclusive*/ false };
}

/**
 * @struct WrapperLimits
 * @brief Storage for the optional min/max bounds of a wrapped value.
 *
 * Specialized so that the members (std::optional< T >) are only instanciated
 * for limitable types. Preventing instantiation of non-limitable types, especially
 * abstract types that can't be instantiated with std::optional< absT >.
 */
template< typename T, bool = is_limitable_v< T > >
struct WrapperLimits
{};

template< typename T >
struct WrapperLimits< T, true >
{
  std::optional< WrapperBound< limit_value_type_t< T > > > min;
  std::optional< WrapperBound< limit_value_type_t< T > > > max;

  string getRangeStr() const
  {
    if( !this->min.has_value() && !this->max.has_value() )
    {
      return string();
    }

    string const lowerRange = this->min.has_value()
      ? GEOS_FMT( "{}{}", this->min->isInclusive ? "[" : "(", this->min->value )
      : string( "(-inf" );
    string const upperRange = this->max.has_value()
      ? GEOS_FMT( "{}{}", this->max->value, this->max->isInclusive ? "]" : ")" )
      : string( "+inf)" );
    return GEOS_FMT( "{}, {}", lowerRange, upperRange );
  }
};


/**
 * @name Helper methods
 */
///@{

/**
 * @brief Compare the given value with the min limit, taking account for the inclusive xor exclusive
 *        property of the limit.
 * @param value The value to compare to the limit
 * @param minLimit The min limit containing the inclusive
 * @return True if the value is below the min limit, false otherwise
 */
template< typename T >
bool isValueBelowMin( T const & value, WrapperBound< T > const & minLimit )
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
bool isValueAboveMax( T const & value, WrapperBound< T > const & maxLimit )
{
  return maxLimit.isInclusive ? ( value >  maxLimit.value )
                              : ( value >= maxLimit.value );
}

///@}

} /* namespace dataRepository */

} /* namespace geos */

#endif /* GEOS_DATAREPOSITORY_WRAPPERLIMITS_HPP_ */
