/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file EnumStrings.hpp
 *
 * Collection of utilities to facilitate I/O of enumeration types.
 * Provides a macro definition that allows associating string names
 * with enumeration constants and a set of functions that make use
 * of these strings, like stream insertion/extraction operators.
 */

#ifndef GEOSX_COMMON_ENUM_HPP_
#define GEOSX_COMMON_ENUM_HPP_

#include "codingUtilities/StringUtilities.hpp"
#include "common/DataTypes.hpp"
#include "common/Logger.hpp"

#include <iostream>
#include <sstream>
#include <array>
#include <type_traits>
#include <algorithm>

namespace geosx
{

namespace internal
{
/**
 * @brief Simple compile-time variadic function that counts its arguments.
 * @tparam ARGS variadic pack of argument types
 * @return the number of arguments passed
 */
template< typename ... ARGS >
constexpr int countArgs( ARGS ... )
{
  return sizeof...( ARGS );
}
}

/**
 * @brief Associate a list of string names with enumeration values.
 * @param ENUM the enumeration type
 * @param ... list of names (C-string literals)
 *
 * Conditions (not enforced but won't work correctly if violated):
 *  - the macro must be called in the same namespace the enumeration type is defined in
 *  - the number and order of string arguments passed must match the enum values
 *  - enumeration constants must not have custom values assigned
 *
 * After the macro has been called, template instantiation EnumStrings<ENUM>
 * may be used to get access to strings at runtime. While not strictly necessary,
 * it is recommended that macro call immediately follows the enum definition
 * (or the class definition, if enum is defined inside a class).
 */
#define ENUM_STRINGS( ENUM, ... )                                     \
  static_assert( std::is_enum< ENUM >::value, "Not an enumeration" ); \
                                                                      \
  inline auto const & getEnumStrings( ENUM const )                    \
  {                                                                   \
    constexpr int N = internal::countArgs( __VA_ARGS__ );             \
    static constexpr std::array< char const *, N > ss{ __VA_ARGS__ }; \
    return ss;                                                        \
  }                                                                   \
                                                                      \
  inline std::ostream & operator<<( std::ostream & os, ENUM const e ) \
  {                                                                   \
    os << EnumStrings< ENUM >::toString( e );                         \
    return os;                                                        \
  }                                                                   \
                                                                      \
  inline std::istream & operator>>( std::istream & is, ENUM & e )     \
  {                                                                   \
    string s; is >> s;                                                \
    e = EnumStrings< ENUM >::fromString( s );                         \
    return is;                                                        \
  }

/**
 * @brief Provides enum <-> string conversion facilities.
 * @tparam ENUM the enumeration type
 */
template< typename ENUM >
struct EnumStrings
{
  /// Alias for the enumeration type
  using enum_type = ENUM;

  /// Alias for enum's underlying fundamental type
  using base_type = std::underlying_type_t< ENUM >;

  /**
   * @brief @return An array of strings associated with enumeration.
   */
  static auto const & get()
  {
    return getEnumStrings( enum_type{} ); // invoke ADL
  }

  /**
   * @brief Get a list of valid options as a delimited string.
   * @param delim delimiter (defaults to single space)
   * @return the string containing all valid strings for this type
   */
  static string concat( string const & delim = " " )
  {
    auto const & strings = get();
    return stringutilities::join( begin( strings ), end( strings ), delim );
  }

  /**
   * @brief Convert enum to string.
   * @param e the enum value to convert
   * @return the corresponding string
   *
   * An error is raised if enum's numerical value is greater of equal than the number of strings.
   */
  static string toString( enum_type const & e )
  {
    auto const & strings = get();
    base_type const index = static_cast< base_type >( e );
    GEOSX_ERROR_IF_GE_MSG( index, LvArray::integerConversion< base_type >( strings.size() ),
                           "Unrecognized value of " << TypeName< ENUM >::brief() << ": " << index << ".\n" <<
                           "Valid range is 0.." << strings.size() - 1 );
    return strings[ index ];
  }

  /**
   * @brief Convert string to enum
   * @param s the string to convert
   * @return the corresponding enum value
   */
  static enum_type fromString( string const & s )
  {
    auto const & strings = get();
    auto const it = std::find( begin( strings ), end( strings ), s );
    GEOSX_ERROR_IF( it == strings.end(),
                    "'" << s << "' is not a recognized value of " << TypeName< enum_type >::brief() << ".\n" <<
                    "Valid options are: \n  " << concat( "\n  " ) );
    enum_type const e = static_cast< enum_type >( LvArray::integerConversion< base_type >( std::distance( begin( strings ), it ) ) );
    return e;
  }
};

namespace internal
{
IS_VALID_EXPRESSION( HasEnumStrings, ENUM, getEnumStrings( std::declval< ENUM >() ) );
}

/**
 * @brief Specialization of TypeRegex for enumeration types with strings attached (pun intended).
 * @tparam ENUM the type of enumeration
 */
template< typename ENUM >
struct TypeRegex< ENUM, std::enable_if_t< internal::HasEnumStrings< ENUM > > >
{
  /**
   * @brief @return Regex for validating enumeration inputs for @p ENUM type.
   */
  static string get()
  {
    return EnumStrings< ENUM >::concat( "|" );
  }
};

} // namespace geosx

#endif //GEOSX_COMMON_ENUM_HPP_
