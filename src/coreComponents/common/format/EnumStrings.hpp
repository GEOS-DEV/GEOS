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
 * @file EnumStrings.hpp
 *
 * Collection of utilities to facilitate I/O of enumeration types.
 * Provides a macro definition that allows associating string names
 * with enumeration constants and a set of functions that make use
 * of these strings, like stream insertion/extraction operators.
 */

#ifndef GEOS_COMMON_FORMAT_ENUMSTRINGS_HPP
#define GEOS_COMMON_FORMAT_ENUMSTRINGS_HPP

#include "EnumStringsCore.hpp"
#include "common/format/StringUtilities.hpp"
// #include "codingUtilities/RTTypes.hpp"
#include "common/DataTypes.hpp"
#include "common/logger/Logger.hpp"
#include "common/format/Format.hpp"

#include <iostream>
#include <type_traits>
#include <algorithm>

namespace geos
{

/**
 * @brief Associate a list of string names with enumeration values.
 * @param ENUM the enumeration type
 * @param ... list of names (C-string literals)
 *
 */
#define ENUM_STRINGS( ENUM, ... )                                     \
  inline auto const & getEnumStrings( ENUM const )                    \
  {                                                                   \
    static constexpr char const * ss[] { __VA_ARGS__ };               \
    return ss;                                                        \
  }                                                                   \
                                                                      \
  inline auto const & getEnumTypeNameString( ENUM const )             \
  {                                                                   \
    return #ENUM;                                                     \
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
  }                                                                   \
                                                                      \
  inline string toString( ENUM const e )                              \
  {                                                                   \
    return EnumStrings< ENUM >::toString( e );                        \
  }                                                                   \
                                                                      \
  static_assert( std::is_enum< ENUM >::value, "Not an enumeration" )

/**
 * @brief Provides enum <-> string conversion facilities.
 * @tparam ENUM the enumeration type
 */
template< typename ENUM >
struct EnumStrings : public EnumStringsCore< ENUM >
{

  using Base = EnumStringsCore<ENUM>;
  /**
   * @brief Get a list of valid options as a delimited string.
   * @param delim delimiter (defaults to single space)
   * @return the string containing all valid strings for this type
   */
  static string concat( string const & delim = " " )
  {
    auto const & strings = Base::get();
    return stringutilities::join( std::begin( strings ), std::end( strings ), delim );
  }

  /**
   * @brief Convert enum to string.
   * @param e the enum value to convert
   * @return the corresponding string
   *
   * An error is raised if enum's numerical value is greater of equal than the number of strings.
   */
  static string toString( typename Base::enum_type const & e )
  {
    auto const & strings = Base::get();
    std::size_t size = std::distance( std::begin( strings ), std::end( strings ) );
    typename Base::base_type const index = static_cast< typename Base::base_type >( e );
    GEOS_THROW_IF( index >= LvArray::integerConversion< typename Base::base_type >( size ),
                   "Invalid value " << index << " of type " << getEnumTypeNameString( typename Base::enum_type{} ) << ". Valid range is 0.." << size - 1,
                   InputError );
    return strings[ index ];
  }

  /**
   * @brief Convert string to enum
   * @param s the string to convert
   * @return the corresponding enum value
   */
  static typename Base::enum_type fromString( string const & s )
  {
    auto const & strings = Base::get();
    auto const it = std::find( std::begin( strings ), std::end( strings ), s );
    GEOS_THROW_IF( it == std::end( strings ),
                   "Invalid value '" << s << "' of type " << getEnumTypeNameString( typename Base::enum_type{} ) << ". Valid options are: " << concat( ", " ),
                   InputError );
    typename Base::enum_type const e = static_cast< typename Base::enum_type >( LvArray::integerConversion< typename  Base::base_type >( std::distance( std::begin( strings ), it ) ) );
    return e;
  }
};

namespace internal
{
IS_VALID_EXPRESSION( HasEnumStrings, ENUM, getEnumStrings( std::declval< ENUM >() ) );
}

} // namespace geos

// Formatter specialization for enums
template< typename Enum >
struct GEOS_FMT_NS::formatter< Enum, std::enable_if_t< std::is_enum< Enum >::value && geos::internal::HasEnumStrings< Enum >, char > >
  : GEOS_FMT_NS::formatter< std::string >
{
  template< typename FormatContext >
  auto format( Enum e, FormatContext & ctx ) const
  {
    return formatter< std::string >::format( toString( e ), ctx );
  }
};

// Formatter specialization for enums
template< typename Enum >
struct GEOS_FMT_NS::formatter< Enum, std::enable_if_t< std::is_enum< Enum >::value && !geos::internal::HasEnumStrings< Enum >, char > >
  : GEOS_FMT_NS::formatter< std::underlying_type_t< Enum > >
{
  template< typename FormatContext >
  auto format( Enum e, FormatContext & ctx ) const
  {
    return GEOS_FMT_NS::formatter< std::underlying_type_t< Enum > >::format( toUnderlying( e ), ctx );
  }
};

#endif //GEOS_COMMON_FORMAT_ENUMSTRINGS_HPP
