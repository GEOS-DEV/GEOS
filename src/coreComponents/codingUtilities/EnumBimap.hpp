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
 * @file EnumBimap.hpp
 *
 * Collection of utilities to facilitate I/O of enumeration types.
 * Provides a macro definition that allows associating string names
 * with enumeration constants and a set of functions that make use
 * of these strings, like stream insertion/extraction operators.
 */

#ifndef GEOS_CODINGUTILITIES_ENUMBIMAP_HPP
#define GEOS_CODINGUTILITIES_ENUMBIMAP_HPP

#include "codingUtilities/StringUtilities.hpp"
#include "common/DataTypes.hpp"
#include "common/Format.hpp"

#include <iostream>
#include <type_traits>
#include <algorithm>

namespace geos
{

/**
 * @brief Associate a list of string names with enumeration values in order to provide
 * enum <-> string conversion through the EnumBimap<ENUM> static struct.
 * Compared to EnumStrings, EnumBimap supports enums that are non-continous and not starting
 * from zero, but the conversions require std::map lookups.
 * @param ENUM the enumeration type
 * @param ... list of pairs of enum values - C-string literals
 *
 * The macro must be called in the same namespace the enumeration type is defined in.
 *
 * It is recommended that macro call immediately follows the enum definition (or the class
 * definition, if enum is defined inside a class).
 *
 * @example
 * enum class LogLevel
 * { Silent = -1, Important = 0, Detailed = 1, Error = -2 };
 *
 * ENUM_BIMAP(LogLevel,
 *   { LogLevel::Silent, "Silent" },
 *   { LogLevel::Important, "Important" },
 *   { LogLevel::Detailed, "Detailed" },
 *   { LogLevel::Error, "Error" });
 */
#define ENUM_BIMAP( ENUM, ... )                                                                   \
  template< typename ENUM >                                                                       \
  inline std::map< ENUM, std::string > const & ::geos::EnumBimap< ENUM >::getToStringMap()        \
  {                                                                                               \
    static const std::map< ENUM, std::string > map( {__VA_ARGS__} );                              \
    return map;                                                                                   \
  }                                                                                               \
                                                                                                  \
  template< typename ENUM >                                                                       \
  inline constexpr size_t ::geos::EnumBimap< ENUM >::size()                                       \
  {                                                                                               \
    return std::initializer_list< std::pair< ENUM, std::string > >( {__VA_ARGS__} ).size();       \
  }                                                                                               \
                                                                                                  \
  inline std::ostream & operator<<( std::ostream & os, ENUM const e )                             \
  {                                                                                               \
    os << EnumBimap< ENUM >::toString( e );                                                       \
    return os;                                                                                    \
  }                                                                                               \
                                                                                                  \
  inline std::istream & operator>>( std::istream & is, ENUM & e )                                 \
  {                                                                                               \
    string s; is >> s;                                                                            \
    e = EnumBimap< ENUM >::fromString( s );                                                       \
    return is;                                                                                    \
  }                                                                                               \
                                                                                                  \
  static_assert( std::is_enum< ENUM >::value, "Not an enumeration" )

/**
 * @brief Provides enum <-> string conversion facilities.
 * Use the ENUM_BIMAP() macro to define the enum string labels.
 * @tparam ENUM the enumeration type
 */
template< typename ENUM >
struct EnumBimap
{
  /// Alias for enum's underlying fundamental type
  using UnderlyingType = std::underlying_type_t< ENUM >;

  /**
   * @return the reference to the the enum -> string map.
   */
  static inline std::map< ENUM, std::string > const & getToStringMap();

  /**
   * @return the reference to the the string -> enum map.
   */
  static inline std::map< std::string, ENUM > const & getToEnumMap()
  {
    auto reverseMap = []( std::map< ENUM, std::string > const & map )
    {
      std::map< std::string, ENUM > rmap;
      for( auto pair : map )
        rmap[pair.second] = pair.first;
      return rmap;
    };
    static const std::map< std::string, ENUM > m( reverseMap( getToStringMap() ) );
    return m;
  }

  /**
   * @return the reference to the strings list declared by .
   */
  static inline std::vector< std::string > const & getStrings()
  {
    auto getMapValues = []( std::map< ENUM, std::string > const & map )
    {
      std::vector< std::string > v;
      v.reserve( map.size() );
      for( auto pair : map )
        v.push_back( pair.second );
      return v;
    };
    static const std::vector< std::string > v( getMapValues( getToStringMap() ) );
    return v;
  }

  /**
   * @return the reference to the strings list.
   */
  static inline std::vector< ENUM > const & getValues()
  {
    auto getMapKeys = []( std::map< ENUM, std::string > const & map )
    {
      std::vector< ENUM > v;
      v.reserve( map.size() );
      for( auto pair : map )
        v.push_back( pair.first );
      return v;
    };
    static const std::vector< ENUM > v( getMapKeys( getToStringMap() ) );
    return v;
  }

  /**
   * @brief Convert an enum value to its string representation.
   * An error is raised if enum's numerical value has no known label.
   * @param e the enum value to convert
   * @return the corresponding string
   */
  static string toString( ENUM const & e )
  {
    try
    {
      return getToStringMap().at( e );
    } catch( std::out_of_range const & )
    {
      throw InputError( wrongValueMsg( static_cast< UnderlyingType >( e ) ) );
    }
  }

  /**
   * @brief Convert a string to its refering enum value.
   * The string can be a numeric value
   * An error is raised if the string is not reconized.
   * @param s the string to convert
   * @return the corresponding enum value
   */
  static ENUM fromString( string const & s )
  {
    // Try to parse a string value
    auto enumIt = getToEnumMap().find( s );
    if( enumIt != getToEnumMap().end() )
    {
      return enumIt->second;
    }
    else
    {
      // Try to parse a numerical value
      UnderlyingType literalValue;
      std::istringstream ss( s );
      ss >> literalValue;
      if( ss.fail() || !ss.eof() )
      {
          throw std::invalid_argument( wrongValueMsg( s ));
      }

      auto enumIt2 = getToStringMap().find( static_cast< ENUM >( literalValue ) );
      if( enumIt2 == getToStringMap().end())
      {
        throw std::invalid_argument( wrongValueMsg( s ));
      }
      return enumIt2->first;
    }
  }

  /**
   * @return a string containing a bullet point list of the valid enum values and
   * strings for this type
   */
  template< typename T >
  static string wrongValueMsg( T wrongValue )
  {
    std::ostringstream oss;
    oss << "Invalid value '" << wrongValue << "' of type " << TypeName< ENUM >::brief()
        << ". Valid values are: " << std::endl;
    for( auto pair : getToStringMap() )
      oss << " - " << static_cast< UnderlyingType >( pair.first ) << " = " << pair.second << std::endl;
    return oss.str();
  }
};


namespace internal
{
IS_VALID_EXPRESSION( HasEnumBimap, ENUM, EnumBimap< ENUM >::getEnumStrings() );
}

/**
 * @brief Specialization of TypeRegex for enumeration types with strings attached (pun intended).
 * @tparam ENUM the type of enumeration
 */
template< typename ENUM >
struct TypeRegex< ENUM, std::enable_if_t< internal::HasEnumBimap< ENUM > > >
{
  /**
   * @brief @return Regex for validating enumeration inputs for @p ENUM type.
   */
  static string get()
  {
    return EnumBimap< ENUM >::getStrings( "|" ) + "|" + EnumBimap< ENUM >::getValues( "|" );
  }
};

} // namespace geos

// Formatter specialization for enums
template< typename ENUM >
struct GEOS_FMT_NS::formatter< ENUM, std::enable_if_t< std::is_enum< ENUM >::value && geos::internal::HasEnumBimap< ENUM >, char > >
  : GEOS_FMT_NS::formatter< std::string >
{
  template< typename FormatContext >
  auto format( ENUM e, FormatContext & ctx ) const
  {
    return formatter< std::string >::format( EnumBimap< ENUM >::toString( e ), ctx );
  }
};

// Formatter specialization for enums
template< typename ENUM >
struct GEOS_FMT_NS::formatter< ENUM, std::enable_if_t< std::is_enum< ENUM >::value && !geos::internal::HasEnumBimap< ENUM >, char > >
  : GEOS_FMT_NS::formatter< std::underlying_type_t< ENUM > >
{
  template< typename FormatContext >
  auto format( ENUM e, FormatContext & ctx ) const
  {
    return GEOS_FMT_NS::formatter< std::underlying_type_t< ENUM > >::format( toUnderlying( e ), ctx );
  }
};

#endif //GEOS_CODINGUTILITIES_ENUMBIMAP_HPP
