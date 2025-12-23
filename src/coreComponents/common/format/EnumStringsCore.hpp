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

#ifndef GEOS_COMMON_FORMAT_ENUMSTRINGSCORE_HPP
#define GEOS_COMMON_FORMAT_ENUMSTRINGSCORE_HPP

#include <type_traits>
#include <ostream>
namespace geos
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
 *
 * enum struct VTKOutputMode
 * {
 *   BINARY,
 *   ASCII
 * };
 * ENUM_STRINGS( VTKOutputMode,
 *               "binary",
 *               "ascii" );
 */
#define ENUM_STRINGS_CORE( ENUM, ... )                                     \
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
  inline std::ostream & operator<<( std::ostream & os, ENUM const e ) \
  {                                                                   \
    os << EnumStringsCore< ENUM >::toRawString( e );                  \
    return os;                                                        \
  }                                                                   \
  inline std::string toString( ENUM const e )                         \
  {                                                                   \
    return EnumStringsCore< ENUM >::toRawString( e );                 \
  }                                                                   \

template< typename ENUM >
struct EnumStringsCore
{
  using enum_type = ENUM;
  using base_type = std::underlying_type_t< ENUM >;

  static auto const & get()
  {
    return getEnumStrings( enum_type{} );
  }

  static const char * toRawString( enum_type const e )
  {
    auto const & strings = get();
    return strings[static_cast< base_type >(e)];
  }
};
}
#endif //GEOS_COMMON_FORMAT_ENUMSTRINGSCORE_HPP
