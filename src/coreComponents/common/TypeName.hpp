/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TypeName.hpp
 *
 * Collection of utilities to query type names in different forms (brief/full).
 */

#ifndef GEOSX_COMMON_TYPENAME_HPP_
#define GEOSX_COMMON_TYPENAME_HPP_

#include "common/DataTypes.hpp"

namespace geosx
{

namespace internal
{
template< typename T >
struct TypeAliasName : std::false_type {};

/**
 * @brief Generate a specialization of TypeAliasName that maps an aliased type to its alias name.
 * @param TYPE the type
 */
#define GEOSX_TYPE_ALIAS_NAME( TYPE )         \
template<>                                    \
struct TypeAliasName< TYPE > : std::true_type \
{                                             \
  static constexpr auto name = #TYPE;         \
}

GEOSX_TYPE_ALIAS_NAME( integer );
GEOSX_TYPE_ALIAS_NAME( localIndex );
GEOSX_TYPE_ALIAS_NAME( globalIndex );
GEOSX_TYPE_ALIAS_NAME( real32 );
GEOSX_TYPE_ALIAS_NAME( real64 );
GEOSX_TYPE_ALIAS_NAME( string );
GEOSX_TYPE_ALIAS_NAME( R1Tensor );
GEOSX_TYPE_ALIAS_NAME( Path );

#undef GEOSX_TYPE_ALIAS_NAME

}

/**
 * @brief Utility class for querying type names at runtime.
 * @tparam T the target type
 * @tparam ENABLE used to conditionally enable specializations
 *
 * This relies on LvArray's demangling facilities and simply
 * adds some convenience methods like getting the brief name.
 */
template< typename T, typename ENABLE = void >
struct TypeName
{
  /**
   * @brief @return Full name of the type.
   */
  static string full()
  {
    return ::LvArray::system::demangle( typeid( T ).name() );
  }

  /**
   * @brief @return brief name of the type (ignoring namespaces).
   */
  static string brief()
  {
    string const full_name = full();
    string::size_type const pos = full_name.find_last_of( "::" );
    return ( pos == string::npos ) ? full_name : full_name.substr( pos );
  }

  /**
   * @brief @return an XML-safe name of the type (no invalid characters)
   */
  static string safe()
  {
    return std::regex_replace( full(), std::regex( R"r(::|<|>|,|\s+)r" ), "_" );
  }
};

template< typename T >
struct TypeName< T, std::enable_if_t< internal::TypeAliasName< T >::value > >
{
  /**
   * @brief @return Full name of the type.
   */
  static string full()
  {
    return internal::TypeAliasName< T >::name;
  }

  /**
   * @brief @return brief name of the type (ignoring namespaces).
   */
  static string brief()
  {
    return full();
  }

  /**
   * @brief @return an XML-safe name of the type (no invalid characters)
   */
  static string safe()
  {
    return full();
  }
};

/**
 * @brief Specialization of TypeName for LvArray::Array.
 * @tparam T array type
 * @tparam NDIM array number of dimensions
 * @tparam PERMUTATION array permutation
 * @tparam INDEX_TYPE array index type
 * @tparam BUFFER_TYPE array buffer type
 */
template< typename T,
          int NDIM,
          typename PERMUTATION,
          typename INDEX_TYPE,
          template< typename > class BUFFER_TYPE >
struct TypeName< LvArray::Array< T, NDIM, PERMUTATION, INDEX_TYPE, BUFFER_TYPE > >
{
  /**
   * @brief @return Full name of the type.
   */
  static string full()
  {
    return ::LvArray::system::demangle( typeid( LvArray::Array< T, NDIM, PERMUTATION, INDEX_TYPE, BUFFER_TYPE > ).name() );
  }

  /**
   * @brief @return brief name of the type (ignoring namespaces).
   */
  static string brief()
  {
    return string( "Array<" ) + TypeName< T >::brief() + ", " + std::to_string( NDIM ) + ">";
  }

  /**
   * @brief @return an XML-safe name of the type (no invalid characters)
   */
  static string safe()
  {
    return TypeName< T >::safe() + "_array" + ( NDIM > 1 ? std::to_string( NDIM ) + "d" : "" );
  }
};

/**
 * @brief Builds a validation regex for a multi-dim array
 * @param subPattern pattern for the array element type
 * @param dimension number of array dimensions
 * @return the validation regex string
 */
inline std::string constructArrayRegex( std::string subPattern, integer const dimension )
{
  if( dimension > 1 )
  {
    subPattern = constructArrayRegex( subPattern, dimension-1 );
  }

  std::string arrayPattern;
  if( dimension == 1 )
  {
    // Allow the bottom-level to be empty
    arrayPattern = R"r(\{\s*((()r" + subPattern + R"r(),\s*)*()r" + subPattern + R"r())?\s*\})r";
  }
  else
  {
    arrayPattern = R"r(\{\s*(()r" + subPattern + R"r(),\s*)*()r" + subPattern + R"r()\s*\})r";
  }

  return arrayPattern;
}

/**
 * @brief Extension point for custom types to provide a validation regexp to schema.
 * @tparam T the type for which the regex is defined
 * @tparam ENABLE used to conditionally enable partial specializations
 *
 * Specializations should define the following method:
 * \code{cpp}
 *   static string get();
 * \endcode
 */
template< typename T, typename ENABLE = void >
struct TypeRegex
{
  /**
   * @brief @return Validation regex the type.
   */
  static string get() { return ".*"; }
};

/**
 * @brief Specialization of TypeRegex for unsigned integer types.
 * @tparam T target type
 */
template< typename T >
struct TypeRegex< T, std::enable_if_t< std::is_integral< T >::value && std::is_unsigned< T >::value > >
{
  /**
   * @brief @return Validation regex the type.
   */
  static string get()
  {
    return R"r([\d]+)r";
  }
};

/**
 * @brief Specialization of TypeRegex for signed integer types.
 * @tparam T target type
 */
template< typename T >
struct TypeRegex< T, std::enable_if_t< std::is_integral< T >::value && std::is_signed< T >::value > >
{
  /**
   * @brief @return Validation regex the type.
   */
  static string get()
  {
    return R"r([+-]?[\d]+)r";
  }
};

/**
 * @brief Specialization of TypeRegex for floating point types.
 * @tparam T target type
 */
template< typename T >
struct TypeRegex< T, std::enable_if_t< std::is_floating_point< T >::value > >
{
  /**
   * @brief @return Validation regex the type.
   */
  static string get()
  {
    // Explanation of parts:
    //  * [+-]?[\\d]*  matches an optional +/- at the beginning, any numbers preceding the decimal
    //  * ([\\d]\\.?|\\.[\\d]) matches the decimal region of the number (0, 1., 2.3, .4)
    //  * [\\d]*  matches any number of numbers following the decimal
    //  * ([eE][-+]?[\\d]+|\\s*)  matches an optional scientific notation number
    // Note: the xsd regex implementation does not allow an empty branch, so use allow whitespace at the end
    return R"r([+-]?[\d]*([\d]\.?|\.[\d])[\d]*([eE][-+]?[\d]+|\s*))r";
  }
};

/**
 * @brief Specialization of TypeRegex for string.
 * @tparam T target type
 */
template< typename T >
struct TypeRegex< T, std::enable_if_t< std::is_base_of< string, T >::value > >
{
  /**
   * @brief @return Validation regex the type.
   */
  static string get()
  {
    return R"r([^,\{\}]*)r";
  }
};

/**
 * @brief Specialization of TypeRegex for R1Tensor.
 * @tparam T target type
 */
template<>
struct TypeRegex< R1Tensor >
{
  /**
   * @brief @return Validation regex the type.
   */
  static string get()
  {
    string const rr = TypeRegex< double >::get();
    return string( R"r(\s*()r" ) + rr + R"r(,\s*){2})r" + rr;
  }
};

/**
 * @brief Specialization of TypeRegex for arrays.
 * @tparam T type of array element
 * @tparam NDIM number of array dimensions
 * @tparam PERMUTATION array permutation
 * @tparam INDEX_TYPE array index type
 * @tparam BUFFER_TYPE array buffer type
 */
template< typename T,
          int NDIM,
          typename PERMUTATION,
          typename INDEX_TYPE,
          template< typename > class BUFFER_TYPE >
struct TypeRegex< LvArray::Array< T, NDIM, PERMUTATION, INDEX_TYPE, BUFFER_TYPE > >
{
  /**
   * @brief @return Validation regex the type, built out of regexes for subtype @p T.
   */
  static string get()
  {
    return constructArrayRegex( TypeRegex< T >::get(), NDIM );
  }
};

} // namespace geosx

#endif //GEOSX_COMMON_TYPENAME_HPP_
