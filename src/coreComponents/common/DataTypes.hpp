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
 * @file DataTypes.hpp
 *
 * This file contains various aliases and functions that provide operations regarding the
 * use of the data types.
 */

#ifndef GEOS_COMMON_DATATYPES_HPP
#define GEOS_COMMON_DATATYPES_HPP

// Source includes
#include "common/GeosxConfig.hpp"
#include "GeosxMacros.hpp"
#include "BufferAllocator.hpp"
#include "DataLayouts.hpp"
#include "Tensor.hpp"
#include "Logger.hpp"
#include "LvArray/src/Macros.hpp"
#include "LvArray/src/Array.hpp"
#include "LvArray/src/ArrayOfArrays.hpp"
#include "LvArray/src/ArrayOfSets.hpp"
#include "LvArray/src/SparsityPattern.hpp"
#include "LvArray/src/CRSMatrix.hpp"
#include "LvArray/src/SortedArray.hpp"
#include "LvArray/src/StackBuffer.hpp"
#include "LvArray/src/ChaiBuffer.hpp"

#include "Path.hpp"

// TPL includes
#include <camp/camp.hpp>

// System includes
#ifdef GEOSX_USE_MPI
  #include <mpi.h>
#endif

#include <cassert>
//#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>
#include <typeindex>
#include <typeinfo>
#include <string>
#include <map>
#include <unordered_map>
#include <vector>
#include <set>
#include <string_view>

/**
 * top level geosx namespace contains all code that is specific to GEOSX
 */
namespace geos
{

/**
 * @brief Perform a type cast of base to derived pointer.
 * @tparam NEW_TYPE      derived pointer type
 * @tparam EXISTING_TYPE base type
 * @param val            base pointer to cast
 * @return               pointer cast to derived type or @p nullptr
 */
template< typename NEW_TYPE, typename EXISTING_TYPE >
NEW_TYPE dynamicCast( EXISTING_TYPE * const val )
{
  static_assert( std::is_pointer< NEW_TYPE >::value, "NEW_TYPE must be a pointer." );
  return dynamic_cast< NEW_TYPE >( val );
}

/**
 * @brief Perform a type cast of base to derived reference.
 * @tparam NEW_TYPE      derived reference type
 * @tparam EXISTING_TYPE base type
 * @param val            base reference to cast
 * @return               reference cast to derived type or @p nullptr
 */
template< typename NEW_TYPE, typename EXISTING_TYPE >
NEW_TYPE dynamicCast( EXISTING_TYPE & val )
{
  static_assert( std::is_reference< NEW_TYPE >::value, "NEW_TYPE must be a reference." );

  using POINTER_TO_NEW_TYPE = std::remove_reference_t< NEW_TYPE > *;
  POINTER_TO_NEW_TYPE ptr = dynamicCast< POINTER_TO_NEW_TYPE >( &val );
  GEOS_ERROR_IF( ptr == nullptr, "Cast from " << LvArray::system::demangleType( val ) << " to " <<
                 LvArray::system::demangleType< NEW_TYPE >() << " failed." );

  return *ptr;
}

/// Global MPI communicator used by GEOSX.
#ifdef GEOSX_USE_MPI
extern MPI_Comm MPI_COMM_GEOSX;
#else
extern int MPI_COMM_GEOSX;
#endif

/**
 * @name Basic data types used in GEOSX.
 */
///@{

/// Unsigned size type.
using size_t      = std::size_t;

/// Signed integer type.
using integer     = std::int32_t;

/// Local index type (for indexing objects within an MPI partition).
using localIndex  = GEOSX_LOCALINDEX_TYPE;

/// Global index type (for indexing objects across MPI partitions).
using globalIndex = GEOSX_GLOBALINDEX_TYPE;

/// String type.
using string      = std::string;

/// String type.
using string_view = std::string_view;

/// 32-bit floating point type.
using real32 = float;
/// 64-bit floating point type.
using real64 = double;

///@}

/**
 * @name Binary buffer data types.
 */
///@{

/// Type stored in communication buffers.
using buffer_unit_type = signed char;

#ifdef GEOSX_USE_CHAI
/// Type of storage for communication buffers.
using buffer_type = std::vector< buffer_unit_type, BufferAllocator< buffer_unit_type > >;
#else
/// Type of storage for communication buffers.
using buffer_type = std::vector< buffer_unit_type >;
#endif

///@}

/**
 * @name Alias for timestamp in GEOSX.
 */
///@{

/// Timestamp type (used to perform actions such a sparsity pattern computation after mesh modifications)
using Timestamp = unsigned long long int;

///@}

//START_SPHINX_INCLUDE_00

/**
 * @name Aliases for LvArray::Array class family.
 */
///@{

/// Multidimensional array type. See LvArray:Array for details.
template< typename T,
          int NDIM,
          typename PERMUTATION=camp::make_idx_seq_t< NDIM > >
using Array = LvArray::Array< T, NDIM, PERMUTATION, localIndex, LvArray::ChaiBuffer >;

/// Multidimensional array view type. See LvArray:ArrayView for details.
template< typename T,
          int NDIM,
          int USD = NDIM - 1 >
using ArrayView = LvArray::ArrayView< T, NDIM, USD, localIndex, LvArray::ChaiBuffer >;

/// Multidimensional array slice type. See LvArray:ArraySlice for details.
template< typename T, int NDIM, int USD = NDIM - 1 >
using ArraySlice = LvArray::ArraySlice< T, NDIM, USD, localIndex >;

/// Multidimensional stack-based array type. See LvArray:StackArray for details.
template< typename T, int NDIM, int MAXSIZE, typename PERMUTATION=camp::make_idx_seq_t< NDIM > >
using StackArray = LvArray::StackArray< T, NDIM, PERMUTATION, localIndex, MAXSIZE >;

///@}

/**
 * @name Short-hand aliases for commonly used array types.
 */
///@{

/// Alias for a local (stack-based) rank-1 tensor type
using R1Tensor = Tensor< real64, 3 >;
/// Alias for a local (stack-based) rank-1 tensor type using 32 bits integers
using R1Tensor32 = Tensor< real32, 3 >;

/// Alias for a local (stack-based) rank-2 Voigt tensor type
using R2SymTensor = Tensor< real64, 6 >;


/// Alias for 1D array.
template< typename T >
using array1d = Array< T, 1 >;

/// Alias for 1D array view.
template< typename T >
using arrayView1d = ArrayView< T, 1 >;

/// Alias for 1D array slice.
template< typename T, int USD = 0 >
using arraySlice1d = ArraySlice< T, 1, USD >;

/// Alias for 1D stack array.
template< typename T, int MAXSIZE >
using stackArray1d = StackArray< T, 1, MAXSIZE >;

/// Alias for 2D array.
template< typename T, typename PERMUTATION=camp::make_idx_seq_t< 2 > >
using array2d = Array< T, 2, PERMUTATION >;

/// Alias for 2D array view.
template< typename T, int USD = 1 >
using arrayView2d = ArrayView< T, 2, USD >;

/// Alias for 2D array slice.
template< typename T, int USD = 1 >
using arraySlice2d = ArraySlice< T, 2, USD >;

/// Alias for 2D stack array.
template< typename T, int MAXSIZE >
using stackArray2d = StackArray< T, 2, MAXSIZE >;

/// Alias for 3D array.
template< typename T, typename PERMUTATION=camp::make_idx_seq_t< 3 > >
using array3d = Array< T, 3, PERMUTATION >;

/// Alias for 3D array view.
template< typename T, int USD=2 >
using arrayView3d = ArrayView< T, 3, USD >;

/// Alias for 3D array slice.
template< typename T, int USD=2 >
using arraySlice3d = ArraySlice< T, 3, USD >;

/// Alias for 3D stack array.
template< typename T, int MAXSIZE >
using stackArray3d = StackArray< T, 3, MAXSIZE >;

/// Alias for 4D array.
template< typename T, typename PERMUTATION=camp::make_idx_seq_t< 4 > >
using array4d = Array< T, 4, PERMUTATION >;

/// Alias for 4D array view.
template< typename T, int USD=3 >
using arrayView4d = ArrayView< T, 4, USD >;

/// Alias for 4D array slice.
template< typename T, int USD=3 >
using arraySlice4d = ArraySlice< T, 4, USD >;

/// Alias for 4D stack array.
template< typename T, int MAXSIZE >
using stackArray4d = StackArray< T, 4, MAXSIZE >;

/// Alias for 5D array.
template< typename T, typename PERMUTATION=camp::make_idx_seq_t< 5 > >
using array5d = Array< T, 5, PERMUTATION >;

/// Alias for 5D array view.
template< typename T, int USD=4 >
using arrayView5d = ArrayView< T, 5, USD >;

/// Alias for 5D array slice.
template< typename T, int USD=4 >
using arraySlice5d = ArraySlice< T, 5, 4 >;

/// Alias for 5D stack array.
template< typename T, int MAXSIZE >
using stackArray5d = StackArray< T, 5, MAXSIZE >;

///@}

/**
 * @name Aliases for sorted arrays and set types.
 */
///@{

/// A set of local indices.
template< typename T >
using set = std::set< T >;

/// A sorted array of local indices.
template< typename T >
using SortedArray = LvArray::SortedArray< T, localIndex, LvArray::ChaiBuffer >;

/// A sorted array view of local indices.
template< typename T >
using SortedArrayView = LvArray::SortedArrayView< T, localIndex, LvArray::ChaiBuffer >;

///@}

/**
 * @name Aliases for LvArray::ArrayOfArrays class family.
 */
///@{

/// Array of variable-sized arrays. See LvArray::ArrayOfArrays for details.
template< typename T, typename INDEX_TYPE=localIndex >
using ArrayOfArrays = LvArray::ArrayOfArrays< T, INDEX_TYPE, LvArray::ChaiBuffer >;

/// View of array of variable-sized arrays. See LvArray::ArrayOfArraysView for details.
template< typename T, typename INDEX_TYPE=localIndex, bool CONST_SIZES=std::is_const< T >::value >
using ArrayOfArraysView = LvArray::ArrayOfArraysView< T, INDEX_TYPE const, CONST_SIZES, LvArray::ChaiBuffer >;

/// Array of variable-sized sets. See LvArray::ArrayOfSets for details.
template< typename T, typename INDEX_TYPE=localIndex >
using ArrayOfSets = LvArray::ArrayOfSets< T, INDEX_TYPE, LvArray::ChaiBuffer >;

/// View of array of variable-sized sets. See LvArray::ArrayOfSetsView for details.
template< typename T, typename INDEX_TYPE=localIndex >
using ArrayOfSetsView = LvArray::ArrayOfSetsView< T, INDEX_TYPE const, LvArray::ChaiBuffer >;

/// Alias for Sparsity pattern class.
template< typename COL_INDEX, typename INDEX_TYPE=localIndex >
using SparsityPattern = LvArray::SparsityPattern< COL_INDEX, INDEX_TYPE, LvArray::ChaiBuffer >;

/// Alias for Sparsity pattern View.
template< typename COL_INDEX, typename INDEX_TYPE=localIndex >
using SparsityPatternView = LvArray::SparsityPatternView< COL_INDEX, INDEX_TYPE const, LvArray::ChaiBuffer >;

/// Alias for CRS Matrix class.
template< typename T, typename COL_INDEX=globalIndex >
using CRSMatrix = LvArray::CRSMatrix< T, COL_INDEX, localIndex, LvArray::ChaiBuffer >;

/// Alias for CRS Matrix View.
template< typename T, typename COL_INDEX=globalIndex >
using CRSMatrixView = LvArray::CRSMatrixView< T, COL_INDEX, localIndex const, LvArray::ChaiBuffer >;

///@}

//END_SPHINX_INCLUDE_00

/**
 * @name Ordered and unordered map types.
 */
///@{

/**
 * @brief Base template for ordered and unordered maps.
 * @tparam TKEY key type
 * @tparam TVAL value type
 * @tparam SORTED a @p std::integral_constant<bool> indicating whether map is ordered
 */
template< typename TKEY, typename TVAL, typename SORTED >
class mapBase
{};

/// @cond DO_NOT_DOCUMENT
template< typename TKEY, typename TVAL >
class mapBase< TKEY, TVAL, std::integral_constant< bool, true > > : public std::map< TKEY, TVAL >
{
  using std::map< TKEY, TVAL >::map; // enable list initialization
};

template< typename TKEY, typename TVAL >
class mapBase< TKEY, TVAL, std::integral_constant< bool, false > > : public std::unordered_map< TKEY, TVAL >
{
  using std::unordered_map< TKEY, TVAL >::unordered_map; // enable list initialization
};
/// @endcond

/**
 * @brief Stream output operator for map types.
 * @tparam K key type
 * @tparam V value type
 * @tparam SORTED
 * @param stream output stream
 * @param map the map to print
 * @return reference to output stream
 */
template< typename K, typename V, typename SORTED >
inline
std::ostream & operator<< ( std::ostream & stream, mapBase< K, V, SORTED > const & map )
{
  stream << "{\n";
  for( auto const & pair : map )
  {
    stream << pair.first << " : " << pair.second << "\n";
  }
  stream << "}";
  return stream;
}

/// Ordered map type.
template< typename TKEY, typename TVAL >
using map = mapBase< TKEY, TVAL, std::integral_constant< bool, true > >;

/// Unordered map type.
template< typename TKEY, typename TVAL >
using unordered_map = mapBase< TKEY, TVAL, std::integral_constant< bool, false > >;

///@}

/**
 * @name Aliases for commonly used array types.
 */
///@{

/// A 1-dimensional array of geos::integer types.
using integer_array = array1d< integer >;

/// A 1-dimensional array of geos::real32 types.
using real32_array = array1d< real32 >;

/// A 1-dimensional array of geos::real64 types.
using real64_array = array1d< real64 >;

/// A 1-dimensional array of geos::string types.
using string_array = array1d< string >;

/// A 1-dimensional array of geos::Path types
using path_array = array1d< Path >;

/// A 1-dimensional array of geos::localIndex types
using localIndex_array = array1d< localIndex >;

/// A 1-dimensional array of geos::globalIndex types
using globalIndex_array = array1d< globalIndex >;


/// A 2-dimensional array of geos::integer types.
using integer_array2d = array2d< integer >;

/// A 2-dimensional array of geos::real32 types.
using real32_array2d = array2d< real32 >;

/// A 2-dimensional array of geos::real64 types.
using real64_array2d = array2d< real64 >;

/// A 2-dimensional array of geos::localIndex types
using localIndex_array2d = array2d< localIndex >;

/// A 2-dimensional array of geos::globalIndex types
using globalIndex_array2d = array2d< globalIndex >;


/// A 3-dimensional array of geos::integer types.
using integer_array3d = array3d< integer >;

/// A 3-dimensional array of geos::real32 types.
using real32_array3d = array3d< real32 >;

/// A 3-dimensional array of geos::real64 types.
using real64_array3d = array3d< real64 >;

/// A 3-dimensional array of geos::localIndex types.
using localIndex_array3d = array3d< localIndex >;

/// A 3-dimensional array of geos::globalIndex types.
using globalIndex_array3d = array3d< globalIndex >;


/// A 4-dimensional array of geos::integer types.
using integer_array4d = array4d< integer >;

/// A 4-dimensional array of geos::real32 types.
using real32_array4d = array4d< real32 >;

/// A 4-dimensional array of geos::real64 types.
using real64_array4d = array4d< real64 >;

/// A 4-dimensional array of geos::localIndex types.
using localIndex_array4d = array4d< localIndex >;

/// A 4-dimensional array of geos::globalIndex types.
using globalIndex_array4d = array4d< globalIndex >;

///@}

/// A variable for the maximum value of a geos::globalIndex.
constexpr static auto GLOBALINDEX_MAX = std::numeric_limits< globalIndex >::max();

/// A variable for the maximum value of a geos::localIndex.
constexpr static auto LOCALINDEX_MAX = std::numeric_limits< localIndex >::max();

/// A global variable for the value of a object that has not been assigned a geos::globalIndex.
constexpr static localIndex unmappedLocalIndexValue = -1;


/**
 * @brief Print a short summary of a few select type aliases.
 */
void printTypeSummary();

/**
 * @brief Class to manage the type selection of types at runtime.
 */
class rtTypes
{
public:

  /**
   * @brief Convert a @p std::type_index to a string.
   * @param key the std::type_index of the type
   * @return a hard coded string that is related to the std::type_index
   */
  static string typeNames( std::type_index const key )
  {
    static const std::unordered_map< std::type_index, string > type_names =
    {
      {std::type_index( typeid(integer)), "integer"},
      {std::type_index( typeid(real32)), "real32"},
      {std::type_index( typeid(real64)), "real64"},
      {std::type_index( typeid(localIndex)), "localIndex"},
      {std::type_index( typeid(globalIndex)), "globalIndex"},
      {std::type_index( typeid(R1Tensor)), "R1Tensor"},
      {std::type_index( typeid(R1Tensor32)), "R1Tensor32"},
      {std::type_index( typeid(R2SymTensor)), "R2SymTensor"},
      {std::type_index( typeid(integer_array)), "integer_array"},
      {std::type_index( typeid(real32_array)), "real32_array"},
      {std::type_index( typeid(real64_array)), "real64_array"},
      {std::type_index( typeid(localIndex_array)), "localIndex_array"},
      {std::type_index( typeid(globalIndex_array)), "globalIndex_array"},
      {std::type_index( typeid(integer_array2d)), "integer_array2d"},
      {std::type_index( typeid(real32_array2d)), "real32_array2d"},
      {std::type_index( typeid(real64_array2d)), "real64_array2d"},
      {std::type_index( typeid(localIndex_array2d)), "localIndex_array2d"},
      {std::type_index( typeid(globalIndex_array2d)), "globalIndex_array2d"},
      {std::type_index( typeid(integer_array3d)), "integer_array3d"},
      {std::type_index( typeid(real32_array3d)), "real32_array3d"},
      {std::type_index( typeid(real64_array3d)), "real64_array3d"},
      {std::type_index( typeid(localIndex_array3d)), "localIndex_array3d"},
      {std::type_index( typeid(globalIndex_array3d)), "globalIndex_array3d"},
      {std::type_index( typeid(real64_array4d)), "real64_array4d"},
      {std::type_index( typeid(string)), "string"},
      {std::type_index( typeid(Path)), "path"},
      {std::type_index( typeid(string_array)), "string_array"},
      {std::type_index( typeid(path_array)), "path_array"},
    };

    // If the data type is not defined here, return type_info.name()
    auto const iter = type_names.find( key );
    if( iter != type_names.end() )
    {
      return iter->second;
    }
    else
    {
      return LvArray::system::demangle( key.name());
    }
  }

  /**
   * @brief Matching regex for data types in xml.
   */
  class typeRegex
  {
public:

    /// The type of map used to store the map of type parsing regular expressions
    using regexMapType = std::map< string, string >;

    /**
     * @brief Get an iterator to the beginning of regex map.
     * @return
     */
    regexMapType::iterator begin(){return regexMap.begin();}

    /**
     * @brief Get an iterator to the end of regex map.
     * @return
     */
    regexMapType::iterator end(){return regexMap.end();}

    /**
     * @brief Get a const iterator to the beginning of regex map.
     * @return
     */
    regexMapType::const_iterator begin() const {return regexMap.begin();}

    /**
     * @brief Get a const iterator to the end of regex map.
     * @return
     */
    regexMapType::const_iterator end() const {return regexMap.end();}

private:

    /**
     * @brief Build Array regexes.
     * @param subPattern
     * @param dimension
     * @return
     *
     * @note The sub pattern is the base object you are targeting.  It can either
     *       be a simple type or a lower-dimensional array. Sub-elements and
     *       axes are given as a comma-separated list enclosed in a curly brace.
     *       For example, a 2D string array would look like: {{"a", "b"}, {"c", "d"}}
     */
    string constructArrayRegex( string subPattern, integer dimension )
    {
      if( dimension > 1 )
      {
        subPattern = constructArrayRegex( subPattern, dimension-1 );
      }

      string arrayPattern;
      if( dimension == 1 )
      {
        // Allow the bottom-level to be empty
        arrayPattern = "\\{\\s*((" + subPattern + ",\\s*)*" + subPattern + ")?\\s*\\}";
      }
      else
      {
        arrayPattern = "\\{\\s*(" + subPattern + ",\\s*)*" + subPattern + "\\s*\\}";
      }

      return arrayPattern;
    }

    // Define the component regexes:
    // Regex to match an unsigned int (123, etc.)
    // TODO c++17: Move to static constexpr std::string_view
    string ru = "[\\d]+";

    // Regex to match an signed int (-123, 455, +789, etc.)
    // TODO c++17: Move to static constexpr std::string_view
    string ri = "[+-]?[\\d]+";

    // Regex to match a float (1, +2.3, -.4, 5.6e7, 8E-9, etc.)
    // Explanation of parts:
    // [+-]?[\\d]*  matches an optional +/- at the beginning, any numbers preceding the decimal
    // ([\\d]\\.?|\\.[\\d]) matches the decimal region of the number (0, 1., 2.3, .4)
    // [\\d]*  matches any number of numbers following the decimal
    // ([eE][-+]?[\\d]+|\\s*)  matches an optional scientific notation number
    // Note: the xsd regex implementation does not allow an empty branch, so use allow whitespace at the end
    // TODO c++17: Move to static constexpr std::string_view
    string rr = "[+-]?[\\d]*([\\d]\\.?|\\.[\\d])[\\d]*([eE][-+]?[\\d]+|\\s*)";

    // Regex to match a string that can't be empty and does not contain any whitespaces nor the characters ,{}
    // TODO c++17: Move to static constexpr std::string_view
    string rs = "[^,\\{\\}\\s]+\\s*";

    // Regex to match a string that does not contain any whitespaces nor the characters ,{}
    // TODO c++17: Move to static constexpr std::string_view
    string rse = "[^,\\{\\}\\s]*\\s*";

    // Regex to match a path: a string that can't be empty and does not contain any space nor the characters *?<>|:",
    // TODO c++17: Move to static constexpr std::string_view
    string rp = "[^*?<>\\|:\";,\\s]+\\s*";

    // Regex to match a path: a string that does not contain any space nor the characters *?<>|:",
    // TODO c++17: Move to static constexpr std::string_view
    string rpe = "[^*?<>\\|:\";,\\s]*\\s*";

    // Regex to match a R1Tensor
    // TODO c++17: Move to static constexpr std::string_view
    string r1 = "\\s*\\{\\s*(" + rr + ",\\s*){2}" + rr + "\\s*\\}";

    // Regex to match a R2SymTensor
    // TODO c++17: Move to static constexpr std::string_view
    string r2s = "\\s*\\{\\s*(" + rr + ",\\s*){5}" + rr + "\\s*\\}";

    // Build master list of regexes
    regexMapType regexMap =
    {
      {"integer", ri},
      {"localIndex", ri},
      {"globalIndex", ri},
      {"real32", rr},
      {"real64", rr},
      {"R1Tensor", r1},
      {"R1Tensor32", r1},
      {"R2SymTensor", r2s},
      {"integer_array", constructArrayRegex( ri, 1 )},
      {"localIndex_array", constructArrayRegex( ri, 1 )},
      {"globalIndex_array", constructArrayRegex( ri, 1 )},
      {"real32_array", constructArrayRegex( rr, 1 )},
      {"real64_array", constructArrayRegex( rr, 1 )},
      {"integer_array2d", constructArrayRegex( ri, 2 )},
      {"localIndex_array2d", constructArrayRegex( ri, 2 )},
      {"globalIndex_array2d", constructArrayRegex( ri, 2 )},
      {"real32_array2d", constructArrayRegex( rr, 2 )},
      {"real64_array2d", constructArrayRegex( rr, 2 )},
      {"integer_array3d", constructArrayRegex( ri, 3 )},
      {"localIndex_array3d", constructArrayRegex( ri, 3 )},
      {"globalIndex_array3d", constructArrayRegex( ri, 3 )},
      {"real32_array3d", constructArrayRegex( rr, 3 )},
      {"real64_array3d", constructArrayRegex( rr, 3 )},
      {"real64_array4d", constructArrayRegex( rr, 4 )},
      {"string", rse},
      {"path", rpe},
      {"string_array", constructArrayRegex( rs, 1 )},
      {"path_array", constructArrayRegex( rp, 1 )},
      {"mapPair", rse},
      {"geos_dataRepository_PlotLevel", ri}
    };
  };
};

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
   * @brief Get the type's regex (default implementation).
   * @return empty string, indicating no custom regex
   */
  static string get() { return {}; }
};

/**
 * @brief Utility class for querying type names at runtime.
 * @tparam T the target type
 *
 * This relies on LvArray's demangling facilities and simply
 * adds some convenience methods like getting the brief name.
 */
template< typename T >
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
};

}



#endif /* GEOS_COMMON_DATATYPES_HPP */
