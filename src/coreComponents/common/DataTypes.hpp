/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
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

#ifndef GEOSX_COMMON_DATATYPES_HPP
#define GEOSX_COMMON_DATATYPES_HPP

// Source includes
#include "common/GeosxConfig.hpp"
#include "common/GeosxMacros.hpp"
#include "common/BufferAllocator.hpp"
#include "common/DataLayouts.hpp"
#include "Logger.hpp"
#include "cxx-utilities/src/Macros.hpp"
#include "cxx-utilities/src/Array.hpp"
#include "cxx-utilities/src/ArrayOfArrays.hpp"
#include "cxx-utilities/src/ArrayOfSets.hpp"
#include "cxx-utilities/src/CRSMatrix.hpp"
#include "cxx-utilities/src/Macros.hpp"
#include "cxx-utilities/src/SortedArray.hpp"
#include "cxx-utilities/src/StackBuffer.hpp"

#include "math/TensorT/TensorT.h"
#include "Path.hpp"

// TPL includes
#include <camp/camp.hpp>

// System includes
#ifdef GEOSX_USE_MPI
  #include <mpi.h>
#endif

#include <cassert>
#include <cmath>
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

/// macro definition to specify whether or not to use dynamic_cast
#ifndef USE_DYNAMIC_CASTING
  #define USE_DYNAMIC_CASTING 1
#endif

/**
 * top level geosx namespace contains all code that is specific to GEOSX
 */
namespace geosx
{

/**
 * @brief Perform a type cast of base to derived pointer.
 * @tparam NEW_TYPE      derived pointer type
 * @tparam EXISTING_TYPE base type
 * @param val            base pointer to cast
 * @return               pointer cast to derived type or @p nullptr
 *
 * Depending on value of @p USE_DYNAMIC_CASTING, will use either
 * @p dynamic_cast or @p static_cast. The latter could result in undefined
 * behavior if the cast is invalid (e.g. EXISTING_TYPE not base of @p NEW_TYPE)
 */
template< typename NEW_TYPE, typename EXISTING_TYPE >
NEW_TYPE dynamicCast( EXISTING_TYPE * const val )
{
  static_assert( std::is_pointer< NEW_TYPE >::value, "NEW_TYPE must be a pointer." );

#if USE_DYNAMIC_CASTING
  return dynamic_cast< NEW_TYPE >( val );
#else
  return static_cast< NEW_TYPE >( val );
#endif
}

/**
 * @brief Perform a type cast of base to derived reference.
 * @tparam NEW_TYPE      derived reference type
 * @tparam EXISTING_TYPE base type
 * @param val            base reference to cast
 * @return               reference cast to derived type or @p nullptr
 *
 * Depending on value of @p USE_DYNAMIC_CASTING, will use either
 * @p dynamic_cast or @p static_cast. The latter could result in undefined
 * behavior if the cast is invalid (e.g. EXISTING_TYPE not base of @p NEW_TYPE)
 */
template< typename NEW_TYPE, typename EXISTING_TYPE >
NEW_TYPE dynamicCast( EXISTING_TYPE & val )
{
  static_assert( std::is_reference< NEW_TYPE >::value, "NEW_TYPE must be a reference." );

  using POINTER_TO_NEW_TYPE = std::remove_reference_t< NEW_TYPE > *;
  POINTER_TO_NEW_TYPE ptr = dynamicCast< POINTER_TO_NEW_TYPE >( &val );
  GEOSX_ERROR_IF( ptr == nullptr, "Cast failed." );

  return *ptr;
}

/// Global MPI communicator used by GEOSX.
#ifdef GEOSX_USE_MPI
extern MPI_Comm MPI_COMM_GEOSX;
#else
constexpr int MPI_COMM_GEOSX = 0;
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
using localIndex  = std::ptrdiff_t;
/// Global index type (for indexing objects across MPI partitions).
using globalIndex = long long int;//std::int64_t;
/// String type.
using string      = std::string;

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

#ifdef USE_CHAI
/// Type of storage for communication buffers.
using buffer_type = std::vector< buffer_unit_type, buffer_allocator< buffer_unit_type > >;
#else
/// Type of storage for communication buffers.
using buffer_type = std::vector< buffer_unit_type >;
#endif

///@}

/**
 * @name Aliases for LvArray::Array class family.
 */
///@{

/// Multidimensional array type. See LvArray:Array for details.
template< typename T,
          int NDIM,
          typename PERMUTATION=camp::make_idx_seq_t< NDIM >,
          template< typename > class DATA_VECTOR_TYPE=LvArray::NewChaiBuffer >
using Array = LvArray::Array< T, NDIM, PERMUTATION, localIndex, DATA_VECTOR_TYPE >;

/// Multidimensional array view type. See LvArray:ArrayView for details.
template< typename T,
          int NDIM,
          int USD = NDIM - 1,
          template< typename > class DATA_VECTOR_TYPE=LvArray::NewChaiBuffer >
using ArrayView = LvArray::ArrayView< T, NDIM, USD, localIndex, DATA_VECTOR_TYPE >;

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
template< typename T >
using array4d = Array< T, 4 >;

/// Alias for 4D array view.
template< typename T >
using arrayView4d = ArrayView< T, 4 >;

/// Alias for 4D array slice.
template< typename T >
using arraySlice4d = ArraySlice< T, 4 >;

/// Alias for 4D stack array.
template< typename T, int MAXSIZE >
using stackArray4d = StackArray< T, 4, MAXSIZE >;

/// Alias for 5D array.
template< typename T >
using array5d = Array< T, 5 >;

/// Alias for 5D array view.
template< typename T >
using arrayView5d = ArrayView< T, 5 >;

/// Alias for 5D array slice.
template< typename T >
using arraySlice5d = ArraySlice< T, 5 >;

/// Alias for 5D stack array.
template< typename T, int MAXSIZE >
using stackArray5d = StackArray< T, 5, MAXSIZE >;


/// Alias for CRS Matrix class.
template< typename T, typename ROWINDEX, typename COLINDEX >
using CRSMatrix = LvArray::CRSMatrix< T, COLINDEX, ROWINDEX >;

/// Alias for CRS Matrix View.
template< typename T, typename COLINDEX, typename LINEEARINDEX >
using CRSMatrixView = LvArray::CRSMatrixView< T, COLINDEX, LINEEARINDEX >;

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
using SortedArray = LvArray::SortedArray< T, localIndex >;

/// A sorted array view of local indices.
template< typename T >
using SortedArrayView = LvArray::SortedArrayView< T, localIndex >;

///@}

/**
 * @name Aliases for LvArray::ArrayOfArrays class family.
 */
///@{

/// Array of variable-sized arrays. See LvArray::ArrayOfArrays for details.
template< typename T >
using ArrayOfArrays = LvArray::ArrayOfArrays< T, localIndex >;

/// View of array of variable-sized arrays. See LvArray::ArrayOfArraysView for details.
template< typename T, bool CONST_SIZES=std::is_const< T >::value >
using ArrayOfArraysView = LvArray::ArrayOfArraysView< T, localIndex const, CONST_SIZES >;

/// Array of variable-sized sets. See LvArray::ArrayOfSets for details.
template< typename T >
using ArrayOfSets = LvArray::ArrayOfSets< T, localIndex >;

/// View of array of variable-sized sets. See LvArray::ArrayOfSetsView for details.
template< typename T >
using ArrayOfSetsView = LvArray::ArrayOfSetsView< T, localIndex const >;

///@}


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
{};

template< typename TKEY, typename TVAL >
class mapBase< TKEY, TVAL, std::integral_constant< bool, false > > : public std::unordered_map< TKEY, TVAL >
{};

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
/// @endcond

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

using integer_array        = array1d< integer >;
using integer_const_array  = array1d< integer const >;

using real32_array        = array1d< real32 >;
using real32_const_array  = array1d< real32 const >;

using real64_array        = array1d< real64 >;
using real64_const_array  = array1d< real64 const >;

using string_array        = array1d< string >;
using string_const_array  = array1d< string const >;

using path_array        = array1d< Path >;
using path_const_array  = array1d< Path const >;

using localIndex_array        = array1d< localIndex >;
using localIndex_const_array  = array1d< localIndex const >;

using globalIndex_array        = array1d< globalIndex >;
using globalIndex_const_array  = array1d< globalIndex const >;

using mpiBuffer = array1d< char >;

using integer_set        = set< integer >;
using integer_const_set  = set< integer const >;

using real32_set        = set< real32 >;
using real32_const_set  = set< real32 const >;

using real64_set        = set< real64 >;
using real64_const_set  = set< real64 const >;

using string_set        = set< string >;
using string_const_set  = set< string const >;

using localIndex_set        = set< localIndex >;
using localIndex_const_set  = set< localIndex const >;

using globalIndex_set        = set< globalIndex >;
using globalIndex_const_set  = set< globalIndex const >;



using integer_array2d       = array2d< integer >;
using integer_const_array2d = array2d< integer const >;

using real32_array2d       = array2d< real32 >;
using real32_const_array2d = array2d< real32 const >;

using real64_array2d       = array2d< real64 >;
using real64_const_array2d = array2d< real64 const >;

using string_array2d       = array2d< string >;
using string_const_array2d = array2d< string const >;

using localIndex_array2d       = array2d< localIndex >;
using localIndex_const_array2d = array2d< localIndex const >;

using globalIndex_array2d       = array2d< globalIndex >;
using globalIndex_const_array2d = array2d< globalIndex const >;



using integer_array3d       = array3d< integer >;
using integer_const_array3d = array3d< integer const >;

using real32_array3d       = array3d< real32 >;
using real32_const_array3d = array3d< real32 const >;

using real64_array3d       = array3d< real64 >;
using real64_const_array3d = array3d< real64 const >;

using string_array3d       = array3d< string >;
using string_const_array3d = array3d< string const >;

using localIndex_array3d       = array3d< localIndex >;
using localIndex_const_array3d = array3d< localIndex const >;

using globalIndex_array3d       = array3d< globalIndex >;
using globalIndex_const_array3d = array3d< globalIndex const >;

///@}

/**
 * @name Legacy typedefs.
 */
///@{

using r1_array = array1d< R1Tensor >;
using r2_array = array1d< R2Tensor >;
using r2Sym_array = array1d< R2SymTensor >;

using r1_array2d= array2d< R1Tensor >;
using r2_array2d= array2d< R2Tensor >;
using r2Sym_array2d= array2d< R2SymTensor >;

//using mapPair = std::pair<integer, localIndex>;
using mapPair_array = std::pair< localIndex_array, localIndex_array >;


constexpr static auto GLOBALINDEX_MAX = std::numeric_limits< globalIndex >::max();
constexpr static auto LOCALINDEX_MAX = std::numeric_limits< localIndex >::max();

constexpr static localIndex unmappedLocalIndexValue = -1;

///@}


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
  static std::string typeNames( std::type_index const key )
  {
    const std::unordered_map< std::type_index, std::string > type_names =
    {
      {std::type_index( typeid(integer)), "integer"},
      {std::type_index( typeid(real32)), "real32"},
      {std::type_index( typeid(real64)), "real64"},
      {std::type_index( typeid(localIndex)), "localIndex"},
      {std::type_index( typeid(globalIndex)), "globalIndex"},
      {std::type_index( typeid(R1Tensor)), "R1Tensor"},
      {std::type_index( typeid(R2Tensor)), "R2Tensor"},
      {std::type_index( typeid(R2SymTensor)), "R2SymTensor"},
      {std::type_index( typeid(integer_array)), "integer_array"},
      {std::type_index( typeid(real32_array)), "real32_array"},
      {std::type_index( typeid(real64_array)), "real64_array"},
      {std::type_index( typeid(localIndex_array)), "localIndex_array"},
      {std::type_index( typeid(globalIndex_array)), "globalIndex_array"},
      {std::type_index( typeid(r1_array)), "r1_array"},
      {std::type_index( typeid(r2_array)), "r2_array"},
      {std::type_index( typeid(r2Sym_array)), "r2Sym_array"},
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
      {std::type_index( typeid(r1_array2d)), "r1_array2d"},
      {std::type_index( typeid(r2_array2d)), "r2_array2d"},
      {std::type_index( typeid(r2Sym_array2d)), "r2Sym_array2d"},
      {std::type_index( typeid(string)), "string"},
      {std::type_index( typeid(Path)), "path"},
      {std::type_index( typeid(string_array)), "string_array"},
      {std::type_index( typeid(path_array)), "path_array"},
      {std::type_index( typeid(mapPair_array)), "mapPair_array"}
    };

    // If the data type is not defined here, return type_info.name()
    auto tmp = type_names.find( key );
    if( tmp != type_names.end())
    {
      return type_names.at( key );
    }
    else
    {
      return cxx_utilities::demangle( key.name());
    }
  }


  /**
   * @brief A set of enums for each geosx defined data type.
   */
  enum class TypeIDs
  {
    integer_id,          //!< integer_id
    localIndex_id,       //!< localIndex_id
    globalIndex_id,      //!< globalIndex_id
    real32_id,           //!< real32_id
    real64_id,           //!< real64_id
    r1Tensor_id,         //!< r1Tensor_id
    r2Tensor_id,         //!< r2Tensor_id
    r2SymTensor_id,      //!< r2SymTensor_id
    integer_array_id,    //!< integer_array_id
    localIndex_array_id, //!< localIndex_array_id
    globalIndex_array_id,//!< globalIndex_array_id
    real32_array_id,     //!< real32_array_id
    real64_array_id,     //!< real64_array_id
    r1_array_id,         //!< r1_array_id
    r2_array_id,         //!< r2_array_id
    r2Sym_array_id,      //!< r2Sym_array_id

    integer_array2d_id,    //!< integer_array2d_id
    localIndex_array2d_id, //!< localIndex_array2d_id
    globalIndex_array2d_id,//!< globalIndex_array2d_id
    real32_array2d_id,     //!< real32_array2d_id
    real64_array2d_id,     //!< real64_array2d_id
    real64_array2d_ji_id,   //!< real64_array2d_ji_id
    r1_array2d_id,         //!< r1_array2d_id
    r2_array2d_id,         //!< r2_array2d_id
    r2Sym_array2d_id,      //!< r2Sym_array2d_id

    integer_array3d_id,    //!< integer_array3d_id
    localIndex_array3d_id, //!< localIndex_array3d_id
    globalIndex_array3d_id,//!< globalIndex_array3d_id
    real32_array3d_id,     //!< real32_array3d_id
    real64_array3d_id,     //!< real64_array3d_id
    real64_array3d_kji_id,  //!< real64_array3d_kji_id

    string_id,           //!< string_id
    Path_id,             //!< Path_id
    string_array_id,     //!< string_array_id
    path_array_id,       //!< path_array_Iid
    mapPair_array_id,    //!< mapPair_array_id
    none_id              //!< none_id
  };

  /**
   * @brief Return a TypeID value given a name.
   * @param name the string of the type
   * @return a TypeIDs value corresponding to the input string
   */
  static TypeIDs typeID( string const & name )
  {
    const std::unordered_map< string, TypeIDs > type_names =
    {
      { "integer", TypeIDs::integer_id },
      { "localIndex", TypeIDs::localIndex_id },
      { "globalIndex", TypeIDs::globalIndex_id },
      { "real32", TypeIDs::real32_id },
      { "real64", TypeIDs::real64_id },
      { "R1Tensor", TypeIDs::r1Tensor_id },
      { "R2Tensor", TypeIDs::r2Tensor_id },
      { "R2SymTensor", TypeIDs::r2SymTensor_id },
      { "integer_array", TypeIDs::integer_array_id },
      { "localIndex_array", TypeIDs::localIndex_array_id },
      { "globalIndex_array", TypeIDs::globalIndex_array_id },
      { "real32_array", TypeIDs::real32_array_id },
      { "real64_array", TypeIDs::real64_array_id },
      { "r1_array", TypeIDs::r1_array_id },
      { "r2_array", TypeIDs::r2_array_id },
      { "r2Sym_array", TypeIDs::r2Sym_array_id },

      { "integer_array2d", TypeIDs::integer_array2d_id },
      { "localIndex_array2d", TypeIDs::localIndex_array2d_id },
      { "globalIndex_array2d", TypeIDs::globalIndex_array2d_id },
      { "real32_array2d", TypeIDs::real32_array2d_id },
      { "real64_array2d", TypeIDs::real64_array2d_id },
      { "real64_array2d_ji", TypeIDs::real64_array2d_ji_id },
      { "r1_array2d", TypeIDs::r1_array2d_id },
      { "r2_array2d", TypeIDs::r2_array2d_id },
      { "r2Sym_array2d", TypeIDs::r2Sym_array2d_id },

      { "integer_array3d", TypeIDs::integer_array3d_id },
      { "localIndex_array3d", TypeIDs::localIndex_array3d_id },
      { "globalIndex_array3d", TypeIDs::globalIndex_array3d_id },
      { "real32_array3d", TypeIDs::real32_array3d_id },
      { "real64_array3d", TypeIDs::real64_array3d_id },
      { "real64_array3d_kji", TypeIDs::real64_array3d_kji_id },

      { "string", TypeIDs::string_id },
      { "path", TypeIDs::Path_id },
      { "string_array", TypeIDs::string_array_id },
      { "path_array", TypeIDs::path_array_id },
      { "map_array", TypeIDs::path_array_id },
      { "mapPair_array", TypeIDs::mapPair_array_id },
      { "", TypeIDs::none_id }
    };
    return type_names.at( name );
  }

  /**
   * @brief Return a TypeID enum given a std::type_index.
   * @param typeIndex the type_index we would to get the TypeID for
   * @return the TypeID associated with the typeIndex
   */
  static TypeIDs typeID( std::type_index typeIndex )
  {
    const std::unordered_map< std::type_index, TypeIDs > type_names =
    {
      { std::type_index( typeid(integer)), TypeIDs::integer_id },
      { std::type_index( typeid(localIndex)), TypeIDs::real32_id },
      { std::type_index( typeid(globalIndex)), TypeIDs::real64_id },
      { std::type_index( typeid(real32)), TypeIDs::real32_id },
      { std::type_index( typeid(real64)), TypeIDs::real64_id },
      { std::type_index( typeid(R1Tensor)), TypeIDs::r1Tensor_id },
      { std::type_index( typeid(R2Tensor)), TypeIDs::r2Tensor_id },
      { std::type_index( typeid(R2SymTensor)), TypeIDs::r2SymTensor_id },
      { std::type_index( typeid(integer_array)), TypeIDs::integer_array_id },
      { std::type_index( typeid(localIndex_array)), TypeIDs::localIndex_array_id },
      { std::type_index( typeid(globalIndex_array)), TypeIDs::globalIndex_array_id },
      { std::type_index( typeid(real32_array)), TypeIDs::real32_array_id },
      { std::type_index( typeid(real64_array)), TypeIDs::real64_array_id },
      { std::type_index( typeid(r1_array)), TypeIDs::r1_array_id },
      { std::type_index( typeid(r2_array)), TypeIDs::r2_array_id },
      { std::type_index( typeid(r2Sym_array)), TypeIDs::r2Sym_array_id },

      { std::type_index( typeid(integer_array2d)), TypeIDs::integer_array2d_id },
      { std::type_index( typeid(localIndex_array2d)), TypeIDs::localIndex_array2d_id },
      { std::type_index( typeid(globalIndex_array2d)), TypeIDs::globalIndex_array2d_id },
      { std::type_index( typeid(real32_array2d)), TypeIDs::real32_array2d_id },
      { std::type_index( typeid(real64_array2d)), TypeIDs::real64_array2d_id },
      { std::type_index( typeid(array2d< real64, RAJA::PERM_JI >)), TypeIDs::real64_array2d_ji_id },
      { std::type_index( typeid(r1_array2d)), TypeIDs::r1_array2d_id },
      { std::type_index( typeid(r2_array2d)), TypeIDs::r2_array2d_id },
      { std::type_index( typeid(r2Sym_array2d)), TypeIDs::r2Sym_array2d_id },

      { std::type_index( typeid(integer_array3d)), TypeIDs::integer_array3d_id },
      { std::type_index( typeid(localIndex_array3d)), TypeIDs::localIndex_array3d_id },
      { std::type_index( typeid(globalIndex_array3d)), TypeIDs::globalIndex_array3d_id },
      { std::type_index( typeid(real32_array3d)), TypeIDs::real32_array3d_id },
      { std::type_index( typeid(real64_array3d)), TypeIDs::real64_array3d_id },
      { std::type_index( typeid(array3d< real64, RAJA::PERM_KJI >)), TypeIDs::real64_array3d_kji_id },

      { std::type_index( typeid(string)), TypeIDs::string_id },
      { std::type_index( typeid(Path)), TypeIDs::Path_id },
      { std::type_index( typeid(string_array)), TypeIDs::string_array_id },
      { std::type_index( typeid(path_array)), TypeIDs::path_array_id },
      { std::type_index( typeid(mapPair_array)), TypeIDs::mapPair_array_id }
    };
    auto iterType = type_names.find( typeIndex );
    if( iterType != type_names.end() )
    {
      return type_names.at( typeIndex );
    }
    else
    {
      return TypeIDs::none_id;
    }
  }

  /**
   * @brief Matching regex for data types in xml.
   */
  class typeRegex
  {
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
    std::string constructArrayRegex( std::string subPattern, integer dimension )
    {
      if( dimension > 1 )
      {
        subPattern = constructArrayRegex( subPattern, dimension-1 );
      }

      std::string arrayPattern;
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
    std::string ru = "[\\d]+";

    // Regex to match an signed int (-123, 455, +789, etc.)
    std::string ri = "[+-]?[\\d]+";

    // Regex to match a float (1, +2.3, -.4, 5.6e7, 8E-9, etc.)
    // Explanation of parts:
    // [+-]?[\\d]*  matches an optional +/- at the beginning, any numbers preceding the decimal
    // ([\\d]\\.?|\\.[\\d]) matches the decimal region of the number (0, 1., 2.3, .4)
    // [\\d]*  matches any number of numbers following the decimal
    // ([eE][-+]?[\\d]+|\\s*)  matches an optional scientific notation number
    // Note: the xsd regex implementation does not allow an empty branch, so use allow whitespace at the end
    std::string rr = "[+-]?[\\d]*([\\d]\\.?|\\.[\\d])[\\d]*([eE][-+]?[\\d]+|\\s*)";

    // Regex to match a string that does not contain the characters  ,{}
    std::string rs = "[^,\\{\\}]*";

    // Regexes to match a R1Tensor, R2Tensor, and R2SymTensor
    // These are identical aside from the number of repetitions in the curly brackets
    std::string r1 = "\\s*(" + rr + ",\\s*){2}" + rr;
    std::string r2 = "\\s*(" + rr + ",\\s*){8}" + rr;
    std::string r2s = "\\s*(" + rr + ",\\s*){5}" + rr;

    // Build master list of regexes
    std::unordered_map< std::string, std::string > regexMap =
    {
      {"integer", ri},
      {"localIndex", ri},
      {"globalIndex", ri},
      {"real32", rr},
      {"real64", rr},
      {"R1Tensor", r1},
      {"R2Tensor", r2},
      {"R2SymTensor", r2s},
      {"integer_array", constructArrayRegex( ri, 1 )},
      {"localIndex_array", constructArrayRegex( ri, 1 )},
      {"globalIndex_array", constructArrayRegex( ri, 1 )},
      {"real32_array", constructArrayRegex( rr, 1 )},
      {"real64_array", constructArrayRegex( rr, 1 )},
      {"r1_array", constructArrayRegex( r1, 1 )},
      {"r2_array", constructArrayRegex( r2, 1 )},
      {"r2Sym_array", constructArrayRegex( r2s, 1 )},
      {"integer_array2d", constructArrayRegex( ri, 2 )},
      {"localIndex_array2d", constructArrayRegex( ri, 2 )},
      {"globalIndex_array2d", constructArrayRegex( ri, 2 )},
      {"real32_array2d", constructArrayRegex( rr, 2 )},
      {"real64_array2d", constructArrayRegex( rr, 2 )},
      {"r1_array2d", constructArrayRegex( r1, 2 )},
      {"r2_array2d", constructArrayRegex( r2, 2 )},
      {"r2Sym_array2d", constructArrayRegex( r2s, 2 )},
      {"integer_array3d", constructArrayRegex( ri, 3 )},
      {"localIndex_array3d", constructArrayRegex( ri, 3 )},
      {"globalIndex_array3d", constructArrayRegex( ri, 3 )},
      {"real32_array3d", constructArrayRegex( rr, 3 )},
      {"real64_array3d", constructArrayRegex( rr, 3 )},
      {"string", rs},
      {"path", rs},
      {"string_array", constructArrayRegex( rs, 1 )},
      {"path_array", constructArrayRegex( rs, 1 )},
      {"mapPair", rs},
      {"mapPair_array", constructArrayRegex( rs, 1 )}
    };

public:

    /**
     * @brief Get an iterator to the beginning of regex map.
     * @return
     */
    std::unordered_map< std::string, std::string >::iterator begin(){return regexMap.begin();}

    /**
     * @brief Get an iterator to the end of regex map.
     * @return
     */
    std::unordered_map< std::string, std::string >::iterator end(){return regexMap.end();}

    /**
     * @brief Get a const iterator to the beginning of regex map.
     * @return
     */
    std::unordered_map< std::string, std::string >::const_iterator begin() const {return regexMap.begin();}

    /**
     * @brief Get a const iterator to the end of regex map.
     * @return
     */
    std::unordered_map< std::string, std::string >::const_iterator end() const {return regexMap.end();}
  };


  /**
   * @brief this function provides a switchyard for the intrinsic supported GEOSX types which calls a generic lambda
   *        that takes in a single argument which may be used to infer type.
   * @tparam LAMBDA the template arg that represents the lambda function
   * @param type the TypeIDs we would like to pass to the lambda function
   * @param lambda the lambda function to call
   * @return the return type of lambda
   */
  template< typename LAMBDA >
  static auto ApplyIntrinsicTypeLambda1( const TypeIDs type,
                                         LAMBDA lambda )
  {
    switch( type )
    {
      case ( TypeIDs::integer_id ):
      {
        return lambda( integer( 1 ) );
      }
      case ( TypeIDs::real32_id ):
      {
        return lambda( real32( 1 ) );
      }
      case ( TypeIDs::real64_id ):
      {
        return lambda( real64( 1 ) );
      }
      case ( TypeIDs::r1Tensor_id ):
      {
        return lambda( R1Tensor() );
      }
      case ( TypeIDs::r2Tensor_id ):
      {
        return lambda( R2Tensor() );
      }
      case ( TypeIDs::r2SymTensor_id ):
      {
        return lambda( R2SymTensor() );
      }
      case ( TypeIDs::integer_array_id ):
      {
        return lambda( integer_array( 1 ) );
      }
      case ( TypeIDs::real32_array_id ):
      {
        return lambda( real32_array( 1 ) );
      }
      case ( TypeIDs::real64_array_id ):
      {
        return lambda( real64_array( 1 ) );
      }
      case ( TypeIDs::string_id ):
      {
        return lambda( string( "" ) );
      }
      case ( TypeIDs::Path_id ):
      {
        return lambda( Path( "" ) );
      }
      default:
      {
        GEOSX_ERROR( "TypeID not recognized." );
      }
    }
  }



  /**
   * @brief this function provides a switchyard for the intrinsic supported GEOSX array types which calls a generic
   *        lambda that takes in a single argument which may be used to infer type.
   * @tparam LAMBDA the template arg that represents the lambda function
   * @param type the TypeIDs we would like to pass to the lambda function
   * @param lambda the lambda function to call
   * @return the return type of lambda
   */
  template< typename LAMBDA >
  static auto ApplyArrayTypeLambda1( const TypeIDs type,
                                     LAMBDA lambda )
  {
    switch( type )
    {
      case ( TypeIDs::integer_array_id ):
      {
        return lambda( integer_array( 1 ) );
      }
      case ( TypeIDs::localIndex_array_id ):
      {
        return lambda( localIndex_array( 1 ) );
      }
      case ( TypeIDs::globalIndex_array_id ):
      {
        return lambda( globalIndex_array( 1 ) );
      }
      case ( TypeIDs::real32_array_id ):
      {
        return lambda( real32_array( 1 ) );
      }
      case ( TypeIDs::real64_array_id ):
      {
        return lambda( real64_array( 1 ) );
      }
      case ( TypeIDs::r1_array_id ):
      {
        return lambda( r1_array( 1 ) );
      }
      case ( TypeIDs::r2_array_id ):
      {
        return lambda( r2_array( 1 ) );
      }
      case ( TypeIDs::r2Sym_array_id ):
      {
        return lambda( r2Sym_array( 1 ) );
      }
      case ( TypeIDs::real64_array2d_id ):
      {
        return lambda( array2d< real64 > {} );
      }
      case ( TypeIDs::real64_array2d_ji_id ):
      {
        return lambda( array2d< real64, RAJA::PERM_JI > {} );
      }
      default:
      {
        GEOSX_ERROR( "TypeID not recognized." );
      }
    }
  }



  /**
   * @brief Provides a switchyard for the intrinsic supported GEOSX array types which calls a generic
   *        lambda that takes in a two arguments argument which may be used to infer array type and underlying type.
   * @tparam LAMBDA the template arg that represents the lambda function
   * @param type the TypeIDs we would like to pass to the lambda function
   * @param errorIfTypeNotFound whether to report an error if the type has not been handled
   * @param lambda the lambda function to call
   * @return the return type of lambda
   */
  template< typename LAMBDA >
  static auto ApplyArrayTypeLambda2( const TypeIDs type,
                                     bool const errorIfTypeNotFound,
                                     LAMBDA && lambda )
  {
    switch( type )
    {
      case ( TypeIDs::integer_array_id ):
      {
        return lambda( integer_array( 1 ), integer( 1 ) );
      }
      case ( TypeIDs::localIndex_array_id ):
      {
        return lambda( localIndex_array( 1 ), localIndex( 1 ) );
      }
      case ( TypeIDs::globalIndex_array_id ):
      {
        return lambda( globalIndex_array( 1 ), globalIndex() );
      }
      case ( TypeIDs::real32_array_id ):
      {
        return lambda( real32_array( 1 ), real32( 1 ) );
      }
      case ( TypeIDs::real64_array_id ):
      {
        return lambda( real64_array( 1 ), real64( 1 ) );
      }
      case ( TypeIDs::r1_array_id ):
      {
        return lambda( r1_array( 1 ), R1Tensor() );
      }
      case ( TypeIDs::r2_array_id ):
      {
        return lambda( r2_array( 1 ), R2Tensor() );
      }
      case ( TypeIDs::r2Sym_array_id ):
      {
        return lambda( r2Sym_array( 1 ), R2SymTensor()  );
      }
      case ( TypeIDs::integer_array2d_id ):
      {
        return lambda( integer_array2d(), integer( 1 ) );
      }
      case ( TypeIDs::localIndex_array2d_id ):
      {
        return lambda( localIndex_array2d(), localIndex( 1 ) );
      }
      case ( TypeIDs::globalIndex_array2d_id ):
      {
        return lambda( globalIndex_array2d(), globalIndex() );
      }
      case ( TypeIDs::real32_array2d_id ):
      {
        return lambda( real32_array2d(), real32( 1 ) );
      }
      case ( TypeIDs::real64_array2d_id ):
      {
        return lambda( real64_array2d(), real64( 1 ) );
      }
      case ( TypeIDs::real64_array2d_ji_id ):
      {
        return lambda( array2d< real64, RAJA::PERM_JI >(), real64( 1 ) );
      }
      case ( TypeIDs::r1_array2d_id ):
      {
        return lambda( r1_array2d(), R1Tensor() );
      }
      case ( TypeIDs::r2_array2d_id ):
      {
        return lambda( r2_array2d(), R2Tensor() );
      }
      case ( TypeIDs::r2Sym_array2d_id ):
      {
        return lambda( r2Sym_array2d(), R2SymTensor()  );
      }
      case ( TypeIDs::integer_array3d_id ):
      {
        return lambda( integer_array3d(), integer( 1 ) );
      }
      case ( TypeIDs::localIndex_array3d_id ):
      {
        return lambda( localIndex_array3d(), localIndex( 1 ) );
      }
      case ( TypeIDs::globalIndex_array3d_id ):
      {
        return lambda( globalIndex_array3d(), globalIndex() );
      }
      case ( TypeIDs::real32_array3d_id ):
      {
        return lambda( real32_array3d(), real32( 1 ) );
      }
      case ( TypeIDs::real64_array3d_id ):
      {
        return lambda( real64_array3d(), real64( 1 ) );
      }
      case ( TypeIDs::real64_array3d_kji_id ):
      {
        return lambda( array3d< real64, RAJA::PERM_KJI >(), real64( 1 ) );
      }
      default:
      {
        if( errorIfTypeNotFound )
        {
          GEOSX_ERROR( "TypeID not recognized." );
        }
      }
    }
  }

  /**
   * @brief this function provides a switchyard for the supported GEOSX types which calls a generic lambda
   *        that takes in a single argument which may be used to infer type.
   * @tparam LAMBDA the template arg that represents the lambda function
   * @param type the TypeIDs we would like to pass to the lambda function
   * @param lambda the lambda function to call
   * @return the return type of lambda
   */
  template< typename LAMBDA >
  static auto ApplyTypeLambda1( const TypeIDs type,
                                LAMBDA lambda )
  {
    switch( type )
    {
      case ( TypeIDs::integer_id ):
      {
        return lambda( integer( 1 ) );
      }
      case ( TypeIDs::real32_id ):
      {
        return lambda( real32( 1 ) );
      }
      case ( TypeIDs::real64_id ):
      {
        return lambda( real64( 1 ) );
      }
      case ( TypeIDs::r1Tensor_id ):
      {
        return lambda( R1Tensor() );
      }
      case ( TypeIDs::r2Tensor_id ):
      {
        return lambda( R2Tensor() );
      }
      case ( TypeIDs::r2SymTensor_id ):
      {
        return lambda( R2SymTensor() );
      }
      case ( TypeIDs::integer_array_id ):
      {
        return lambda( integer_array( 1 ) );
      }
      case ( TypeIDs::localIndex_array_id ):
      {
        return lambda( localIndex_array( 1 ) );
      }
      case ( TypeIDs::globalIndex_array_id ):
      {
        return lambda( globalIndex_array( 1 ) );
      }
      case ( TypeIDs::real32_array_id ):
      {
        return lambda( real32_array( 1 ) );
      }
      case ( TypeIDs::real64_array_id ):
      {
        return lambda( real64_array( 1 ) );
      }
      case ( TypeIDs::r1_array_id ):
      {
        return lambda( r1_array( 1 ) );
      }
      case ( TypeIDs::r2_array_id ):
      {
        return lambda( r2_array( 1 ) );
      }
      case ( TypeIDs::r2Sym_array_id ):
      {
        return lambda( r2Sym_array( 1 ) );
      }
      case ( TypeIDs::integer_array2d_id ):
      {
        return lambda( integer_array2d( 1, 1 ) );
      }
      case ( TypeIDs::localIndex_array2d_id ):
      {
        return lambda( localIndex_array2d( 1, 1 ) );
      }
      case ( TypeIDs::globalIndex_array2d_id ):
      {
        return lambda( globalIndex_array2d( 1, 1 ) );
      }
      case ( TypeIDs::real32_array2d_id ):
      {
        return lambda( real32_array2d( 1, 1 ) );
      }
      case ( TypeIDs::real64_array2d_id ):
      {
        return lambda( real64_array2d( 1, 1 ) );
      }
      case ( TypeIDs::real64_array2d_ji_id ):
      {
        return lambda( array2d< real64, RAJA::PERM_JI >( 1, 1 ) );
      }
      case ( TypeIDs::r1_array2d_id ):
      {
        return lambda( r1_array2d( 1, 1 ) );
      }
      case ( TypeIDs::r2_array2d_id ):
      {
        return lambda( r2_array2d( 1, 1 ) );
      }
      case ( TypeIDs::r2Sym_array2d_id ):
      {
        return lambda( r2Sym_array2d( 1, 1 ) );
      }
      case ( TypeIDs::integer_array3d_id ):
      {
        return lambda( integer_array3d( 1, 1, 1 ) );
      }
      case ( TypeIDs::localIndex_array3d_id ):
      {
        return lambda( localIndex_array3d( 1, 1, 1 ) );
      }
      case ( TypeIDs::globalIndex_array3d_id ):
      {
        return lambda( globalIndex_array3d( 1, 1, 1 ) );
      }
      case ( TypeIDs::real32_array3d_id ):
      {
        return lambda( real32_array3d( 1, 1, 1 ) );
      }
      case ( TypeIDs::real64_array3d_id ):
      {
        return lambda( real64_array3d( 1, 1, 1 ) );
      }
      case ( TypeIDs::real64_array3d_kji_id ):
      {
        return lambda( array3d< real64, RAJA::PERM_KJI >( 1, 1, 1 ) );
      }
      case ( TypeIDs::string_id ):
      {
        return lambda( string( "" ) );
      }
      case ( TypeIDs::Path_id ):
      {
        return lambda( Path( "" ) );
      }
      case ( TypeIDs::string_array_id ):
      {
        return lambda( string_array( 1 ) );
      }
      case ( TypeIDs::path_array_id ):
      {
        return lambda( path_array( 1 ) );
      }
      case ( TypeIDs::mapPair_array_id ):
      {
        return lambda( mapPair_array() );
      }
      default:
      {
        GEOSX_ERROR( "TypeID not recognized." );
        return lambda( double(1) );
      }
    }
  }


  /**
   * @brief this function provides a switchyard for the intrinsic supported GEOSX types which calls a generic
   *        lambda that takes in a two arguments argument which may be used to infer array type and underlying type.
   * @tparam LAMBDA the template arg that represents the lambda function
   * @param type the TypeIDs we would like to pass to the lambda function
   * @param lambda the lambda function to call
   * @return the return type of lambda
   */
  template< typename LAMBDA >
  static auto ApplyTypeLambda2( const TypeIDs type,
                                LAMBDA lambda )
  {
    switch( type )
    {
      case ( TypeIDs::integer_id ):
      {
        return lambda( integer( 1 ), integer( 1 ) );
      }
      case ( TypeIDs::real32_id ):
      {
        return lambda( real32( 1 ), real32( 1 ) );
      }
      case ( TypeIDs::real64_id ):
      {
        return lambda( real64( 1 ), real64( 1 ) );
      }
      case ( TypeIDs::r1Tensor_id ):
      {
        return lambda( R1Tensor(), R1Tensor() );
      }
      case ( TypeIDs::r2Tensor_id ):
      {
        return lambda( R2Tensor(), R2Tensor() );
      }
      case ( TypeIDs::r2SymTensor_id ):
      {
        return lambda( R2SymTensor(), R2SymTensor() );
      }
      case ( TypeIDs::integer_array_id ):
      {
        return lambda( integer_array( 1 ), integer( 1 ) );
      }
      case ( TypeIDs::localIndex_array_id ):
      {
        return lambda( localIndex_array( 1 ), localIndex( 1 ) );
      }
      case ( TypeIDs::globalIndex_array_id ):
      {
        return lambda( globalIndex_array( 1 ), globalIndex( 1 ) );
      }
      case ( TypeIDs::real32_array_id ):
      {
        return lambda( real32_array( 1 ), real32( 1 ) );
      }
      case ( TypeIDs::real64_array_id ):
      {
        return lambda( real64_array( 1 ), real64( 1 ) );
      }
      case ( TypeIDs::r1_array_id ):
      {
        return lambda( r1_array( 1 ), R1Tensor() );
      }
      case ( TypeIDs::r2_array_id ):
      {
        return lambda( r2_array( 1 ), R2Tensor() );
      }
      case ( TypeIDs::r2Sym_array_id ):
      {
        return lambda( r2Sym_array( 1 ), R2SymTensor() );
      }
      case ( TypeIDs::string_id ):
      {
        return lambda( string( "" ), string( "" ) );
      }
      case ( TypeIDs::Path_id ):
      {
        return lambda( Path( "" ), Path( "" ) );
      }
      case ( TypeIDs::string_array_id ):
      {
        return lambda( string_array( 1 ), string( "" ) );
      }
      case ( TypeIDs::path_array_id ):
      {
        return lambda( path_array( 1 ), Path( "" ) );
      }
      // case ( TypeIDs::mapPair_array_id ):
      // {
      //   return lambda( mapPair_array(1), mapPair({}) );
      // }
      default:
      {
        GEOSX_ERROR( "TypeID not recognized." );
      }
    }
  }


  /**
   * @brief this function provides a switchyard for the supported GEOSX types which calls a generic
   *        lambda that takes in a two arguments argument which may be used to infer array type and underlying type.
   * @tparam LAMBDA the template arg that represents the lambda function
   * @param type the TypeIDs we would like to pass to the lambda function
   * @param lambda the lambda function to call
   * @return the return type of lambda
   */
  template< typename LAMBDA >
  static auto ApplyIntrinsicTypeLambda2( const TypeIDs type,
                                         LAMBDA lambda )
  {
    switch( type )
    {
      case ( TypeIDs::integer_id ):
      {
        return lambda( integer( 1 ), integer( 1 ) );
      }
      case ( TypeIDs::real32_id ):
      {
        return lambda( real32( 1 ), real32( 1 ) );
      }
      case ( TypeIDs::real64_id ):
      {
        return lambda( real64( 1 ), real64( 1 ) );
      }
      case ( TypeIDs::r1Tensor_id ):
      {
        return lambda( R1Tensor(), R1Tensor() );
      }
      case ( TypeIDs::r2Tensor_id ):
      {
        return lambda( R2Tensor(), R2Tensor() );
      }
      case ( TypeIDs::r2SymTensor_id ):
      {
        return lambda( R2SymTensor(), R2SymTensor() );
      }
      case ( TypeIDs::integer_array_id ):
      {
        return lambda( integer_array( 1 ), integer( 1 ) );
      }
      case ( TypeIDs::real32_array_id ):
      {
        return lambda( real32_array( 1 ), real32( 1 ) );
      }
      case ( TypeIDs::real64_array_id ):
      {
        return lambda( real64_array( 1 ), real64( 1 ) );
      }
      case ( TypeIDs::string_id ):
      {
        return lambda( string( "" ), string( "" ) );
      }
      case ( TypeIDs::Path_id ):
      {
        return lambda( Path( "" ), Path( "" ) );
      }
      case ( TypeIDs::string_array_id ):
      {
        return lambda( string_array( 1 ), string( "" ) );
      }
      case ( TypeIDs::path_array_id ):
      {
        return lambda( path_array( 1 ), Path( "" ) );
      }
      case ( TypeIDs::integer_array2d_id ):
      {
        return lambda( integer_array2d(), integer( 1 ) );
      }
      case ( TypeIDs::localIndex_array2d_id ):
      {
        return lambda( localIndex_array2d(), localIndex( 1 ) );
      }
      case ( TypeIDs::globalIndex_array2d_id ):
      {
        return lambda( globalIndex_array2d(), globalIndex() );
      }
      case ( TypeIDs::real32_array2d_id ):
      {
        return lambda( real32_array2d(), real32( 1 ) );
      }
      case ( TypeIDs::real64_array2d_id ):
      {
        return lambda( real64_array2d(), real64( 1 ) );
      }
      case ( TypeIDs::real64_array2d_ji_id ):
      {
        return lambda( array2d< real64, RAJA::PERM_JI >(), real64( 1 ) );
      }
      //  case ( TypeIDs::mapPair_array_id ):
      //  {
      //    return lambda( mapPair_array(1), mapPair({}) );
      //  }
      default:
      {
        GEOSX_ERROR( "TypeID not recognized." );
      }
    }
  }

};

}



#endif /* GEOSX_COMMON_DATATYPES_HPP */
