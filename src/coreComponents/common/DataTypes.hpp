/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file DataTypes.hpp
 * This file contains various aliases and functions that provide operations regarding the
 * use of the data types.
 */

#ifndef DATATYPES_HPP
#define DATATYPES_HPP

#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>
#include <vector>

#include "common/GeosxConfig.hpp"
#include "Macros.hpp"
#include "Logger.hpp"
#include "Array.hpp"
#include "StackArrayWrapper.hpp"
#include "SortedArray.hpp"
#include "ArrayOfArrays.hpp"
#include "math/TensorT/TensorT.h"

#ifdef GEOSX_USE_ATK
#include "axom/sidre/core/SidreTypes.hpp"
#endif

#ifdef GEOSX_USE_MPI
#include <mpi.h>
#endif


/**
 * top level geosx namespace contains all code that is specific to GEOSX
 */
namespace geosx
{


#ifdef GEOSX_USE_MPI
extern MPI_Comm MPI_COMM_GEOSX;
#endif

// underlying types not for general use!!
//using int32 = std::int32_t;
//using int64 = std::int64_t;
//using uint32 = std::uint32_t;
//using uint64 = std::uint64_t;

/// alias for std::size_t
using size_t      = std::size_t;
using integer     = std::int32_t;
using localIndex  = std::ptrdiff_t;
using globalIndex = long long int;//std::int64_t;
using string      = std::string;

using real32 = float;
using real64 = double;
//using real   = double;


template< typename T >
using ptr = T*;

template< typename T >
using const_ptr = T const *;


using integer_ptr        = ptr< integer >;
using integer_const_ptr  = const_ptr< integer >;

using localIndex_ptr         = ptr< localIndex >;
using localIndex_const_ptr   = const_ptr< localIndex >;

using real32_ptr        = ptr< real32 >;
using real32_const_ptr  = const_ptr< real32 >;

using real64_ptr        = ptr< real64 >;
using real64_const_ptr  = const_ptr< real64 >;


using buffer_unit_type = char;
using buffer_type = std::vector< buffer_unit_type >;

//***** BEGIN ARRAY TYPEDEFS *****

template< typename T, int NDIM=1 >
using array_view = LvArray::ArrayView< T, NDIM, localIndex >;

template< typename T, int NDIM=1 >
using array_slice = LvArray::ArraySlice< T, NDIM, localIndex >;

template< typename T, int NDIM, int MAXSIZE >
using stack_array = LvArray::Array< T, NDIM, localIndex, LvArray::StackArrayWrapper< T, MAXSIZE > >;

template< typename T >
using array1d = LvArray::Array< T, 1, localIndex >;

template< typename T >
using arrayView1d = array_view< T, 1 >;

template< typename T >
using arraySlice1d = LvArray::ArraySlice1d< T, localIndex >;

template< typename T >
using arraySlice1dRval = LvArray::ArraySlice1d_rval< T, localIndex >;

template< typename T, int MAXSIZE >
using stackArray1d = stack_array< T, 1, MAXSIZE >;

template< typename T >
using array2d = LvArray::Array< T, 2, localIndex >;

template< typename T >
using arrayView2d = array_view< T, 2 >;

template< typename T >
using arraySlice2d = LvArray::ArraySlice< T, 2, localIndex >;

template< typename T, int MAXSIZE >
using stackArray2d = stack_array< T, 2, MAXSIZE >;

template< typename T, bool CONST_SIZES=false >
using ArrayOfArraysView = LvArray::ArrayOfArraysView< T, localIndex const, CONST_SIZES >;

template< typename T >
using ArrayOfArrays = LvArray::ArrayOfArrays< T, localIndex >;

template< typename T >
using array3d = LvArray::Array< T, 3, localIndex >;

template< typename T >
using arrayView3d = array_view< T, 3 >;

template< typename T >
using arraySlice3d = LvArray::ArraySlice< T, 3, localIndex >;

template< typename T, int MAXSIZE >
using stackArray3d = stack_array< T, 3, MAXSIZE >;

template< typename T >
using array4d = LvArray::Array< T, 4, localIndex >;

template< typename T >
using arrayView4d = LvArray::ArrayView< T, 4, localIndex >;

template< typename T >
using arraySlice4d = LvArray::ArrayView< T, 4, localIndex >;

template< typename T, int MAXSIZE >
using stackArray4d = stack_array< T, 4, MAXSIZE >;

template< typename T >
using array5d = LvArray::Array< T, 5, localIndex >;

template< typename T >
using arrayView5d = LvArray::ArrayView< T, 5, localIndex >;

template< typename T >
using arraySlice5d = LvArray::ArrayView< T, 5, localIndex >;

template< typename T, int MAXSIZE >
using stackArray5d = stack_array< T, 5, MAXSIZE >;

template< typename T >
using set = LvArray::SortedArray< T, localIndex >;

template< typename TKEY, typename TVAL >
using map = std::map< TKEY, TVAL >;

template< typename TKEY, typename TVAL >
using unordered_map = std::unordered_map< TKEY, TVAL >;

using integer_array        = array1d< integer >;
using integer_const_array  = array1d< integer const >;

using real32_array        = array1d< real32 >;
using real32_const_array  = array1d< real32 const >;

using real64_array        = array1d< real64 >;
using real64_const_array  = array1d< real64 const >;

using string_array        = array1d< string >;
using string_const_array  = array1d< string const >;

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


//***** BEGIN LEGACY TYPEDEFS *****

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

//***** END LEGACY TYPEDEFS *****


void printTypeSummary();
/**
 * @class rtTypes
 * class to manage the type selection of types at runtime
 */
class rtTypes
{
public:

  /**
   * @brief function to return a typename as a string
   * @param key the std::type_index of the type
   * @return a hard coded string taht is related to the std::type_index
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
      {std::type_index( typeid(string_array)), "string_array"},
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
   * @enum TypeIDs
   * @brief a set of enums for each geosx defined data type
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
    r1_array2d_id,         //!< r1_array2d_id
    r2_array2d_id,         //!< r2_array2d_id
    r2Sym_array2d_id,      //!< r2Sym_array2d_id

    integer_array3d_id,    //!< integer_array3d_id
    localIndex_array3d_id, //!< localIndex_array3d_id
    globalIndex_array3d_id,//!< globalIndex_array3d_id
    real32_array3d_id,     //!< real32_array3d_id
    real64_array3d_id,     //!< real64_array3d_id

    string_id,           //!< string_id
    string_array_id,     //!< string_array_id
    mapPair_array_id,    //!< mapPair_array_id
    none_id              //!< none_id
  };

  /**
   * @brief function to return a TypeID value given a name
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
      { "r1_array2d", TypeIDs::r1_array2d_id },
      { "r2_array2d", TypeIDs::r2_array2d_id },
      { "r2Sym_array2d", TypeIDs::r2Sym_array2d_id },

      { "integer_array3d", TypeIDs::integer_array3d_id },
      { "localIndex_array3d", TypeIDs::localIndex_array3d_id },
      { "globalIndex_array3d", TypeIDs::globalIndex_array3d_id },
      { "real32_array3d", TypeIDs::real32_array3d_id },
      { "real64_array3d", TypeIDs::real64_array3d_id },

      { "string", TypeIDs::string_id },
      { "string_array", TypeIDs::string_array_id },
      { "mapPair_array", TypeIDs::mapPair_array_id },
      { "", TypeIDs::none_id }
    };
    return type_names.at( name );
  }

  /**
   * @brief function to return a TypeID enum given a std::type_index
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
      { std::type_index( typeid(r1_array2d)), TypeIDs::r1_array2d_id },
      { std::type_index( typeid(r2_array2d)), TypeIDs::r2_array2d_id },
      { std::type_index( typeid(r2Sym_array2d)), TypeIDs::r2Sym_array2d_id },

      { std::type_index( typeid(integer_array3d)), TypeIDs::integer_array3d_id },
      { std::type_index( typeid(localIndex_array3d)), TypeIDs::localIndex_array3d_id },
      { std::type_index( typeid(globalIndex_array3d)), TypeIDs::globalIndex_array3d_id },
      { std::type_index( typeid(real32_array3d)), TypeIDs::real32_array3d_id },
      { std::type_index( typeid(real64_array3d)), TypeIDs::real64_array3d_id },

      { std::type_index( typeid(string)), TypeIDs::string_id },
      { std::type_index( typeid(string_array)), TypeIDs::string_array_id },
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

#ifdef GEOSX_USE_ATK

  static axom::sidre::TypeID toSidreType( std::type_index typeIndex )
  {
    const axom::sidre::TypeID integer_id = axom::sidre::detail::SidreTT< integer >::id;
    const axom::sidre::TypeID localIndex_id = axom::sidre::detail::SidreTT< localIndex >::id;

    /* We can't use SidreTT<globalIndex>::id here because that returns NO_TYPE_ID.
     * This is due to a mismatch between globalIndex (long long int) and std::int64_t */
    const axom::sidre::TypeID globalIndex_id = axom::sidre::detail::SidreTT< axom::common::int64 >::id;

    const axom::sidre::TypeID real32_id = axom::sidre::detail::SidreTT< real32 >::id;
    const axom::sidre::TypeID real64_id = axom::sidre::detail::SidreTT< real64 >::id;
    const axom::sidre::TypeID char_id = axom::sidre::TypeID::UINT8_ID;

    const std::unordered_map< std::type_index, axom::sidre::TypeID > sidre_types =
    {
      { std::type_index( typeid(integer)), integer_id },
      { std::type_index( typeid(localIndex)), localIndex_id },
      { std::type_index( typeid(globalIndex)), globalIndex_id },
      { std::type_index( typeid(real32)), real32_id },
      { std::type_index( typeid(real64)), real64_id },
      { std::type_index( typeid(R1Tensor)), real64_id },
      { std::type_index( typeid(R2Tensor)), real64_id },
      { std::type_index( typeid(R2SymTensor)), real64_id },
      { std::type_index( typeid(char)), char_id }
    };

    auto it = sidre_types.find( typeIndex );
    if( it == sidre_types.end())
    {
      return axom::sidre::TypeID::NO_TYPE_ID;
    }
    return it->second;
  }

  static localIndex getSidreSize( std::type_index typeIndex )
  {
    const std::unordered_map< std::type_index, localIndex > sidre_sizes =
    {
      { std::type_index( typeid(integer)), sizeof(integer) },
      { std::type_index( typeid(localIndex)), sizeof(localIndex) },
      { std::type_index( typeid(globalIndex)), sizeof(globalIndex) },
      { std::type_index( typeid(real32)), sizeof(real32) },
      { std::type_index( typeid(real64)), sizeof(real64) },
      { std::type_index( typeid(R1Tensor)), sizeof(real64) },
      { std::type_index( typeid(R2Tensor)), sizeof(real64) },
      { std::type_index( typeid(R2SymTensor)), sizeof(real64) },
      { std::type_index( typeid(char)), sizeof(char) }
    };

    auto it = sidre_sizes.find( typeIndex );
    if( it == sidre_sizes.end())
    {
      GEOS_ERROR( "Unsupported type of with type index name: " << typeIndex.name());
    }
    return it->second;
  }

#endif /* GEOSX_USE_ATK */


  // Matching regex for data types in xml
  class typeRegex
  {
private:
    std::string ru = "[0-9]*";
    std::string ri = "[+-]?[0-9]*";
    std::string rr = "[+-]?[0-9]*\\.?([0-9]*)?[eE]?[-+]?([0-9]*)?";
    std::string rs = "[a-zA-Z0-9_,\\(\\)+-/\\*]*";
    std::string r1 = rr + ",? " + rr + ",? " + rr;
    std::string r2 = rr + ",? " + rr + ",? " + rr + ",? " + rr + ",? " + rr + ",? " + rr + ",? " + rr + ",? " + rr + ",? " + rr;
    std::string r2s = rr + ",? " + rr + ",? " + rr + ",? " + rr + ",? " + rr + ",? " + rr;

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
      {"integer_array", "((" + ri + ",? )*)?" + ri},
      {"localIndex_array", "((" + ri + ",? )*)?" + ri},
      {"globalIndex_array", "((" + ri + ",? )*)?" + ri},
      {"real32_array", "((" + rr + ",? )*)?" + rr},
      {"real64_array", "((" + rr + ",? )*)?" + rr},
      {"r1_array", "((" + r1 + "; )*)?" + r1},
      {"r2_array", "((" + r2 + "; )*)?" + r2},
      {"r2Sym_array", "((" + r2s + "; )*)?" + r2s},
      {"integer_array2d", ""},
      {"localIndex_array2d", ""},
      {"globalIndex_array2d", ""},
      {"real32_array2d", ""},
      {"real64_array2d", ""},
      {"r1_array2d", ""},
      {"r2_array2d", ""},
      {"r2Sym_array2d", ""},
      {"integer_array3d", ""},
      {"localIndex_array3d", ""},
      {"globalIndex_array3d", ""},
      {"real32_array3d", ""},
      {"real64_array3d", ""},
      {"string", rs},
      {"string_array", "((" + rs + ",? )*)?" + rs},
      {"mapPair", rs},
      {"mapPair_array", "((" + rs + ",? )*)?" + rs}
    };

public:
    std::unordered_map< std::string, std::string >::iterator begin(){return regexMap.begin();}
    std::unordered_map< std::string, std::string >::iterator end(){return regexMap.end();}
    std::unordered_map< std::string, std::string >::const_iterator begin() const {return regexMap.begin();}
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
      default:
      {
        GEOS_ERROR( "TypeID not recognized." );
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
      default:
      {
        GEOS_ERROR( "TypeID not recognized." );
      }
    }
  }



  /**
   * @brief this function provides a switchyard for the intrinsic supported GEOSX array types which calls a generic
   *        lambda that takes in a two arguments argument which may be used to infer array type and underlying type.
   * @tparam LAMBDA the template arg that represents the lambda function
   * @param type the TypeIDs we would like to pass to the lambda function
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
      default:
      {
        if( errorIfTypeNotFound )
        {
          GEOS_ERROR( "TypeID not recognized." );
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
      case ( TypeIDs::string_id ):
      {
        return lambda( string( "" ) );
      }
      case ( TypeIDs::string_array_id ):
      {
        return lambda( string_array( 1 ) );
      }
      case ( TypeIDs::mapPair_array_id ):
      {
        return lambda( mapPair_array() );
      }
      default:
      {
        GEOS_ERROR( "TypeID not recognized." );
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
      case ( TypeIDs::string_array_id ):
      {
        return lambda( string_array( 1 ), string( "" ) );
      }
      // case ( TypeIDs::mapPair_array_id ):
      // {
      //   return lambda( mapPair_array(1), mapPair({}) );
      // }
      default:
      {
        GEOS_ERROR( "TypeID not recognized." );
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
      case ( TypeIDs::string_array_id ):
      {
        return lambda( string_array( 1 ), string( "" ) );
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
      //  case ( TypeIDs::mapPair_array_id ):
      //  {
      //    return lambda( mapPair_array(1), mapPair({}) );
      //  }
      default:
      {
        GEOS_ERROR( "TypeID not recognized." );
      }
    }
  }

  // TODO depricate
  inline static void equate( R1Tensor & lhs, integer const component, real64 const & rhs )
  {
    lhs[component] = rhs;
  }

  // TODO depricate
  template< typename TLHS, typename TRHS >
  inline static void equate( TLHS & lhs,
                             integer const,//component,
                             TRHS const & rhs )
  {
    lhs = rhs;
  }


  // TODO depricate
  inline static void add( R1Tensor & lhs,
                          integer const component,
                          real64 const & rhs )
  {
    lhs[component] += rhs;
  }

  // TODO depricate
  template< typename TLHS, typename TRHS >
  inline static void add( TLHS & lhs,
                          integer const,// component,
                          TRHS const & rhs )
  {
    lhs += rhs;
  }



  // TODO depricate
  inline static real64 value( R1Tensor & lhs, integer const component )
  {
    return lhs[component];
  }

  // TODO depricate
  inline static real64 value( R2Tensor & lhs, integer const component )
  {
    return lhs.Data()[component];
  }

  // TODO depricate
  inline static real64 value( R2SymTensor & lhs, integer const component )
  {
    return lhs.Data()[component];
  }

  // TODO depricate
  template< typename TLHS >
  inline static TLHS value( TLHS const & lhs, integer const )
  {
    return lhs;
  }

};

}



#endif /* COMPONENTS_CORE_SRC_DATAREPOSITORY_DATATYPES_HPP_ */
