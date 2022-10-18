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
 * @file TypeDispatch.hpp
 *
 * Collection of tools to help dispatch templated functions on types
 */

#ifndef GEOSX_COMMON_TYPEDISPATCH_HPP
#define GEOSX_COMMON_TYPEDISPATCH_HPP

#include "common/DataTypes.hpp"

#include <camp/camp.hpp>

#include <unordered_map>

namespace geosx
{

/**
 * @brief Namespace containing type dispatching facilities.
 */
namespace types
{

namespace internal
{

template< camp::idx_t NDIM >
struct LayoutList
{
  static_assert( NDIM < 5, "Dispatching on 5D arrays and higher not implemented" );
  using types = camp::list<>;
};

template<>
struct LayoutList< 1 >
{
  using types = camp::list< RAJA::PERM_I >;
};

template<>
struct LayoutList< 2 >
{
  using types = camp::list< RAJA::PERM_IJ,
                            RAJA::PERM_JI >;
};

template<>
struct LayoutList< 3 >
{
  using types = camp::list< RAJA::PERM_IJK,
                            RAJA::PERM_IKJ,
                            RAJA::PERM_JIK,
                            RAJA::PERM_JKI,
                            RAJA::PERM_KIJ,
                            RAJA::PERM_KJI >;
};

template<>
struct LayoutList< 4 >
{
  using types = camp::list< RAJA::PERM_IJKL,
                            RAJA::PERM_IJLK,
                            RAJA::PERM_IKJL,
                            RAJA::PERM_IKLJ,
                            RAJA::PERM_ILJK,
                            RAJA::PERM_ILKJ,
                            RAJA::PERM_JIKL,
                            RAJA::PERM_JILK,
                            RAJA::PERM_JKIL,
                            RAJA::PERM_JKLI,
                            RAJA::PERM_JLIK,
                            RAJA::PERM_JLKI,
                            RAJA::PERM_KIJL,
                            RAJA::PERM_KILJ,
                            RAJA::PERM_KJIL,
                            RAJA::PERM_KJLI,
                            RAJA::PERM_KLIJ,
                            RAJA::PERM_KLJI,
                            RAJA::PERM_LIJK,
                            RAJA::PERM_LIKJ,
                            RAJA::PERM_LJIK,
                            RAJA::PERM_LJKI,
                            RAJA::PERM_LKIJ,
                            RAJA::PERM_LKJI >;
};

template< int NDIM >
using Layouts = typename LayoutList< NDIM >::types;

// Expands a list of dimensions with all corresponding permutations
template< typename T >
struct AddLayoutsImpl;

template< camp::idx_t ... NDIMS >
struct AddLayoutsImpl< camp::list< camp::num< NDIMS >... > >
{
  using types = typename camp::join< camp::cartesian_product< camp::list< camp::num< NDIMS > >, Layouts< NDIMS > > ... >::type;
};

template< typename T >
using AddLayouts = typename AddLayoutsImpl< T >::types;

// Helper to convert a constructed camp::list of type tuples into an actual Array type
template< typename T >
struct ArrayType;

template< typename T, typename NDIM, typename LAYOUT >
struct ArrayType< camp::list< T, camp::list< NDIM, LAYOUT > > >
{
  using type = ::geosx::Array< T, NDIM::value, LAYOUT >;
};

// Helper to apply a template to all types in a list
template< template< typename > class F, typename T >
struct ApplyImpl;

template< template< typename > class F, typename ... Ts >
struct ApplyImpl< F, camp::list< Ts ... > >
{
  using types = camp::list< typename F< Ts >::type ... >;
};

template< template< typename > class F, typename T >
using Apply = typename ApplyImpl< F, T >::types;

// Helper to increment values in an compile-time integer sequence
template< typename T, int N >
struct IncrementImpl;

template< camp::idx_t ... Is, int N >
struct IncrementImpl< camp::idx_seq< Is... >, N >
{
  using type = camp::idx_seq< Is + N ... >;
};

template< typename T, int N >
using Increment = typename IncrementImpl< T, N >::type;

} // namespace detail

/**
 * @brief Construct a list of types.
 */
template< typename ... Ts >
using TypeList = camp::list< Ts... >;

/**
 * @brief Concatenate multiple type lists.
 */
template< typename ... Ls >
using Join = typename camp::join< Ls ... >::type;

/**
 * @brief Construct a list of GEOSX multidimensional array types (geosx::Array), containing all
 *        value types in type list @p TYPES and all dimensionalities in list of integral constants @p NDIMS.
 */
template< typename TYPES, typename NDIMS >
using ArrayTypes = internal::Apply< internal::ArrayType, camp::cartesian_product< TYPES, internal::AddLayouts< NDIMS > > >;

/**
 * @brief List of major integral types used in GEOSX type system.
 */
using IntegralTypes = TypeList< integer, localIndex, globalIndex >;

/**
 * @brief List of major real-valued (floating point) types used in GEOSX type system.
 */
using RealTypes = TypeList< real32, real64 >;

/**
 * @brief Generate a list of types representing array dimensionalities from M up to (and including) @p N.
 */
template< int M, int N >
using DimsRange = camp::as_list< internal::Increment< camp::make_idx_seq_t< N - M + 1 >, M > >;

/**
 * @brief Generate a list of types representing array dimensionality exactly @p N.
 */
template< int N >
using DimsSingle = DimsRange< N, N >;

/**
 * @brief Generate a list of types representing array dimensionalities up to (and including) @p N.
 */
template< int N >
using DimsUpTo = DimsRange< 1, N >;

/**
 * @brief List of real-valued array types typically used in GEOSX (dimensions up to 4).
 */
using RealArrays = ArrayTypes< RealTypes, DimsUpTo< 4 > >;

/**
 * @brief List of integer-valued array types typically used in GEOSX (dimensions 1 and 2).
 */
using IntegralArrays = ArrayTypes< IntegralTypes, DimsUpTo< 2 > >;

/**
 * @brief List of all array types (real- and integer-valued) typically used in GEOSX.
 */
using StandardArrays = Join< RealArrays, IntegralArrays >;

namespace internal
{

template< typename ... Ts, std::size_t ... Is >
std::unordered_map< std::type_index, std::size_t > const &
getTypeIndexMap( TypeList< Ts... >,
                 std::index_sequence< Is... > )
{
  static std::unordered_map< std::type_index, std::size_t > const result( { { typeid( Ts ), Is } ... } );
  return result;
}

template< typename ... Ts, typename LAMBDA >
bool dispatchViaTable( TypeList< Ts... > const types,
                       std::type_index const type,
                       LAMBDA && lambda )
{
  static_assert( sizeof...(Ts) > 0, "Dispatching on empty type list not supported" );

  // Initialize a table of handlers, once per unique combination of type list and lambda
  using Handler = void (*)( LAMBDA && );
  static Handler const handlers[] = { []( LAMBDA && f ){ f( Ts{} ); } ... };

  // Initialize a hashmap of std::type_index to contiguous indices, once per unique type list
  auto const & typeIndexMap = getTypeIndexMap( types, std::index_sequence_for< Ts... >{} );
  auto const it = typeIndexMap.find( type );
  if( it != typeIndexMap.end() )
  {
    GEOSX_ASSERT_GT( sizeof...( Ts ), it->second ); // sanity check
    handlers[ it->second ]( std::forward< LAMBDA >( lambda ) );
    return true;
  }
  return false;
}

} // namespace internal

/**
 * @brief Dispatch a generic worker function @p lambda based on runtime type.
 * @tparam Ts list of types
 * @tparam LAMBDA type of user-provided function or lambda, must have one auto argument
 * @param types list of types as an object (for deduction)
 * @param type type index of the runtime type
 * @param errorIfTypeNotFound flag indicating whether dispatch should issue an error when type not matched
 * @param lambda user-provided callable, will be called with a single prvalue of type indicated by @p type
 * @return @p true iff type has been dispatch
 *
 * @todo Do we want @p errorIfTypeNotFound parameter? Options:
 *       - make it a template parameter (caller always knows whether or not it wants hard errors)
 *       - make the caller process return value and raise error if needed (and however they want)
 */
template< typename ... Ts, typename LAMBDA >
bool dispatch( TypeList< Ts... > const types,
               std::type_index const type,
               bool const errorIfTypeNotFound,
               LAMBDA && lambda )
{
  bool const success = internal::dispatchViaTable( types, type, std::forward< LAMBDA >( lambda ) );
  if( !success && errorIfTypeNotFound )
  {
    GEOSX_ERROR( "Type " << LvArray::system::demangle( type.name() ) << " was not dispatched.\n" <<
                 "Check the stack trace below and revise the type list passed to dispatch().\n" <<
                 "If you are unsure about this error, please report it to GEOSX issue tracker." );
  }
  return success;
}

} // namespace types

} // namespace geosx

#endif //GEOSX_COMMON_TYPEDISPATCH_HPP
