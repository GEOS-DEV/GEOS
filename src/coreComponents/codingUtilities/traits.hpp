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
 * @file traits.hpp
 */

#ifndef GEOSX_CODINGUTILITIES_TRAITS_HPP_
#define GEOSX_CODINGUTILITIES_TRAITS_HPP_

// Source includes
#include "common/DataTypes.hpp"
#include "LvArray/src/typeManipulation.hpp"
#include "LvArray/src/bufferManipulation.hpp"
#include "SFINAE_Macros.hpp"

// System includes
#include <type_traits>

namespace geosx
{

namespace traits
{

/**
 * @brief Defines a static constexpr bool HasMemberFunction_data< @p CLASS >
 *        that is true iff the method @p CLASS ::data() exists and the return value is convertable to a pointer.
 * @tparam CLASS The type to test.
 */
HAS_MEMBER_FUNCTION( data, void const *, );

/**
 * @brief Defines a static constexpr bool HasMemberFunction_move< @p CLASS >
 *        that is true iff the method @p CLASS ::move(LvArray::MemorySpace, bool) exists.
 * @tparam CLASS The type to test.
 */
template< typename CLASS >
static constexpr bool HasMemberFunction_move = LvArray::bufferManipulation::HasMemberFunction_move< CLASS >;

/**
 * @brief Defines a static constexpr bool HasMemberFunction_setName< @p CLASS >
 *        that is true iff the method @p CLASS ::setName( string ) exists.
 * @tparam CLASS The type to test.
 */
HAS_MEMBER_FUNCTION_NO_RTYPE( setName, string() );

/**
 * @brief Defines a static constexpr bool HasMemberFunction_size< @p CLASS >
 *        that is true iff the method @p CLASS ::size() exists and the return value is convertable to a localIndex.
 * @tparam CLASS The type to test.
 */
HAS_MEMBER_FUNCTION( size, localIndex, );

/**
 * @brief Defines a static constexpr bool HasMemberFunction_capacity< @p CLASS >
 *        that is true iff the method @p CLASS ::capacity() exists and the return value is convertable to a localIndex.
 * @tparam CLASS The type to test.
 */
HAS_MEMBER_FUNCTION( capacity, localIndex, );

/**
 * @brief Defines a static constexpr bool HasMemberFunction_resize< @p CLASS >
 *        that is True iff the method @p CLASS ::resize( int ) exists.
 * @tparam CLASS The type to test.
 */
HAS_MEMBER_FUNCTION_NO_RTYPE( resize, 0 );

/**
 * @brief Defines a static constexpr bool HasMemberFunction_reserve< @p CLASS >
 *        that is true iff the method @p CLASS ::reserve( localIndex ) exists.
 * @tparam CLASS The type to test.
 */
HAS_MEMBER_FUNCTION_NO_RTYPE( reserve, localIndex( 55 ) );

/**
 * @brief Defines a static constexpr bool HasMemberFunction_toView< @p CLASS >
 *        that is true iff the method @p CLASS ::toView() exists.
 * @tparam CLASS The type to test.
 */
template< typename CLASS >
static constexpr bool HasMemberFunction_toView = LvArray::typeManipulation::HasMemberFunction_toView< CLASS >;

/**
 * @brief Defines a static constexpr bool with two template parameter CanStreamInto
 *        that is true iff the operation std::declval< SRC & >() >> src::declval< DST & >() exists.
 * @tparam SRC The type of the source.
 * @tparam DST The type of the destination.
 */
IS_VALID_EXPRESSION_2( CanStreamInto, SRC, DST, std::declval< SRC & >() >> std::declval< DST & >() );

/**
 * @brief Defines a static constexpr bool HasMemberType_value_type< @p CLASS >
 *        that is true iff @p CLASS ::value_type is valid and not an enum.
 * @tparam CLASS The type to test.
 */
HAS_ALIAS( value_type );

namespace internal
{
template< class T,
          bool HAS_DATA_METHOD = HasMemberFunction_data< T > >
struct GetPointerType
{
  using Pointer = T *;
  using ConstPointer = T const *;
};

template< class T >
struct GetPointerType< T, true >
{
  using Pointer = decltype( std::declval< T >().data() );
  using ConstPointer = std::remove_pointer_t< Pointer > const *;
};

} // namespace internal

/// Type aliased to whatever T::data() returns or T * if that method doesn't exist.
template< typename T >
using Pointer = typename internal::GetPointerType< T >::Pointer;

/// The const version of Pointer.
template< typename T >
using ConstPointer = typename internal::GetPointerType< T >::ConstPointer;

/// Type aliased to whatever T::toView() returns or T & if that method doesn't exist.
template< typename T >
using ViewType = LvArray::typeManipulation::ViewType< T >;

/// Type aliased to whatever T::toViewConst() returns or T const & if that method doesn't exist.
template< typename T >
using ViewTypeConst = LvArray::typeManipulation::ViewTypeConst< T >;

/// True if T is or inherits from string.
template< typename T >
constexpr bool is_string = std::is_base_of< string, T >::value;

/// True if T is an instantiation of LvArray::Array.
template< typename T >
constexpr bool is_array = LvArray::isArray< T >;

/// True if T is an instantiation of LvArray::ArrayView.
template< typename T >
constexpr bool is_array_view = LvArray::isArrayView< T >;

/// True if T is an instantiation of LvArray::SortedArray.
template< typename T >
constexpr bool is_sorted_array = LvArray::isSortedArray< T >;

/// True if T is an instantiation of LvArray:SortedArrayView.
template< typename T >
constexpr bool is_sorted_array_view = LvArray::isSortedArrayView< T >;

/// True if T is an instantiation of LvArray::ArrayView or LvArray::Array
template< typename T >
constexpr bool is_array_type = traits::is_array_view< T > || traits::is_array< T >;

/// True if T is an instantiation of LvArray::SortedArrayView or LvArray::SortedArray
template< typename T >
constexpr bool is_sorted_array_type = traits::is_sorted_array_view< T > || traits::is_sorted_array< T >;

/// True if T is a Tensor class.
template< typename T >
constexpr bool is_tensorT = std::is_same< std::remove_const_t< T >, R1Tensor >::value;

/// True of T has operator=() defined.
template< typename _T >
struct hasCopyAssignmentOperatorImpl
{
private:
  template< typename T > static constexpr auto test( int )->decltype( T()=T(), bool () )
  { return true; }

  template< typename T > static constexpr auto test( ... )->bool
  { return false; }
public:
  static constexpr bool value = test< _T >( 0 );
};
/// True if T has operator= defined, or it is arithmetic or an enum.
template< typename T >
static constexpr bool hasCopyAssignmentOp = hasCopyAssignmentOperatorImpl< T >::value ||
                                            std::is_arithmetic< T >::value ||
                                            std::is_enum< T >::value;


} /* namespace traits */

template< typename T, bool COND >
struct add_const_if
{
  using type = typename std::conditional< COND, typename std::add_const< T >::type, T >::type;
};

template< typename T, bool COND >
using add_const_if_t = typename add_const_if< T, COND >::type;

} /* namespace geosx */

#endif /* GEOSX_CODINGUTILITIES_TRAITS_HPP_ */
