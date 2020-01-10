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
 * @file traits.hpp
 */

#ifndef GEOSX_CODINGUTILITIES_TRAITS_HPP_
#define GEOSX_CODINGUTILITIES_TRAITS_HPP_

// Source includes
#include "common/DataTypes.hpp"
#include "cxx-utilities/src/templateHelpers.hpp"
#include "SFINAE_Macros.hpp"

// System includes
#include <type_traits>

namespace geosx
{

namespace traits
{

namespace internal
{
  HAS_ALIAS( value_type )

  HAS_ALIAS( pointer )

  template< class T,
            bool HASPOINTERTYPE = has_alias_pointer< T >::value >
  struct PointerHelper
  {
    using Pointer = T *;
    using ConstPointer = T const *;
  };

  template< class T >
  struct PointerHelper< T, true >
  {
    using Pointer = typename T::pointer;
    using ConstPointer = typename T::const_pointer;
  };

  template< typename T >
  struct has_data_method
  {
    HAS_MEMBER_FUNCTION_VARIANT( data, nonconst, typename PointerHelper< T >::Pointer, , , )
    HAS_MEMBER_FUNCTION_VARIANT( data, const,    typename PointerHelper< T >::Pointer, const, , )

    static constexpr bool value = has_memberfunction_vnonconst_data< T >::value ||
                                  has_memberfunction_vconst_data< T >::value;
  };

  template< typename T >
  struct has_chai_move_method
  {
    HAS_MEMBER_FUNCTION( move,
                         void,
                         ,
                         VA_LIST( chai::ExecutionSpace, bool ),
                         VA_LIST( chai::CPU, true ) )
    static constexpr bool value = has_memberfunction_move< T >::value;
  };

  template< typename T >
  struct has_empty_method
  {
    HAS_MEMBER_FUNCTION( empty, bool, const, , )
    static constexpr bool value = has_memberfunction_empty< T >::value;
  };

  template< typename T, typename INDEX_TYPE >
  struct has_size_method
  {
    HAS_MEMBER_FUNCTION( size, INDEX_TYPE, const, , )
    static constexpr bool value = has_memberfunction_size< T >::value;
  };

  template< typename T, typename INDEX_TYPE >
  struct has_dimension_size_method
  {
    HAS_MEMBER_FUNCTION( size, INDEX_TYPE, const, VA_LIST( int ), VA_LIST( 0 ) )
    static constexpr bool value = has_memberfunction_size< T >::value;
  };

  template< typename T, typename INDEX_TYPE >
  struct has_resize_method
  {
    HAS_MEMBER_FUNCTION( resize, void, , VA_LIST( INDEX_TYPE ), VA_LIST( INDEX_TYPE( 0 ) ) )
    static constexpr bool value = has_memberfunction_resize< T >::value;
  };

  template< typename T, typename DVT, typename INDEX_TYPE >
  struct has_resize_default_method
  {
    HAS_MEMBER_FUNCTION( resizeDefault, void, , VA_LIST( INDEX_TYPE, DVT const & ), VA_LIST( INDEX_TYPE( 0 ), std::declval< DVT const & >() ) )
    static constexpr bool value = has_memberfunction_resizeDefault< T >::value;
  };

  template< typename T >
  struct has_resize_dimensions_method
  {
    HAS_MEMBER_FUNCTION( resize, void, , VA_LIST( int, localIndex const * ),
                         VA_LIST( 0, static_cast< localIndex const * >( nullptr ) ) )
    static constexpr bool value = has_memberfunction_resize< T >::value;
  };

} // namespace internal

template< typename T >
using Pointer = typename internal::PointerHelper< T >::Pointer;

template< typename T >
using ConstPointer = typename internal::PointerHelper< T >::ConstPointer;

template< typename T >
constexpr bool has_alias_value_type = internal::has_alias_value_type< T >::value;

template< typename T >
constexpr bool has_data_method = internal::has_data_method< T >::value;

template< typename T >
constexpr bool has_chai_move_method = internal::has_chai_move_method< T >::value;

template< typename T >
constexpr bool has_empty_method = internal::has_empty_method< T >::value;

template< typename T >
constexpr bool has_size_method = internal::has_size_method< T, int >::value ||
                                 internal::has_size_method< T, unsigned int >::value ||
                                 internal::has_size_method< T, long >::value ||
                                 internal::has_size_method< T, unsigned long >::value ||
                                 internal::has_size_method< T, long long >::value ||
                                 internal::has_size_method< T, unsigned long long >::value;

template< typename T >
constexpr bool has_dimension_size_method = internal::has_dimension_size_method< T, int >::value ||
                                           internal::has_dimension_size_method< T, unsigned int >::value ||
                                           internal::has_dimension_size_method< T, long >::value ||
                                           internal::has_dimension_size_method< T, unsigned long >::value ||
                                           internal::has_dimension_size_method< T, long long >::value ||
                                           internal::has_dimension_size_method< T, unsigned long long >::value;

template< typename T >
constexpr bool has_resize_method = internal::has_resize_method< T, int >::value ||
                                   internal::has_resize_method< T, unsigned int >::value ||
                                   internal::has_resize_method< T, long >::value ||
                                   internal::has_resize_method< T, unsigned long >::value ||
                                   internal::has_resize_method< T, long long >::value ||
                                   internal::has_resize_method< T, unsigned long long >::value;

template< typename T, typename DVT >
constexpr bool has_resize_default_method = internal::has_resize_default_method< T, DVT, int >::value ||
                                           internal::has_resize_default_method< T, DVT, unsigned int >::value ||
                                           internal::has_resize_default_method< T, DVT, long >::value ||
                                           internal::has_resize_default_method< T, DVT, unsigned long >::value ||
                                           internal::has_resize_default_method< T, DVT, long long >::value ||
                                           internal::has_resize_default_method< T, DVT, unsigned long long >::value;

template< typename T >
constexpr bool has_resize_default_method< T, void > = false;

template< typename T >
constexpr bool has_resize_dimensions_method = internal::has_resize_dimensions_method< T >::value;

template< typename T >
constexpr bool is_string = is_base_of_v< std::string, T >;

template< typename T >
constexpr bool is_std_vector = is_instantiation_of< std::vector, T >;

template< typename T >
constexpr bool is_pair = is_instantiation_of< std::pair, T >;

template< typename T >
constexpr bool is_map = is_instantiation_of< mapBase, T >;

template< typename T >
constexpr bool is_set = is_instantiation_of< LvArray::SortedArray, T >;

template< typename T >
constexpr bool is_array = LvArray::isArray< T >;

template< typename T >
constexpr bool is_tensorT = is_same_v< std::remove_const_t< T >, R1Tensor > ||
                            is_same_v< std::remove_const_t< T >, R2Tensor > ||
                            is_same_v< std::remove_const_t< T >, R2SymTensor >;

} /* namespace traits */

template<typename T, bool COND>
struct add_const_if
{
  using type = typename std::conditional<COND, typename std::add_const<T>::type, T>::type;
};

template<typename T, bool COND>
using add_const_if_t = typename add_const_if<T, COND>::type;

} /* namespace geosx */

#endif /* GEOSX_CODINGUTILITIES_TRAITS_HPP_ */
