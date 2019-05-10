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

/*
 * GeosxTraits.hpp
 *
 *  Created on: Feb 22, 2018
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_CODINGUTILITIES_GEOSXTRAITS_HPP_
#define SRC_COMPONENTS_CORE_SRC_CODINGUTILITIES_GEOSXTRAITS_HPP_

#include <type_traits>
#include "common/DataTypes.hpp"
#include "templateHelpers.hpp"

namespace geosx
{

namespace traits
{

template <class T>
constexpr bool is_string = is_instance_of_v<std::string, T>;

template <class T>
constexpr bool is_std_vector = is_instantiation_of_v<std::vector, T>;

template <class T>
constexpr bool is_pair = is_instantiation_of_v<std::pair, T>;

template <class T>
constexpr bool is_map = is_instantiation_of_v<std::map, T> || is_instantiation_of_v<std::unordered_map, T>;

template <class T>
constexpr bool is_set = is_instantiation_of_v<LvArray::SortedArray, T>;

template<typename>
constexpr bool is_array = false;

template< typename T, int NDIM, typename INDEX_TYPE >
constexpr bool is_array< LvArray::Array<T,NDIM,INDEX_TYPE> > = true;

template <class T>
constexpr bool is_tensorT = is_instance_of_v<R1Tensor, T> ||
                            is_instance_of_v<R2Tensor, T> ||
                            is_instance_of_v<R2SymTensor, T>;

} /* namespace traits */



namespace bufferOps
{

/* Forward declaration of is_packable */
template< typename T >
struct is_packable;


template< typename T >
struct is_noncontainer_type_packable
{
  static constexpr bool value = std::is_trivial<T>::value ||
                                std::is_arithmetic<T>::value ||
                                traits::is_tensorT<T> ||
                                traits::is_string<T>;
};
template< typename T >
constexpr bool is_noncontainer_type_packable<T>::value;

template<typename>
struct is_packable_array : std::false_type {};

template<typename T, int NDIM, typename INDEX_TYPE>
struct is_packable_array< LvArray::Array<T,NDIM,INDEX_TYPE> > : is_packable<T> {};

template<typename T, int NDIM, typename INDEX_TYPE>
struct is_packable_array< LvArray::ArrayView<T,NDIM,INDEX_TYPE> > : is_packable<T> {};

template<typename T, int NDIM, typename INDEX_TYPE>
struct is_packable_array< LvArray::ArraySlice<T,NDIM,INDEX_TYPE> > : is_packable<T> {};


template<typename>
struct is_packable_set : std::false_type {};

template< typename T >
struct is_packable_set< set<T> >
{
  static constexpr bool value = is_packable<T>::value;
};
template< typename T>
constexpr bool is_packable_set< set<T> >::value;


template<typename>
struct is_packable_map : std::false_type {};

template<typename T_KEY, typename T_VAL>
struct is_packable_map< map<T_KEY,T_VAL> >
{
  static constexpr bool value = is_packable<T_KEY>::value &&
                                is_packable<T_VAL>::value;
};
template< typename T_KEY, typename T_VAL>
constexpr bool is_packable_map< map<T_KEY,T_VAL> >::value;


template< typename T >
struct is_packable
{
  static constexpr bool value = is_noncontainer_type_packable<T>::value ||
                                is_packable_array<T>::value ||
                                is_packable_map<T>::value ||
                                is_packable_set<T>::value ;

};
template< typename T >
constexpr bool is_packable<T>::value;


template< typename T >
struct is_packable_by_index
{
  static constexpr bool value = is_packable_array<T>::value  ;

};
template< typename T >
constexpr bool is_packable_by_index<T>::value;

} /* namespace bufferOps */

template<typename T, bool COND>
struct add_const_if
{
  using type = typename std::conditional<COND, typename std::add_const<T>::type, T>::type;
};

template<typename T, bool COND>
using add_const_if_t = typename add_const_if<T, COND>::type;

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CODINGUTILITIES_GEOSXTRAITS_HPP_ */
