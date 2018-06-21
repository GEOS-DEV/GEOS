// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
/*
 * GeosxTraits.hpp
 *
 *  Created on: Feb 22, 2018
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_CODINGUTILITIES_GEOSXTRAITS_HPP_
#define SRC_COMPONENTS_CORE_SRC_CODINGUTILITIES_GEOSXTRAITS_HPP_
#include "common/DataTypes.hpp"
namespace geosx
{

namespace traits
{

template<typename>
struct is_string : std::false_type {};

template<>
struct is_string< string > : std::true_type {};




template<typename>
struct is_std_vector : std::false_type {};

template<typename T>
struct is_std_vector<std::vector<T> > : std::true_type {};


template<typename>
struct is_pair : std::false_type {};

template<typename T1, typename T2>
struct is_pair< std::pair<T1,T2> > : std::true_type{};


template<typename>
struct is_map : std::false_type {};

template<typename T_KEY, typename T_VAL>
struct is_map< map<T_KEY,T_VAL> > : std::true_type{};


template<typename>
struct is_set : std::false_type {};

template<typename T>
struct is_set< set<T> > : std::true_type{};


template<typename>
struct is_array : std::false_type {};

template< typename T, int NDIM, typename INDEX_TYPE >
struct is_array< multidimensionalArray::ManagedArray<T,NDIM,INDEX_TYPE> > : std::true_type{};



template<typename>
struct is_tensorT : std::false_type {};

template<>
struct is_tensorT< R1TensorT<3> > : std::true_type{};

template<>
struct is_tensorT< R2TensorT<3> > : std::true_type{};

template<>
struct is_tensorT< R2SymTensorT<3> > : std::true_type{};

}



namespace bufferOps
{


template< typename T >
struct is_noncontainer_type_packable
{
  static constexpr bool value = std::is_trivial<T>::value ||
                                std::is_arithmetic<T>::value ||
                                traits::is_tensorT<T>::value ||
                                traits::is_string<T>::value;
};
template< typename T >
constexpr bool is_noncontainer_type_packable<T>::value;


template<typename>
struct is_packable_array : std::false_type {};

template<typename T, int NDIM, typename INDEX_TYPE>
struct is_packable_array< multidimensionalArray::ManagedArray<T,NDIM,INDEX_TYPE> >
{
  static constexpr bool value = is_noncontainer_type_packable<T>::value;
};
template< typename T, int NDIM, typename INDEX_TYPE>
constexpr bool is_packable_array< multidimensionalArray::ManagedArray<T,NDIM,INDEX_TYPE> >::value;


template<typename>
struct is_packable_set : std::false_type {};

template< typename T >
struct is_packable_set< set<T> >
{
  static constexpr bool value = is_noncontainer_type_packable<T>::value;
};
template< typename T>
constexpr bool is_packable_set< set<T> >::value;






template<typename>
struct is_packable_map : std::false_type {};

template<typename T_KEY, typename T_VAL>
struct is_packable_map< map<T_KEY,T_VAL> >
{
  static constexpr bool value = is_noncontainer_type_packable<T_KEY>::value &&
                                ( is_noncontainer_type_packable<T_VAL>::value ||
                                  is_packable_array<T_VAL>::value );
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

}


}


#endif /* SRC_COMPONENTS_CORE_SRC_CODINGUTILITIES_GEOSXTRAITS_HPP_ */
