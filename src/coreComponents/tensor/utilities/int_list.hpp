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
 * @file int_list.hpp
 */

#ifndef GEOSX_TENSOR_INT_LIST
#define GEOSX_TENSOR_INT_LIST

namespace geosx
{

namespace tensor
{

/// A struct to store a list of integers
template <int... Vs>
struct int_list { };

/// Append two int_list together into another int_list
template <typename, typename>
struct append_t;

template <int... VLs, int... VRs>
struct append_t<int_list<VLs...>,int_list<VRs...>>
{
   using type = int_list<VLs...,VRs...>;
};

template <typename L, typename R>
using append = typename append_t<L,R>::type;

/// Create an int_list containing Val N times
template <int Val, int N>
struct int_repeat_t
{
   using type = append<
                  int_list<Val>,
                  typename int_repeat_t<Val,N-1>::type
                >;
};

template <int Val>
struct int_repeat_t<Val,0>
{
   using type = int_list<>;
};

template <int Val, int N>
using int_repeat = typename int_repeat_t<Val,N>::type;

/// Instatiate TT with T, i.e. TT<T>
template <template<int...> class TT, typename T>
struct instantiate_t;

template <template<int...> class TT, int... Vals>
struct instantiate_t<TT, int_list<Vals...>>
{
   using type = TT<Vals...>;
};

template <template<int...> class TT, typename T>
using instantiate = typename instantiate_t<TT,T>::type;

/// Append an int value to a type list
template<int, typename>
struct append_to_type_seq_t { };

template<int V, int... Vs, template<int...> class TT>
struct append_to_type_seq_t<V, TT<Vs...>>
{
   using type = TT<Vs..., V>;
};

template<int Val, typename T>
using append_to_type_seq = typename append_to_type_seq_t<Val, T>::type;

/// Append the value V N times
template<int V, int N, template<int...> class TT>
struct repeat
{
   using type = append_to_type_seq<V, typename repeat<V, N-1, TT>::type>;
};

template<int V, template<int...> class TT>
struct repeat<V, 0, TT>
{
   using type = TT<>;
};

/// Append the value V1 N1 times, and then the value V2 N2 times.
template<int V1, int N1, int V2, int N2, template<int...> class TT>
struct rerepeat
{
   using type = append_to_type_seq<
      V1,
      typename rerepeat<V1, N1-1, V2, N2, TT>::type>;
};

template<int V1, int V2, int N2, template<int...> class TT>
struct rerepeat<V1, 0, V2, N2, TT>
{
   using type = append_to_type_seq<
      V2,
      typename rerepeat<V1, 0, V2, N2-1, TT>::type
      >;
};

template<int V1, int V2, template<int...> class TT>
struct rerepeat<V1, 0, V2, 0, TT>
{
   using type = TT<>;
};

/// Remove the Nth value in an int_list
template <int N, typename Tail, typename Head = int_list<>>
struct remove_t;

template <int N, int Val, int... TailVals, int... HeadVals>
struct remove_t<N, int_list<Val,TailVals...>, int_list<HeadVals...>>
{
   using type = typename remove_t<N-1,
                                  int_list<TailVals...>,
                                  int_list<HeadVals...,Val>>::type;
};

template <int Val, int... TailVals, int... HeadVals>
struct remove_t<0, int_list<Val,TailVals...>, int_list<HeadVals...>>
{
   using type = append< int_list<HeadVals...>, int_list<TailVals...> >;
};

template <int N, typename Vals>
using remove = typename remove_t<N, Vals>::type;

} // namespace tensor

} // geosx namespace

#endif // GEOSX_TENSOR_INT_LIST
