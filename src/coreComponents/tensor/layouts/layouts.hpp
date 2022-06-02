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
 * @file layouts.hpp
 */

#ifndef GEOSX_LAYOUTS_HPP_
#define GEOSX_LAYOUTS_HPP_

/**
 * Layouts are simple data structures (often empty) representing the mapping
 * from a rank N index to a linear index. The linear index representing the
 * tensor's associated index in the container for the value.
 * Their main purpose is to differentiate statically/dynamically known sizes
 * for the tensors, and to abstract threading models when utilizing tensors.
 * Layouts also defines the thread block sizes to run a kernel, and handle it
 * automatically for the user.
 * */

/// A dynamically sized layout
#include "dynamic_layout.hpp"
/// A dynamically sized layout where the first dimension is threaded
#include "dynamic_1dthread_layout.hpp"
/// A dynamically sized layout where the two first dimensions are threaded
#include "dynamic_2dthread_layout.hpp"
/// A dynamically sized layout where the three first dimensions are threaded
#include "dynamic_3dthread_layout.hpp"
/// A statically sized layout
#include "static_layout.hpp"
/// A statically sized layout, except for the last dimension (used for E-vectors)
#include "static_E_layout.hpp"
/// A statically sized layout where the first dimension is threaded
#include "static_1dthread_layout.hpp"
/// A statically sized layout where the two first dimensions are threaded
#include "static_2dthread_layout.hpp"
/// A statically sized layout where the three first dimensions are threaded
#include "static_3dthread_layout.hpp"
/// A layout that removes a chosen dimension to a given layout (used in Get)
#include "restricted_layout.hpp"

namespace mfem
{

// /// Strided Layout
// template <int Rank>
// class StridedLayout
// {
// private:
//    int strides[Rank];
//    int offsets[Rank];
//    int sizes[Rank];

// public:
//    template <typename... Idx> GEOSX_HOST_DEVICE inline
//    constexpr int index(Idx... idx) const
//    {
//       static_assert(sizeof...(Idx)==Rank,"Wrong number of argumets.");
//       return StridedIndex<1>::eval(offsets, strides, idx...);
//    }

//    // Can be constexpr if Tensor inherit from Layout
//    template <int N> GEOSX_HOST_DEVICE inline
//    int Size() const
//    {
//       static_assert(N>=0 && N<Rank,"Accessed size is higher than the rank of the Tensor.");
//       return sizes[N];
//    }

// private:
//    template <int N>
//    struct StridedIndex
//    {
//       template <typename... Idx>
//       static inline int eval(int* offsets, int* strides, int first, Idx... args)
//       {
//          return (offsets[N-1]+first)*strides[N-1] + StridedIndex<N+1>::eval(args...);
//       }
//    };

//    template <>
//    struct StridedIndex<Rank>
//    {
//       template <typename... Idx>
//       static inline int eval(int* offsets, int* strides, int first)
//       {
//          return (offsets[Rank-1]+first)*strides[Rank-1];
//       }
//    };
// };

} // namespace mfem

#endif // GEOSX_LAYOUTS_HPP_
