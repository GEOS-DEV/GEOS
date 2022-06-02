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
 * @file tensor_types_impl.hpp
 */

#ifndef GEOSX_TENSOR_TYPES_IMPL
#define GEOSX_TENSOR_TYPES_IMPL

#include "pow.hpp"
#include "prod.hpp"

namespace geosx
{

namespace tensor
{

/// 1D
/// Dynamic1dThreadLayout
constexpr int get_Dynamic1dThreadLayout_size(int MaxSize, int Rank)
{
   return Rank>2 ? pow(MaxSize,Rank-1) : 1;
}

/// Static1dThreadTensor
constexpr int get_Static1dThreadTensor_size(int Size0)
{
   GEOSX_UNUSED_VAR(Size0);
   return 1;
}

template <typename... Sizes>
constexpr int get_Static1dThreadTensor_size(int Size0, Sizes... sizes)
{
   GEOSX_UNUSED_VAR(Size0);
   return prod(sizes...);
}

/// 2D
/// Dynamic2dThreadLayout
constexpr int get_Dynamic2dThreadLayout_size(int MaxSize, int Rank)
{
   return Rank>2 ? pow(MaxSize,Rank-2) : 1;
}

/// Static2dThreadTensor
constexpr int get_Static2dThreadTensor_size(int Size0)
{
   GEOSX_UNUSED_VAR(Size0);
   return 1;
}

constexpr int get_Static2dThreadTensor_size(int Size0, int Size1)
{
   GEOSX_UNUSED_VAR(Size0, Size1);
   return 1;
}

template <typename... Sizes>
constexpr int get_Static2dThreadTensor_size(int Size0, int Size1, Sizes... sizes)
{
   GEOSX_UNUSED_VAR(Size0, Size1);
   return prod(sizes...);
}

/// 3D
/// Dynamic3dThreadLayout
constexpr int get_Dynamic3dThreadLayout_size(int MaxSize, int Rank)
{
   return Rank>3 ? pow(MaxSize,Rank-3) : 1;
}

/// Static3dThreadTensor
constexpr int get_Static3dThreadTensor_size(int Size0)
{
   GEOSX_UNUSED_VAR(Size0);
   return 1;
}

constexpr int get_Static3dThreadTensor_size(int Size0, int Size1)
{
   GEOSX_UNUSED_VAR(Size0, Size1);
   return 1;
}

constexpr int get_Static3dThreadTensor_size(int Size0, int Size1, int Size2)
{
   GEOSX_UNUSED_VAR(Size0, Size1, Size2);
   return 1;
}

template <typename... Sizes>
constexpr int get_Static3dThreadTensor_size(int Size0, int Size1, int Size2,
                                            Sizes... sizes)
{
   GEOSX_UNUSED_VAR(Size0, Size1, Size2);
   return prod(sizes...);
}

} // namespace tensor

} // namespace geosx

#endif // GEOSX_TENSOR_TYPES_IMPL
