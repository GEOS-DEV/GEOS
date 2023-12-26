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
 * @file Layouts.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_LAYOUTS_HPP
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_LAYOUTS_HPP

#include "common/DataTypes.hpp"
#include "common/GeosxConfig.hpp"

#include "RAJA/RAJA.hpp"

namespace geos
{
namespace constitutive
{
namespace multifluid
{

/// Indices of pressure, temperature, and composition derivatives
struct DerivativeOffset
{
  /// index of derivative wrt pressure
  static integer constexpr dP = 0;
  /// index of derivative wrt temperature
  static integer constexpr dT = 1;
  /// index of first derivative wrt compositions
  static integer constexpr dC = 2;
};

// Data layout depends on the number of dimensions of the data
template< int NDIM >
struct Layout
{
  using PERM = RAJA::PERM_I;
};

#if defined( GEOS_USE_DEVICE )

template<> struct Layout< 2 >
{
  using PERM = RAJA::PERM_JI;
};
template<> struct Layout< 3 >
{
  using PERM = RAJA::PERM_JKI;
};
template<> struct Layout< 4 >
{
  using PERM = RAJA::PERM_JKLI;
};
template<> struct Layout< 5 >
{
  using PERM = RAJA::PERM_JKLMI;
};

#else

template<> struct Layout< 2 >
{
  using PERM = RAJA::PERM_IJ;
};
template<> struct Layout< 3 >
{
  using PERM = RAJA::PERM_IJK;
};
template<> struct Layout< 4 >
{
  using PERM = RAJA::PERM_IJKL;
};
template<> struct Layout< 5 >
{
  using PERM = RAJA::PERM_IJKLM;
};

#endif

namespace internal
{

template< typename T, int DIM, int USD >
struct ArraySliceOrRefHelper
{
  using type = geos::ArraySlice< T, DIM, USD >;
};

// an array slice of DIM=0 decays to a reference to scalar
template< typename T, int USD >
struct ArraySliceOrRefHelper< T, 0, USD >
{
  using type = T &;
};

template< typename T, int DIM, int USD=DIM-1 >
using ArraySliceOrRef = typename ArraySliceOrRefHelper< T, DIM, USD >::type;

} // namespace internal

template< typename T, int NDIM >
using Array = geos::Array< T, NDIM, typename Layout< NDIM >::PERM >;

template< typename T, int NDIM, localIndex CAPACITY >
using StackArray = geos::StackArray< T, NDIM, CAPACITY, typename Layout< NDIM >::PERM >;

template< typename T, int NDIM >
using ArrayView = geos::ArrayView< T, NDIM, getUSD< typename Layout< NDIM >::PERM > >;

template< typename T, int NDIM >
using ArraySlice = internal::ArraySliceOrRef< T, NDIM, getUSD< typename Layout< NDIM >::PERM > >;

} // namespace multifluid
} // namespace constitutive
} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_LAYOUTS_HPP
