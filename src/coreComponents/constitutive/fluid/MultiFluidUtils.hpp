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
 * @file MultiFluidUtils.hpp
 */
#ifndef GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDUTILS_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDUTILS_HPP_

#include "common/DataTypes.hpp"

namespace geosx
{

namespace constitutive
{

namespace internal
{

template< typename T, int DIM, int USD >
struct ArraySliceOrRefHelper
{
  using type = ArraySlice< T, DIM, USD >;
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

/**
 * @brief Helper struct used to represent a variable and its compositional derivatives
 * @tparam DIM number of dimensions
 */
template< int DIM, int USD, int USD_DC >
struct CompositionalVarContainer
{
  internal::ArraySliceOrRef< real64, DIM, USD > const & value; // variable value
  internal::ArraySliceOrRef< real64, DIM, USD > const & dPres; // derivative w.r.t. pressure
  internal::ArraySliceOrRef< real64, DIM, USD > const & dTemp; // derivative w.r.t. temperature
  internal::ArraySliceOrRef< real64, DIM + 1, USD_DC > const & dComp; // derivative w.r.t. composition
};

} // namespace constitutive

} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDUTILS_HPP_
