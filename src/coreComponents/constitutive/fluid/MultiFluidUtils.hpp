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
template< typename T, int DIM, int USD, int USD_DC >
struct MultiFluidVarSlice
{
  internal::ArraySliceOrRef< T, DIM, USD > value;        /// variable value
  internal::ArraySliceOrRef< T, DIM, USD > dPres;        /// derivative w.r.t. pressure
  internal::ArraySliceOrRef< T, DIM, USD > dTemp;        /// derivative w.r.t. temperature
  internal::ArraySliceOrRef< T, DIM + 1, USD_DC > dComp; /// derivative w.r.t. composition
};

/**
 * @brief Struct holding views into fluid data, used to simplify parameter passing in kernel wrapper constructors.
 * @tparam NDIM number of dimensions
 * @tparam USD unit-stride-dim of primary property and derivatives
 * @tparam USD_DC unit-stride-dim of compositional derivatives
 */
template< typename T, int NDIM, int USD, int USD_DC >
struct MultiFluidVarView
{
  ArrayView< T, NDIM, USD > value;        ///< View into property values
  ArrayView< T, NDIM, USD > dPres;        ///< View into property pressure derivatives
  ArrayView< T, NDIM, USD > dTemp;        ///< View into property temperature derivatives
  ArrayView< T, NDIM + 1, USD_DC > dComp; ///< View into property compositional derivatives

  using SliceType = MultiFluidVarSlice< T, NDIM - 2, USD - 2, USD_DC - 2 >;

  GEOSX_HOST_DEVICE
  SliceType operator()( localIndex const k, localIndex const q ) const
  {
    return { value[k][q], dPres[k][q], dTemp[k][q], dComp[k][q] };
  }
};

/**
 * @brief Struct holding views into fluid data, used to simplify parameter passing in kernel wrapper constructors.
 * @tparam NDIM number of dimensions
 * @tparam PERM unit-stride-dim of primary property and derivatives
 * @tparam PERM_DC unit-stride-dim of compositional derivatives
 */
template< typename T, int NDIM, typename PERM, typename PERM_DC >
struct MultiFluidVar
{
  Array< real64, NDIM, PERM > value;        ///< Property values
  Array< real64, NDIM, PERM > dPres;        ///< Property pressure derivatives
  Array< real64, NDIM, PERM > dTemp;        ///< Property temperature derivatives
  Array< real64, NDIM + 1, PERM_DC > dComp; ///< Property compositional derivatives

  using ViewType = MultiFluidVarView< T, NDIM, getUSD< PERM >, getUSD< PERM_DC > >;
  using ViewTypeConst = MultiFluidVarView< T const, NDIM, getUSD< PERM >, getUSD< PERM_DC > >;

  using SliceType = typename ViewType::SliceType;
  using SliceTypeConst = typename ViewTypeConst::SliceType;

  ViewType toView()
  {
    return { value.toView(), dPres.toView(), dTemp.toView(), dComp.toView() };
  }

  ViewTypeConst toViewConst() const
  {
    return { value.toViewConst(), dPres.toViewConst(), dTemp.toViewConst(), dComp.toViewConst() };
  }
};

} // namespace constitutive

} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDUTILS_HPP_
