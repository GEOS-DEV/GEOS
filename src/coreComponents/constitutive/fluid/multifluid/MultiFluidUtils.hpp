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
#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_MULTIFLUIDUTILS_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_MULTIFLUIDUTILS_HPP_

#include "common/DataTypes.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @brief Helper struct used to represent a variable and its compositional derivatives
 * @tparam DIM number of dimensions
 */
template< typename T, int DIM >
struct MultiFluidVarSlice
{
  GEOS_HOST_DEVICE
  MultiFluidVarSlice( multifluid::ArraySlice< T, DIM > inputValue,
                      multifluid::ArraySlice< T, DIM+1 > inputDerivs ):
    value( inputValue ),
    derivs( inputDerivs )
  {}

  multifluid::ArraySlice< T, DIM > value;        /// variable value
  multifluid::ArraySlice< T, DIM + 1 > derivs; /// derivative w.r.t. pressure, temperature, compositions
};

/**
 * @brief Struct holding views into fluid data, used to simplify parameter passing in kernel wrapper constructors.
 * @tparam NDIM number of dimensions
 */
template< typename T, int NDIM >
struct MultiFluidVarView
{
  MultiFluidVarView() = default;

  GEOS_HOST_DEVICE
  MultiFluidVarView ( MultiFluidVarView const & src ):
    value( src.value ),
    derivs( src.derivs )
  {}

  GEOS_HOST_DEVICE
  MultiFluidVarView ( multifluid::ArrayView< T, NDIM > const & valueSrc,
                      multifluid::ArrayView< T, NDIM + 1 > const & derivsSrc ):
    value( valueSrc ),
    derivs( derivsSrc )
  {};

  multifluid::ArrayView< T, NDIM > value;       ///< View into property values
  multifluid::ArrayView< T, NDIM + 1 > derivs;  ///< View into property derivatives w.r.t. pressure, temperature, compositions

  using SliceType = MultiFluidVarSlice< T, NDIM - 2 >;

  GEOS_HOST_DEVICE
  SliceType operator()( localIndex const k, localIndex const q ) const
  {
    return { value[k][q], derivs[k][q] };
  }
};

/**
 * @brief Struct holding views into fluid data, used to simplify parameter passing in kernel wrapper constructors.
 * @tparam NDIM number of dimensions
 */
template< typename T, int NDIM >
struct MultiFluidVar
{
  multifluid::Array< real64, NDIM > value;      ///< Property values
  multifluid::Array< real64, NDIM + 1 > derivs; ///< Property derivatives w.r.t. pressure, temperature, compositions

  using ViewType = MultiFluidVarView< T, NDIM >;
  using ViewTypeConst = MultiFluidVarView< T const, NDIM >;

  using SliceType = typename ViewType::SliceType;
  using SliceTypeConst = typename ViewTypeConst::SliceType;

  ViewType toView()
  {
    return { value.toView(), derivs.toView() };
  }

  ViewTypeConst toViewConst() const
  {
    return { value.toViewConst(), derivs.toViewConst() };
  }
};

namespace detail
{
/**
 * @brief Utility function to convert mass fractions to mole fractions
 * @tparam ARRAY1 the type of array storing the component molar weights
 * @tparam ARRAY2 the type of array storing the component mass fractions
 * @tparam ARRAY3 the type of array storing the component mole fractions
 * @param[in] componentMolarWeight the component molar weights
 * @param[in] massFractions the component mass fractions
 * @param[out] moleFractions the newly converted component mole fractions
 */
template< typename ARRAY1, typename ARRAY2, typename ARRAY3 >
GEOS_HOST_DEVICE
void convertToMoleFractions( integer numComponents,
                             ARRAY1 const & componentMolarWeight,
                             ARRAY2 const & massFractions,
                             ARRAY3 & moleFractions )
{
  real64 totalMolality = 0.0;
  for( integer ic = 0; ic < numComponents; ++ic )
  {
    moleFractions[ic] = massFractions[ic] / componentMolarWeight[ic];
    totalMolality += moleFractions[ic];
  }

  constexpr real64 epsilon = LvArray::NumericLimits< real64 >::epsilon;

  real64 const totalMolalityInv = epsilon < totalMolality ? 1.0 / totalMolality : 0.0;
  for( integer ic = 0; ic < numComponents; ++ic )
  {
    moleFractions[ic] *= totalMolalityInv;
  }
}

} // namespace detail

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_MULTIFLUIDUTILS_HPP_
