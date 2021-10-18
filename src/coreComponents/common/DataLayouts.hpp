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
 * @file DataLayouts.hpp
 */

#ifndef GEOSX_COMMON_DATALAYOUTS_HPP_
#define GEOSX_COMMON_DATALAYOUTS_HPP_

#include "common/GeosxConfig.hpp"

#include "LvArray/src/Array.hpp"
#include "RAJA/RAJA.hpp"

namespace geosx
{

/**
 * @brief Defaul permutation type for @p NDIM-dimensional array.
 * @tparam NDIM number of dimensions
 * @note useful when providing default template parameters
 */
template< int NDIM >
using defaultLayout = camp::make_idx_seq_t< NDIM >;

/**
 * @brief Just a handy shortcut for LvArray::typeManipulation::getStrideOneDimension.
 * @tparam PERM sequence type containing layout permutation
 * Value: corresponding unit stride dimension
 */
template< typename PERM >
static constexpr int getUSD = LvArray::typeManipulation::getStrideOneDimension( PERM {} );

namespace nodes
{

#if defined( GEOSX_USE_CUDA )

/// Node reference position permutation when using cuda.
using REFERENCE_POSITION_PERM = RAJA::PERM_JI;

/// Node total displacement permutation when using cuda.
using TOTAL_DISPLACEMENT_PERM = RAJA::PERM_JI;

/// Node incremental displacement permutation when using cuda.
using INCR_DISPLACEMENT_PERM = RAJA::PERM_JI;

/// Node velocity permutation when using cuda.
using VELOCITY_PERM = RAJA::PERM_JI;

/// Node acceleration permutation when using cuda.
using ACCELERATION_PERM = RAJA::PERM_JI;

#else

/// Node reference position permutation when not using cuda.
using REFERENCE_POSITION_PERM = RAJA::PERM_IJ;

/// Node total displacement permutation when not using cuda.
using TOTAL_DISPLACEMENT_PERM = RAJA::PERM_IJ;

/// Node incremental displacement permutation when not using cuda.
using INCR_DISPLACEMENT_PERM = RAJA::PERM_IJ;

/// Node velocity permutation when not using cuda.
using VELOCITY_PERM = RAJA::PERM_IJ;

/// Node acceleration permutation when not using cuda.
using ACCELERATION_PERM = RAJA::PERM_IJ;

#endif

/// Node reference position unit stride dimension.
static constexpr int REFERENCE_POSITION_USD = LvArray::typeManipulation::getStrideOneDimension( REFERENCE_POSITION_PERM {} );

/// Node total displacement unit stride dimension.
static constexpr int TOTAL_DISPLACEMENT_USD = LvArray::typeManipulation::getStrideOneDimension( TOTAL_DISPLACEMENT_PERM {} );

/// Node incremental displacement unit stride dimension.
static constexpr int INCR_DISPLACEMENT_USD = LvArray::typeManipulation::getStrideOneDimension( INCR_DISPLACEMENT_PERM {} );

/// Node velocity unit stride dimension.
static constexpr int VELOCITY_USD = LvArray::typeManipulation::getStrideOneDimension( VELOCITY_PERM {} );

/// Node acceleration unit stride dimension.
static constexpr int ACCELERATION_USD = LvArray::typeManipulation::getStrideOneDimension( ACCELERATION_PERM {} );

} // namespace nodes

namespace cells
{

#if defined( GEOSX_USE_CUDA )

/// Cell node map permutation when using cuda.
using NODE_MAP_PERMUTATION = RAJA::PERM_JI;

#else

/// Cell node map permutation when not using cuda.
using NODE_MAP_PERMUTATION = RAJA::PERM_IJ;

#endif

/// Cell node map unit stride dimension.
static constexpr int NODE_MAP_USD = LvArray::typeManipulation::getStrideOneDimension( NODE_MAP_PERMUTATION {} );

} // namespace cells

namespace solid
{

#if defined( GEOSX_USE_CUDA )

/// Constitutive model stress permutation when using cuda.
using STRESS_PERMUTATION = RAJA::PERM_KJI;

/// Constitutive model stiffness permutation when using cuda.
using STIFFNESS_PERMUTATION = RAJA::PERM_KJI;

#else

/// Constitutive model stress permutation when not using cuda.
using STRESS_PERMUTATION = RAJA::PERM_IJK;

/// Constitutive model stiffness permutation when not using cuda.
using STIFFNESS_PERMUTATION = RAJA::PERM_IJK;

#endif

/// Constitutive model stress unit stride dimension.
static constexpr int STRESS_USD = LvArray::typeManipulation::getStrideOneDimension( STRESS_PERMUTATION {} );

/// Constitutive model stiffness unit stride dimension.
static constexpr int STIFFNESS_USD = LvArray::typeManipulation::getStrideOneDimension( STIFFNESS_PERMUTATION {} );

} // namespace solid

namespace compflow
{

#if defined(GEOSX_USE_CUDA)

/// Component global density/fraction array layout
using LAYOUT_COMP = RAJA::PERM_JI;

/// Component global fraction compositional derivative array layout
using LAYOUT_COMP_DC = RAJA::PERM_JKI;

/// Phase property array layout
using LAYOUT_PHASE = RAJA::PERM_JI;

/// Phase property array layout
using LAYOUT_PHASE_DC = RAJA::PERM_JKI;

/// Phase composition array layout
using LAYOUT_PHASE_COMP = RAJA::PERM_JKI;

/// Fluid property compositional derivative array layout
using LAYOUT_FLUID_DC = RAJA::PERM_JI;

#else

/// Component global density/fraction array layout
using LAYOUT_COMP = RAJA::PERM_IJ;

/// Component global fraction compositional derivative array layout
using LAYOUT_COMP_DC = RAJA::PERM_IJK;

/// Phase property array layout
using LAYOUT_PHASE = RAJA::PERM_IJ;

/// Phase property array layout
using LAYOUT_PHASE_DC = RAJA::PERM_IJK;

/// Phase composition array layout
using LAYOUT_PHASE_COMP = RAJA::PERM_IJK;

/// Fluid property compositional derivative array layout
using LAYOUT_FLUID_DC = RAJA::PERM_IJ;

#endif

/// Component global density/fraction unit stride dimension
static constexpr int USD_COMP = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_COMP{} );

/// Component global fraction compositional derivative unit stride dimension
static constexpr int USD_COMP_DC = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_COMP_DC{} );

/// Phase property unit stride dimension
static constexpr int USD_PHASE = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_PHASE{} );

/// Phase property compositional derivative unit stride dimension
static constexpr int USD_PHASE_DC = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_PHASE_DC{} );

/// Phase property unit stride dimension
static constexpr int USD_PHASE_COMP = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_PHASE_COMP{} );

/// Fluid property compositional derivative unit stride dimension
static constexpr int USD_FLUID_DC = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_FLUID_DC{} );

} // namespace compflow

} // namespace geosx

#endif // GEOSX_COMMON_DATALAYOUTS_HPP_
