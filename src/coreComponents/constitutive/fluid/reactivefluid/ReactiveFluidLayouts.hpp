/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file Layouts.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_REACTIVEFLUID_REACTIVEFLUIDLAYOUTS_HPP
#define GEOS_CONSTITUTIVE_FLUID_REACTIVEFLUID_REACTIVEFLUIDLAYOUTS_HPP

#include "common/DataTypes.hpp"
#include "common/GeosxConfig.hpp"

#include "LvArray/src/typeManipulation.hpp"
#include "RAJA/RAJA.hpp"

namespace geos
{
namespace constitutive
{

namespace reactivefluid
{
struct DerivativeOffset
{
  /// index of derivative wrt pressure
  static integer constexpr dP = 0;
  /// index of derivative wrt temperature
  static integer constexpr dT = 1;
  /// index of first derivative wrt compositions
  static integer constexpr dC = 2;
};

/// indices of pressure, temperature, and composition derivatives
template< integer NC, integer IS_THERMAL >
struct DerivativeOffsetC {};

template< integer NC >
struct DerivativeOffsetC< NC, 1 >
{
  /// index of derivative wrt pressure
  static integer constexpr dP = 0;
  /// index of derivative wrt temperature
  static integer constexpr dT = dP + 1;
  /// index of first derivative wrt compositions
  static integer constexpr dC = dP+2;
  /// number of derivatives
  static integer constexpr nDer =  NC + 2;
};
template< integer NC >
struct DerivativeOffsetC< NC, 0 >
{
  /// index of derivative wrt pressure
  static integer constexpr dP = 0;
  /// index of first derivative wrt compositions
  static integer constexpr dC = dP+1;
  /// number of derivatives
  static integer constexpr nDer =  NC + 1;
};

  #if defined( GEOS_USE_DEVICE )

/// Constitutive model species variable array array layout
using LAYOUT_SPECIES = RAJA::PERM_JKI;
/// Constitutive model species derivative of species variable array layout
using LAYOUT_SPECIES_DC = RAJA::PERM_JKLI;

/// Constitutive model fluid property array layout
using LAYOUT_FLUID = RAJA::PERM_JI;
/// Constitutive model fluid property species derivative array layout
using LAYOUT_FLUID_DC = RAJA::PERM_JKI;

  #else

/// Constitutive model species variable array array layout
using LAYOUT_SPECIES = RAJA::PERM_IJK;
/// Constitutive model species derivative of species variable array layout
using LAYOUT_SPECIES_DC = RAJA::PERM_IJKL;

/// Constitutive model fluid property array layout
using LAYOUT_FLUID = RAJA::PERM_IJ;
/// Constitutive model fluid property species derivative array layout
using LAYOUT_FLUID_DC = RAJA::PERM_IJK;

  #endif


/// Constitutive model phase composition unit stride dimension
static constexpr int USD_SPECIES = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_SPECIES{} );
/// Constitutive model phase composition compositional derivative unit stride dimension
static constexpr int USD_SPECIES_DC = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_SPECIES_DC{} );

/// Constitutive model fluid property unit stride dimension
static constexpr int USD_FLUID = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_FLUID{} );
/// Constitutive model fluid property compositional derivative unit stride dimension
static constexpr int USD_FLUID_DC = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_FLUID_DC{} );

} // namespace reactivefluid
} // namespace constitutive
} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_REACTIVEFLUID_REACTIVEFLUIDLAYOUTS_HPP
