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

#include "LvArray/src/typeManipulation.hpp"
#include "RAJA/RAJA.hpp"

namespace geos
{
namespace constitutive
{
namespace multifluid
{

/// indices of pressure, temperature, and composition derivatives
struct DerivativeOffset
{
  /// index of derivative wrt pressure
  static integer constexpr dP = 0;
  /// index of derivative wrt temperature
  static integer constexpr dT = 1;
  /// index of first derivative wrt compositions
  static integer constexpr dC = 2;
};

#if defined( GEOS_USE_DEVICE )

/// Constitutive model phase property array layout
using LAYOUT_PHASE = RAJA::PERM_JKI;
/// Constitutive model phase property compositional derivative array layout
using LAYOUT_PHASE_DC = RAJA::PERM_JKLI;

/// Constitutive model phase composition array layout
using LAYOUT_PHASE_COMP = RAJA::PERM_JKLI;
/// Constitutive model phase composition compositional derivative array layout
using LAYOUT_PHASE_COMP_DC = RAJA::PERM_JKLMI;

/// Constitutive model fluid property array layout
using LAYOUT_FLUID = RAJA::PERM_JI;
/// Constitutive model fluid property compositional derivative array layout
using LAYOUT_FLUID_DC = RAJA::PERM_JKI;

#else

/// Constitutive model phase property array layout
using LAYOUT_PHASE = RAJA::PERM_IJK;
/// Constitutive model phase property compositional derivative array layout
using LAYOUT_PHASE_DC = RAJA::PERM_IJKL;

using LAYOUT_PHASE_VELOCITY = RAJA::PERM_IJKL;

/// Constitutive model phase composition array layout
using LAYOUT_PHASE_COMP = RAJA::PERM_IJKL;
/// Constitutive model phase composition compositional derivative array layout
using LAYOUT_PHASE_COMP_DC = RAJA::PERM_IJKLM;

/// Constitutive model fluid property array layout
using LAYOUT_FLUID = RAJA::PERM_IJ;
/// Constitutive model fluid property compositional derivative array layout
using LAYOUT_FLUID_DC = RAJA::PERM_IJK;

#endif

/// Constitutive model phase property unit stride dimension
static constexpr int USD_PHASE = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_PHASE{} );
/// Constitutive model phase property compositional derivative unit stride dimension
static constexpr int USD_PHASE_DC = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_PHASE_DC{} );

/// Constitutive model phase composition unit stride dimension
static constexpr int USD_PHASE_COMP = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_PHASE_COMP{} );
/// Constitutive model phase composition compositional derivative unit stride dimension
static constexpr int USD_PHASE_COMP_DC = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_PHASE_COMP_DC{} );
/// Constitutive model phase velocity unit stride dimension
static constexpr int USD_PHASE_VELOCITY = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_PHASE_VELOCITY{} );

/// Constitutive model fluid property unit stride dimension
static constexpr int USD_FLUID = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_FLUID{} );
/// Constitutive model fluid property compositional derivative unit stride dimension
static constexpr int USD_FLUID_DC = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_FLUID_DC{} );

} // namespace multifluid
} // namespace constitutive
} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_LAYOUTS_HPP
