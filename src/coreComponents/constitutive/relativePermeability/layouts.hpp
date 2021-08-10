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
 * @file layouts.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_RELATIVEPERMEABILITY_LAYOUTS_HPP
#define GEOSX_CONSTITUTIVE_RELATIVEPERMEABILITY_LAYOUTS_HPP

#include "common/GeosxConfig.hpp"

#include "LvArray/src/typeManipulation.hpp"
#include "RAJA/RAJA.hpp"

namespace geosx
{
namespace constitutive
{
namespace relperm
{

#if defined( GEOSX_USE_CUDA )

/// Constitutive model phase relative permeability array layout
using LAYOUT_RELPERM = RAJA::PERM_JKI;
/// Constitutive model phase relative permeability saturation derivative array layout
using LAYOUT_RELPERM_DS = RAJA::PERM_JKLI;

#else

/// Constitutive model phase relative permeability array layout
using LAYOUT_RELPERM = RAJA::PERM_IJK;
/// Constitutive model phase relative permeability saturation derivative array layout
using LAYOUT_RELPERM_DS = RAJA::PERM_IJKL;

#endif

/// Constitutive model phase relative permeability unit stride dimension
static constexpr int USD_RELPERM = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_RELPERM{} );
/// Constitutive model phase relative permeability saturation derivative unit stride dimension
static constexpr int USD_RELPERM_DS = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_RELPERM_DS{} );

} // namespace relperm
} // namespace constitutive
} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_RELATIVEPERMEABILITY_LAYOUTS_HPP
