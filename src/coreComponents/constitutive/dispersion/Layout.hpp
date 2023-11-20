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

#ifndef GEOS_CONSTITUTIVE_DISPERSION_LAYOUTS_HPP
#define GEOS_CONSTITUTIVE_DISPERSION_LAYOUTS_HPP

#include "common/DataTypes.hpp"
#include "common/GeosxConfig.hpp"

#include "LvArray/src/typeManipulation.hpp"
#include "RAJA/RAJA.hpp"

namespace geos
{
namespace constitutive
{
namespace dispersion
{

#if defined( GEOS_USE_DEVICE )
//TODO check
using LAYOUT_PHASE_VELOCITY = RAJA::PERM_JKIL;
#else
using LAYOUT_PHASE_VELOCITY = RAJA::PERM_IJKL;
#endif

/// Constitutive model phase velocity unit stride dimension
static constexpr int USD_PHASE_VELOCITY = LvArray::typeManipulation::getStrideOneDimension(
  LAYOUT_PHASE_VELOCITY{} );


}
}
}

#endif //GEOS_CONSTITUTIVE_DISPERSION_LAYOUTS_HPP
