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

#ifndef GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_LAYOUTS_HPP
#define GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_LAYOUTS_HPP

#include "common/DataTypes.hpp"
#include "common/GeosxConfig.hpp"

#include "LvArray/src/typeManipulation.hpp"
#include "RAJA/RAJA.hpp"

namespace geos
{
namespace constitutive
{
namespace singlefluid
{

struct DerivativeOffset
{
  /// index of derivative wrt pressure
  static integer constexpr dP = 0;
  /// index of derivative wrt temperature
  static integer constexpr dT = 1;
};

}
}
}

#endif //GEOS_CONSTITUTIVE_FLUID_LAYOUTS_HPP
